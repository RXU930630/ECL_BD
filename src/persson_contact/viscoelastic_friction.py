"""Isothermal Persson rigid-rough / smooth-viscoelastic, uniform nominal pressure.

Reference: Persson, Rubber friction: Theory, mechanisms, and challenges,
JCP 163, 141001 (2025), equations 17--20, DOI 10.1063/5.0293616.
No flash heating, nonlinear softening, adhesion force, fluid or macro rolling loss.
"""
import numpy as np
from scipy.special import erf
from numpy.polynomial.legendre import leggauss


class PerssonViscoelasticFriction:
    def __init__(self,psd,material,*,q_intervals=512,angular_points=64,kappa=.5):
        if type(q_intervals) is not int or not 32<=q_intervals<=20000:raise ValueError('q_intervals must be 32..20000')
        if type(angular_points) is not int or not 8<=angular_points<=1024:raise ValueError('angular_points must be 8..1024')
        if isinstance(kappa,bool) or not np.isfinite(kappa) or not 0<=kappa<=1:raise ValueError('kappa must lie in [0,1]')
        self.psd=psd;self.material=material;self.kappa=kappa
        self.q_intervals=q_intervals;self.angular_points=angular_points
        knots=10**np.asarray(psd.log10_q)
        self.edges=np.unique(np.r_[np.geomspace(psd.q0,psd.q1,q_intervals+1),knots[(knots>psd.q0)&(knots<psd.q1)]])
        self.mid=np.sqrt(self.edges[:-1]*self.edges[1:])
        self.weights=psd.moments(3,self.edges[:-1],self.edges[1:])
        self.left_weights=psd.moments(3,self.edges[:-1],self.mid)
        x,w=leggauss(angular_points);theta=(x+1)*np.pi/4
        self.cos=np.cos(theta);self.angular_weights=w*np.pi/4
        self._speed_cache={}

    def _integrals(self,speed):
        v=abs(float(speed))
        if not np.isfinite(v) or v<=0:raise ValueError('A positive sliding-speed magnitude is required; static friction is not computed')
        if v not in self._speed_cache:
            ep,epp=self.material.storage_loss(self.mid[:,None]*v*self.cos[None,:])
            denominator=1-self.material.nu**2
            # Four equal quadrants; loss changes sign at negative frequency.
            angular_square=4*((ep**2+epp**2)@self.angular_weights)/denominator**2
            angular_loss=4*(epp@(self.angular_weights*self.cos))/denominator
            dG=self.weights*angular_square/8
            Gmid=np.r_[0.,np.cumsum(dG)[:-1]]+self.left_weights*angular_square/8
            self._speed_cache[v]=(Gmid,float(np.sum(dG)),angular_loss)
        return self._speed_cache[v]

    def evaluate(self,pressure_pa,speed_m_s,shear_strength_pa=0.):
        p=float(pressure_pa);tau=float(shear_strength_pa)
        if not np.isfinite(p) or p<=0 or not np.isfinite(tau) or tau<0:raise ValueError('pressure>0 and shear strength>=0 required')
        Gmid,Gend,loss=self._integrals(speed_m_s)
        area=erf(p/(2*np.sqrt(Gmid)));correction=self.kappa+(1-self.kappa)*area**2
        mu=float(np.sum(self.weights*area*correction*loss)/(2*p))
        a=float(erf(p/(2*np.sqrt(Gend))))
        adhesive=tau*a/p
        return {'pressure_pa':p,'speed_m_s':float(speed_m_s),'relative_contact_area':a,
            'mu_viscoelastic':mu,'mu_adhesive':adhesive,'mu_fluid':0.,'mu_total':mu+adhesive}
