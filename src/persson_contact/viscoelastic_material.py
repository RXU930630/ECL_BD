"""Persson-format master curves; Hz/MPa input, explicit SI and temperature shifts."""
from dataclasses import dataclass
import hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator


def read_three_columns(path):
    path=Path(path);rows=[];header='';seen_header=False
    for number,line in enumerate(path.read_text(encoding='utf-8-sig').splitlines(),1):
        text=line.strip()
        if not text or text.startswith(('#','!')):continue
        try:
            row=[float(v.replace('D','E').replace('d','e')) for v in text.split()]
        except ValueError as exc:
            if not rows and not seen_header:
                header=text;seen_header=True;continue
            raise ValueError(f'{path.name}:{number}: expected three numeric columns') from exc
        if len(row)!=3 or not np.all(np.isfinite(row)):
            raise ValueError(f'{path.name}:{number}: expected three finite numeric columns')
        rows.append(row)
    data=np.asarray(rows,float)
    if data.ndim!=2 or len(data)<2:raise ValueError(f'{path.name}: at least two data rows required')
    if np.any(np.diff(data[:,0])<=0):raise ValueError(f'{path.name}: first column must strictly increase')
    return header,data


class MasterCurve:
    """Returns positive storage/loss magnitudes. Does not assert Kramers-Kronig consistency.

    Complex API uses exp(+i omega t): E=E'+i E'' for omega>0.
    Zero frequency is not inferred from finite-band DMA measurements.
    """
    def __init__(self,path,poisson_ratio,*,reference_temperature_c,temperature_c,
                 shift_factors_file=None,shift_factors=None,modulus_type='young',frequency_policy='error'):
        if isinstance(poisson_ratio,bool) or not np.isfinite(poisson_ratio) or not -1<float(poisson_ratio)<=.5:
            raise ValueError('Invalid poisson_ratio')
        self.nu=float(poisson_ratio);self.path=Path(path)
        self.name,data=read_three_columns(path)
        if np.any(data[:,0]<=0) or np.any(data[:,1]<=0) or np.any(data[:,2]<0):
            raise ValueError('Master curve requires f>0, storage>0, loss>=0')
        if modulus_type not in ('young','shear','reduced'):raise ValueError('modulus_type must be young/shear/reduced')
        if frequency_policy not in ('error','constant_endpoint'):raise ValueError('Invalid frequency_policy')
        for value in (reference_temperature_c,temperature_c):
            if isinstance(value,bool) or not np.isfinite(value) or value<=-273.15:raise ValueError('Invalid temperature in Celsius')
        self.reference_temperature_c=float(reference_temperature_c);self.temperature_c=float(temperature_c)
        if shift_factors_file is not None and shift_factors is not None:raise ValueError('Choose shift table OR manual shift factors')
        self.aT=self.bT=1.;shift_sha=None
        if shift_factors_file is not None:
            _,table=read_three_columns(shift_factors_file)
            if np.any(table[:,1:]<=0):raise ValueError('Shift factors must be positive linear values, not logarithms')
            if not table[0,0]<=temperature_c<=table[-1,0] or not table[0,0]<=reference_temperature_c<=table[-1,0]:
                raise ValueError('Temperature outside shift table; extrapolation is disabled')
            ref=np.exp([np.interp(reference_temperature_c,table[:,0],np.log(table[:,i])) for i in (1,2)])
            if not np.allclose(ref,1.,rtol=.005,atol=0):raise ValueError('Shift table reference is inconsistent: aT=bT=1 required at reference temperature')
            self.aT,self.bT=np.exp([np.interp(temperature_c,table[:,0],np.log(table[:,i])) for i in (1,2)])
            shift_sha=hashlib.sha256(Path(shift_factors_file).read_bytes()).hexdigest()
        elif shift_factors is not None:
            if not isinstance(shift_factors,dict) or set(shift_factors)!={'aT','bT'}:raise ValueError('Provide manual aT and bT')
            self.aT=float(shift_factors['aT']);self.bT=float(shift_factors['bT'])
            if not np.isfinite(self.aT*self.bT) or min(self.aT,self.bT)<=0:raise ValueError('aT and bT must be positive finite')
            if temperature_c==reference_temperature_c and not np.allclose([self.aT,self.bT],1.,rtol=1e-8):raise ValueError('Reference temperature requires unit shifts')
        elif temperature_c!=reference_temperature_c:
            raise ValueError('Different temperature requires a shift table or explicit aT,bT')
        conversion={'young':1.,'shear':2*(1+self.nu),'reduced':1-self.nu**2}[modulus_type]
        self.frequency_hz=data[:,0].copy();self.storage_pa=data[:,1]*1e6*conversion;self.loss_pa=data[:,2]*1e6*conversion
        x=np.log(self.frequency_hz)
        self._storage=PchipInterpolator(x,np.log(self.storage_pa),extrapolate=False)
        self._loss_log=bool(np.all(self.loss_pa>0))
        self._loss=PchipInterpolator(x,np.log(self.loss_pa) if self._loss_log else self.loss_pa,extrapolate=False)
        self.frequency_policy=frequency_policy;self.query_min_hz=np.inf;self.query_max_hz=0.;self.out_of_band=False
        self.metadata={'material_name':self.name,'master_curve_sha256':hashlib.sha256(self.path.read_bytes()).hexdigest(),
            'master_curve_file':str(self.path),'shift_table_sha256':shift_sha,'modulus_type_input':modulus_type,
            'input_units':'Hz, MPa, MPa; loss is a nonnegative magnitude','internal_modulus':'Young modulus, Pa',
            'reference_temperature_c':self.reference_temperature_c,'temperature_c':self.temperature_c,
            'aT':float(self.aT),'bT':float(self.bT),'frequency_policy':frequency_policy,
            'shift_convention':'E(omega,T)=bT E(aT omega,Tref)',
            'frequency_range_hz':[float(data[0,0]),float(data[-1,0])],
            'complex_convention':'exp(+i omega t); E(-omega)=conjugate(E(omega))',
            'causality_validation':'Not established by interpolation of measured storage/loss alone'}
        for a in (self.frequency_hz,self.storage_pa,self.loss_pa):a.flags.writeable=False

    def storage_loss(self,omega_rad_s):
        omega=np.asarray(omega_rad_s,float)
        if np.any(~np.isfinite(omega)) or np.any(omega==0):raise ValueError('Finite nonzero omega required; zero-frequency modulus is not inferred')
        f=self.aT*abs(omega)/(2*np.pi)
        self.query_min_hz=min(self.query_min_hz,float(f.min()));self.query_max_hz=max(self.query_max_hz,float(f.max()))
        outside=(f<self.frequency_hz[0])|(f>self.frequency_hz[-1]);self.out_of_band|=bool(np.any(outside))
        if np.any(outside) and self.frequency_policy=='error':
            raise ValueError(f'Master curve frequency coverage insufficient: requested {f.min():.6g}..{f.max():.6g} Hz after aT shift')
        f=np.clip(f,self.frequency_hz[0],self.frequency_hz[-1]);x=np.log(f)
        storage=np.exp(self._storage(x))*self.bT
        loss=self._loss(x);loss=np.exp(loss) if self._loss_log else loss
        return storage,np.maximum(loss,0)*self.bT

    def complex_modulus(self,omega_rad_s):
        a,b=self.storage_loss(omega_rad_s)
        return a+1j*np.sign(omega_rad_s)*b

    @property
    def diagnostics(self):
        return {**self.metadata,'queried_master_frequency_range_hz':None if not np.isfinite(self.query_min_hz) else [self.query_min_hz,self.query_max_hz],
                'out_of_band':self.out_of_band}
