program reference_5000
implicit real*8 (a-h,o-z)
real*8 pc(5000),sigmaeff(5000),sigmaeff1(5000),pcontx(20000),pcontx0(20000),xintegral(20000)
open(10,file='input.dat')
read(10,*) Nx,Npa,Np,xL,p00,RMShc,barv,ETemp,patm,paMax,old
if(Nx.gt.20000.or.Np.gt.5000) stop 2
delx=xL/(Nx+1d0)
delpa=paMax/Npa
do i=1,Np
read(10,*) pc(i),sigmaeff(i)
enddo
do ix=1,Nx
read(10,*) pcontx0(ix)
pcontx(ix)=pcontx0(ix)
enddo
close(10)
do ipa=1,Npa
pa=ipa*delpa
      s=0.0
      do 2 ix=1,Nx
      pcont1=pcontx(ix)/p00
      call splint(pc,sigmaeff,sigmaeff1,Np,pcont1,seff1) 
      s=s+1.0/seff1
      xintegral(ix)=s 
  2   continue
      do 3 ix=1,Nx
      pcontxx=pcontx0(ix)-pa+pa*xintegral(ix)/xintegral(Nx)
      pcontx(ix)=old*pcontx(ix)+(1.0-old)*pcontxx
      pfluidxx=pcontx0(ix)-pcontx(ix)
      if (pcontx(ix).lt.0.0) pfluidxx=pa
      if (pcontx(ix).lt.0.0) pcontx(ix)=0.0
   3 continue
enddo
      dotQ=pa*xL*barv*(RMShc**2.0)/(2.0*ETemp*delx*xintegral(Nx))
      xV=dotQ*Etemp/patm
open(20,file='output.dat')
write(20,'(2ES26.17E3)') dotQ,xV
do ix=1,Nx
pf=pcontx0(ix)-pcontx(ix)
if(pcontx(ix).le.0d0) pf=pa
write(20,'(3ES26.17E3)') ix*delx,pcontx(ix),pf
enddo
close(20)
end program
subroutine splint(xa,ya,y2,n,x,y)
implicit real*8 (a-h,o-z)
real*8 xa(n),ya(n),y2(n)
klo=1
khi=n
do while(khi-klo.gt.1)
k=(khi+klo)/2
if(xa(k).gt.x) then
khi=k
else
klo=k
endif
enddo
h=xa(khi)-xa(klo)
y=ya(klo)+(x-xa(klo))*(ya(khi)-ya(klo))/h
return
end
