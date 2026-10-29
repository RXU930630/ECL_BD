program pressure_flow_reference
implicit real*8 (a-h,o-z)
real*8 Probabilityrms(100000),vecSeparationu1(1),vecAria(1)
open(10,file='input.dat')
read(10,*) nurange1,averagesep,relArea,percolation,uend,NittSeal,feedback1,Pumin,Palpha
if(nurange1.gt.100000) stop 2
xnDstart=2d0
xnDstop=1d0/(1d0-percolation)
Aratioperc=relArea/percolation
xnDrun=xnDstop
if(Aratioperc.lt.1d0) xnDrun=xnDstart*(1d0-Aratioperc)+xnDstop*Aratioperc
ellipsgamma=1d0
phip0=1d0
p0=0d0
Deltau=1d0
istop=1
! The extracted loop reads vecSeparationu1(istop); declaration below is scalar-array bridge.
vecSeparationu1(1)=uend
vecAria(1)=relArea
flowsigmaeff1=averagesep**3
flowsigmaeff1y=flowsigmaeff1
AreaContact=1d0
 do jj=1,nurange1
 read(10,*) dummy,Probabilityrms(jj)
 AreaContact=AreaContact-Probabilityrms(jj)
 enddo
close(10)
open(3718,file='iteration.log')
      do 2277 itt=1,NittSeal
!      do 2277 itt=1,200000
      ellipsgammac=ellipsgamma*sqrt((flowsigmaeff1y/flowsigmaeff1))
!      ellipsgamma1=ellipsgammac+xnDrun-2.0
      ellipsgamma1=ellipsgammac*(xnDrun-1.0)
      ellipsgamma11=ellipsgammac/(xnDrun-1.0)
!      if (itt.eq.3) write(5227,*) ellipsgamma11,ellipsgamma1 
!      if (itt.eq.3) write(5227,*) flowsigmaeffy, flowsigmaeff
!      if (itt.eq.3) write(5227,*) 
      


      
      sa=0.0
      say=0.0
      sb=0.0
      uend=vecSeparationu1(istop)
      flowsigma1end=uend**3.0

      flowsigma1end=flowsigma1end*phip0




!      do 3377 jj=1,nurange*4
      do 3377 jj=1,nurange1
!      u=xMinvecSeparationu1+jj*Deltau
      u=Pumin*exp(Palpha*jj)
      Pr=Probabilityrms(jj) 

      flowsigma1=u**3.0-flowsigma1end

!      bond percolation:
!      xa=flowsigma1+2.0*flowsigmaeff1

!      2D Bruggeman:  
      xa=flowsigma1+ellipsgamma1*flowsigmaeff1
      xay=ellipsgamma11*flowsigma1+flowsigmaeff1y
      sa=sa+((ellipsgamma1+1.0)/xa)*Pr
      say=say+((ellipsgamma11+1.0)/xay)*Pr
      sb=sb+Pr
 3377 continue
      sa=sa*Deltau+AreaContact*(ellipsgamma1+1.0)/(ellipsgamma1*flowsigmaeff1) 
      say=say*Deltau+AreaContact*(ellipsgamma11+1.0)/(flowsigmaeff1y) 
      sb=sb*Deltau 
!      2D Bruggeman:  

      sa=1.0/sa      
      say=1.0/say      

      flowsigmaeff1=feedback1*flowsigmaeff1+(1.0-feedback1)*sa
      flowsigmaeff1y=feedback1*flowsigmaeff1y+(1.0-feedback1)*say
      if (itt.eq.100*(itt/100)) write(3718,404) p0, itt, flowsigma1, flowsigmaeff1, sb, vecAria(istop)
 2277 continue
404 format(ES24.16,1X,I8,4(1X,ES24.16))
close(3718)
open(20,file='output.dat')
write(20,'(4ES27.17E3)') flowsigmaeff1,flowsigmaeff1y,AreaContact,sb
close(20)
end program
