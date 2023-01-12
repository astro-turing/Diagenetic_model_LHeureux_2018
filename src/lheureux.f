
!     Last change:  IL   19 Mar 2014    5:25 pm

!  Rythmite limestone/marl - semiimplicit diffusion on c,po and explicit advection on M with upwind scheme
!  Main code
       program marl_PDE
           implicit none
           real*8 pho(0:1000),cao(0:1000),coo(0:1000),ARAo(0:1000),CALo(0:1000)
       REAL*8  ph(0:1000),ca(0:1000),co(0:1000),ARA(0:1000),CAL(0:1000),U(0:1000),W(0:1000)
       REAL*8  phalf(0:1000),ARAhalf(0:1000),CALhalf(0:1000),cahalf(0:1000),cohalf(0:1000),whalf(0:1000)
       REAL*8 Rca(0:1000),Rco(0:1000), RAR(0:1000),RCAL(0:1000),RARdum(0:1000),RCALdum(0:1000)
       REAL*8 OC(0:1000),OA(0:1000),t,dt,dx, P(35),dca(0:1000),dco(0:1000),sigca(0:1000),sigco(0:1000)
       REAL*8 Rp(0:1000),te,wap,S(0:1000),Sdum(0:1000),Udum(0:1000),sigpo(0:1000)
           integer tmax,j ,outx,i,N,outt
           COMMON/general/ pho,cao,coo, ARAo,CALo,tmax,outx,outt, N
       COMMON/par/P
           open (8,file='amarlt1')
       OPEN(9,FILE='amarlt2')
       OPEN(10,FILE='amarlt3')
       OPEN(11,FILE='amarlt4')
       open (12,FILE='amarlx')
       OPEN(13,FILE='marlstuff')
       call init
       dt=P(15)
       dx=P(16)
       do 11 i=0,N
             ph(i)=pho(i)
             ca(i)=cao(i)
             co(i)=coo(i)
             ARA(i)=ARAo(i)
             CAL(i)=CALo(i)
11     continue
           t=0.0
       WRITE(12,101)'t ','AR','AR','AR','AR','CA','CA','CA','CA','P ','P ','P ','P ','ca','ca','ca','ca','co','co','co','co'&
 &             ,'U ','W '          
       WRITE(8,103) 'x ','AR','CA','Po','Ca','CO','TE','U ','W ','OC','OA','Rp'
       WRITE(9,103) 'x ','AR','CA','Po','Ca','CO','TE','U ','W ','OC','OA','Rp'
       WRITE(10,103) 'x ','AR','CA','Po','Ca','CO','TE','U ','W ','OC','OA','Rp'
       WRITE(11,103) 'x ','AR','CA','Po','Ca','CO','TE','U ','W ','OC','OA','Rp'
       write(13,110) 't ','p1','k1','u1','p2','k2','u2','p3','k3','u3'
       call auxf(t,n,ph,ca,co,ARA,CAL,U,S,W,OC,OA,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,RAR,RCAL)
           write (12,100) t,P(18),P(18),P(18),P(18),P(19),P(19),P(19),P(19),P(8),P(8),P(8),P(8),&
&           ca(N/4),ca(N/2),ca(3*N/4),ca(N),co(N/4),co(N/2),co(3*N/4),co(N),U(N/4),W(N/4)         
           wap=w(0)
           wap=wap    

       call projectX(n,ARA,CAL,U,S,RAR,RCAL,ARAhalf,CALhalf)
       call projectY(n,ph,ca,co,W,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,phalf,cahalf,cohalf)
           do 12 j=1,tmax
           IF(j/50000 *50000.eq.j) WRITE(6,*) 'doing t=',j*dt
 !          WRITE(6,*) 'doing j=',j
           call auxf(t,n,phalf,cahalf,cohalf,ARAhalf,CALhalf,Udum,Sdum,Whalf,OC,OA,dca,dco,sigpo,&
 &                      sigca,sigco,Rp,Rca,Rco,RARdum,RCALdum)   
                  call CRANK(t,n,ph,ca,co,Whalf,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,phalf)
           call upwind(n,ARA,CAL,U,S,RAR,RCAL)
           t=t+dt
           call auxf(t,n,ph,ca,co,ARA,CAL,U,S,W,OC,OA,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,RAR,RCAL)
 !  Output          
           IF(j/outt*outt.eq.j) then
              write (12,100) t,ara(N/4),ara(N/2),ara(3*N/4),ara(N),cal(N/4),cal(N/2),cal(3*N/4),cal(N),&
 &            ph(N/4),ph(N/2),ph(3*N/4),ph(N),ca(N/4),ca(N/2),ca(3*N/4),ca(N),&
 &            co(N/4),co(N/2),co(3*N/4),co(N),U(N/4),W(N/4)
           endif
               IF(j.eq.outx) then
              
               do 51 i=0,N
                 TE=1.-ara(i)-cal(i)
                 write (8,102) i*dx,ara(i),cal(i),ph(i),ca(i),co(i),TE,U(i),W(i),OC(i), OA(i),Rp(i)
51             continue
           endif
           IF(j.eq.2*outx) then
              
               do 52 i=0,N
                 TE=1-ara(i)-cal(i)
                 write (9,102) i*dx,ara(i),cal(i),ph(i),ca(i),co(i),TE,U(i),W(i),OC(i), OA(i),Rp(i)
52             continue
           endif
           IF(j.eq.3*outx) then
              
               do 53 i=0,N      
                 TE=1-ara(i)-cal(i)
                 write (10,102) i*dx,ara(i),cal(i),ph(i),ca(i),co(i),TE,U(i),W(i),OC(i), OA(i),Rp(i)
53         continue
           endif
           IF(j.eq.4*outx) then
              
               do 54 i=0,N     
                 TE=1-ara(i)-cal(i) 
                 write (11,102) i*dx,ara(i),cal(i),ph(i),ca(i),co(i),TE,U(i),W(i),OC(i), OA(i),Rp(i)
54         continue
          
           endif
           call projectX(n,ARA,CAL,U,S,RAR,RCAL,ARAhalf,CALhalf)
           call projectY(n,ph,ca,co,W,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,phalf,cahalf,cohalf)
12         continue
100        format (23(d12.5,1x))
101        format (13x,23(a2,11x))
102        format (12(d12.5,1x))
103        FORMAT(13x,12(a2,11x))
110        format(13x,10(a2,11x))
           close (8)
           CLOSE(9)
           CLOSE(10)
           CLOSE(11)
           CLOSE(12)
           CLOSE(13)
           write (6,*) 'fini'
           !pause
           stop 'marl'
           end
!
!
!
!
           subroutine auxf(t,n,ph,ca,co,ARA,CAL,U,S,W,OC,OA,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,RAR,RCAL)
! This routine calculates all auxiliairy functions: Peclet numbers, oversaturations, reaction rates and others
           implicit none
               real*8 P(35),ph(0:1000),ca(0:1000),co(0:1000),ARA(0:1000),CAL(0:1000),U(0:1000),W(0:1000)
           real *8 OC(0:1000),OA(0:1000),k(0:1000),dca(0:1000),dco(0:1000),sigca(0:1000),sigco(0:1000)
           REAL*8 Rp(0:1000),betasV,S(0:1000),sigpo(0:1000)
           REAL*8 Rca(0:1000),Rco(0:1000),RAR(0:1000),RCAL(0:1000)
           REAL*8 Pe1,Pe2,Pe3,sc,sa,t,rhosw0,F
           real*8 dpdx,dx
!           real*8 xpeak,xwidth,dx
              integer icem,icemf
          integer idis,idisf
          integer i,n
              COMMON/par/P 
          dx=P(16)
 !         eps=P(29)
          t=t 
          rhosw0=P(20)
!           dx=P(16)
!            xpeak=P(22)+0.5*P(14)
!           xwidth=0.5*P(14)
            idis=NINT(P(22)/P(16))
            idisf=NINT((P(22)+P(14))/P(16))
            icem=NINT(P(30)/P(16))
            icemf=NINT(P(31)/P(16))
 !           idisf=idisf+1
 !           Ko=phi0**3*(1-dexp(-10.*(1-phi0)/phi0))/(1-phi0)**2
             betasV=P(32)         
!    i = 0 
!  In fact, k is actually k*(1-ph)                
                    u(0)=P(23)
                    S(0)=1.
                    k(0)=ph(0)**3*betasV/(1-ph(0))
                    k(0)=k(0)*(1-dexp(-10.d0*(1-ph(0))/ph(0)))
                    dpdx=(ph(1)-ph(0))/dx
                    dpdx=0.
                    dpdx=dpdx
                    w(0)=u(0)-k(0)*(rhosw0-1)/ph(0)
 !                   w(0)=u(0)-k(0)*(rhosw0-1)/ph(0)*(1+dpdx*P(33)/(ph(0)-P(34)))
 !                   
 
 !                                    
                OA(0)=0.
                 sa=1-ca(0)*co(0)*P(24)
                 sc=ca(0)*co(0)-1.
              
!                 F=dexp(-(xpeak)**2/(2*xwidth**2))
                F=1.d0
                IF( (0.gt.idis).and. (0.lt.idisf)) then                         
                   if (sa.gt.0.) then
                         OA(0)=F*sa**P(1)
                   else
                         OA(0)=-F*(-sa)**P(1)
                         OA(0)=0. 
                   endif           
                endif
                IF( (0.eq.idis).or. (0.eq.idisf)) then                
                   if (sa.gt.0.) then 
                       OA(0)=sa**P(1)
                    else
                      OA(0)=-(-sa)**P(1)
                      OA(0)=0.
                    endif      
                endif
                IF( (0.lt.idis).or. (0.gt.idisf)) then 
                   if (sa.gt.0.) then 
                       OA(0)=0.d0
                    else
                       OA(0)=-(-sa)**P(1)
                       OA(0)=0.
                    endif      
                endif       
              if((0.gt.icem).and.(0.lt.icemf)) then               
                if (sc.gt.0.) then
                           OC(0)=sc**P(2)
                 else
                           OC(0)=-(-sc)**P(2)
                           OC(0)=0.
                 end if
              endif
              if((0.eq.icem).or.(0.eq.icemf)) then               
                if (sc.gt.0.) then
                           OC(0)=sc**P(2)
                 else
                           OC(0)=-(-sc)**P(2)
                           OC(0)=0.
                 end if
              endif 
              if((0.lt.icem).or.(0.gt.icemf)) then               
                if (sc.gt.0.) then
                           OC(0)=0.d0
                 else
                           OC(0)=-(-sc)**P(2)
                           OC(0)=0.
                 end if
              endif         
                       dca(0)=1./(1-2*dlog(ph(0)))
                       dco(0)=dca(0)*P(25)
           Rca(0)=P(26)*(P(6)/P(7)-ca(0))*(OA(0)*ARA(0)-P(3)*OC(0)*CAL(0))*(1-ph(0))/ph(0) 
           Rco(0)=P(26)*(P(6)/P(7)-co(0))*(OA(0)*ARA(0)-P(3)*OC(0)*CAL(0))*(1-ph(0))/ph(0) 
                      Pe1=w(0)*P(16)/(2*dca(0))
                      Pe2=w(0)*P(16)/(2*dco(0))
                      Pe3=w(0)*P(16)/(2*P(35))
                      if(Pe1.gt.0.) then
                         sigca(0)=(1.+dexp(-2*Pe1))/(1.-dexp(-2*Pe1))-1./Pe1
                      else
                         sigca(0)=(1.+dexp(2*Pe1))/(-1.+dexp(2*Pe1))-1./Pe1 
                      endif      
                      if(Pe2.gt.0.) then
                         sigco(0)=(1.+dexp(-2*Pe2))/(1.-dexp(-2*Pe2))-1./Pe2
                      else
                         sigco(0)=(1.+dexp(2*Pe2))/(-1.+dexp(2*Pe2))-1./Pe2 
                      endif    
                      if(Pe3.gt.0.) then
                         sigpo(0)=(1.+dexp(-2*Pe3))/(1.-dexp(-2*Pe3))-1./Pe3
                      else
                         sigpo(0)=(1.+dexp(2*Pe3))/(-1.+dexp(2*Pe3))-1./Pe3 
                      endif            
 !  other reaction rates
             if (1-ARA(0)-CAL(0).lt.1.d-70)then
               RAR(0)=-P(4)*P(26)*(1-CAL(0))-P(26)*OA(0)*(1-CAL(0))*CAL(0)-P(26)*P(3)*OC(0)*(1-CAL(0))*CAL(0)
               RCAL(0)=P(4)*P(26)*(1-CAL(0))+P(3)*P(26)*OC(0)*CAL(0)*(1-CAL(0))+P(26)*OA(0)*CAL(0)*(1-CAL(0))
               Rp(0)=P(26)*((1-CAL(0))*OA(0)-P(3)*CAL(0)*OC(0))*(1-ph(0))
             else
               RAR(0)=-P(4)*P(26)*ARA(0)-P(26)*OA(0)*ARA(0)*(1-ARA(0))-P(26)*P(3)*OC(0)*ARA(0)*CAL(0)
               RCAL(0)=P(4)*P(26)*ARA(0)+P(3)*P(26)*OC(0)*CAL(0)*(1-CAL(0))+P(26)*OA(0)*CAL(0)*ARA(0)
               Rp(0)=P(26)*(ARA(0)*OA(0)-P(3)*CAL(0)*OC(0))*(1-ph(0))
            endif                       
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                         
             do 82 i=1,N
!  supersaturations
!                F=dexp(-(i*dx-xpeak)**2/(2*xwidth**2))
                if(i.lt.N) then
                  dpdx=(ph(i+1)-ph(i-1))/(2*dx)
                else
                  dpdx=(ph(N)-ph(N-1))/dx
                endif    
                  dpdx=0.
                F=1.0d0
                OA(i)=0.
                sa=1-ca(i)*co(i)*P(24)
                 sc=ca(i)*co(i)-1
 
               
                IF( (i.gt.idis).and. (i.lt.idisf)) then        
                   if (sa.gt.0.) then
                         OA(i)=F*sa**P(1)
                   else
                         OA(i)=-F*(-sa)**P(1) 
                         OA(i)=0.
                   endif           
                endif
                IF( (i.eq.idis).or. (i.eq.idisf)) then
                   if (sa.gt.0.) then 
                       OA(i)=sa**P(1)
                    else
                      OA(i)=-(-sa)**P(1)
                      OA(i)=0.
                    endif      
                endif
                IF( (i.lt.idis).or. (i.gt.idisf)) then 
                   if (sa.gt.0.) then 
                       OA(i)=0.d0
                    else
                      OA(i)=-(-sa)**P(1)
                       OA(i)=0.
                    endif      
               endif       
              if((i.gt.icem).and.(i.lt.icemf)) then               
                if (sc.gt.0.) then
                           OC(i)=sc**P(2)
                 else
                           OC(i)=-(-sc)**P(2)
                           OC(i)=0.
                 end if
              endif
              if((i.eq.icem).or.(i.eq.icemf)) then               
                if (sc.gt.0.) then
                           OC(i)=sc**P(2)
                 else
                           OC(i)=-(-sc)**P(2)
                           OC(i)=0.
                 end if
              endif 
              if((i.lt.icem).or.(i.gt.icemf)) then               
                if (sc.gt.0.) then
                           OC(i)=0.d0
                 else
                           OC(i)=-(-sc)**P(2)
                           OC(i)=0.
                 end if
              endif        

! velocity
          if ((1-ph(i)).le.0.05)  then
!          write(6,*) 'p>0.95 for t,i',t,i                 
                        
 ! k *(1-p) actually                          
                           k(i)=betasV*10.*ph(i)**2
                           u(i)=k(i)*(rhosw0-1)+P(23)-P(13)*k(0)
 
                           
 !  approx u(i)=S                          
 !                          u(i)=P(23)
                           w(i)=u(i)-k(i)*(rhosw0-1)/ph(i)
 !                          w(i)=u(i)-k(i)*(rhosw0-1)/ph(i)*(1+dpdx*P(33)/(ph(i)-P(34)))
 !
    
  !                        
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                          
                          
                           dca(i)=1./(1-2*dlog(ph(i)))
                           dco(i)=dca(i)*P(25)
           Rca(i)=P(26)*(P(6)/P(7)-ca(i))*(OA(i)*ARA(i)-P(3)*OC(i)*CAL(i))*(1-ph(i))/ph(i) 
           Rco(i)=P(26)*(P(6)/P(7)-co(i))*(OA(i)*ARA(i)-P(3)*OC(i)*CAL(i))*(1-ph(i))/ph(i) 
                             
                      Pe1=w(i)*P(16)/(2*dca(i))
                      Pe2=w(i)*P(16)/(2*dco(i))
                      Pe3=w(i)*P(16)/(2*P(35))
                      if(Pe1.gt.0.) then
                         sigca(i)=(1.+dexp(-2*Pe1))/(1.-dexp(-2*Pe1))-1./Pe1
                      else
                         sigca(i)=(1.+dexp(2*Pe1))/(-1.+dexp(2*Pe1))-1./Pe1 
                      endif      
                      if(Pe2.gt.0.) then
                         sigco(i)=(1.+dexp(-2*Pe2))/(1.-dexp(-2*Pe2))-1./Pe2
                      else
                         sigco(i)=(1.+dexp(2*Pe2))/(-1.+dexp(2*Pe2))-1./Pe2 
                      endif 
                      If (Pe3.gt.0.) then
                         sigpo(i)=(1.+dexp(-2*Pe3))/(1.-dexp(-2*Pe3))-1./Pe3
                      else
                         sigpo(i)=(1.+dexp(2*Pe3))/(-1.+dexp(2*Pe3))-1./Pe3 
                      endif            
             
          else
                     
                   
                     if (ph(i).le.P(34)) THEN
 !                write(6,*) 'phi = eps at i, t=',ph(i),i,t
 !                pause
                       k(i)=0.
                       dca(i)=0.
                       dco(i)=0.
                       Rca(i)=0.
                       Rco(i)=0.
                       u(i)=k(i)*(rhosw0-1)+P(23)-P(13)*k(0)
                       
 ! approx  u(i)=S       
  !                     u(i)=P(23)             
                       w(i)=u(i)
 !!!!!!!!!!!!!!!!!!                      
                     
                       if(w(i).gt.0.) then
                          sigca(i)=1.
                          sigco(i)=1.
               
                       else  
                          sigca(i)=-1.
                          sigco(i)=-1.
                       endif
                       if(w(i).gt.0.) then
                          sigpo(i)=1.
                       else  
                          sigpo(i)=-1.
                       endif
                    else
                       k(i)=betasV*ph(i)**3/(1-ph(i))
                       k(i)=k(i)*(1-dexp(-10.*(1-ph(i))/ph(i)))
 ! diffusion coefficient
                       dca(i)=1./(1-2*dlog(ph(i)))
                       dco(i)=dca(i)*P(25)
 ! reaction rates 
           Rca(i)=P(26)*(P(6)/P(7)-ca(i))*(OA(i)*ARA(i)-P(3)*OC(i)*CAL(i))*(1-ph(i))/ph(i) 
           Rco(i)=P(26)*(P(6)/P(7)-co(i))*(OA(i)*ARA(i)-P(3)*OC(i)*CAL(i))*(1-ph(i))/ph(i) 
 !          write(13,*)i,Rca(i),Rco(i)
           
!  velocity u,w
                              u(i)=k(i)*(rhosw0-1)+P(23)-P(13)*k(0)
                      w(i)=u(i)-k(i)*(rhosw0-1)/ph(i)
!  approx  u(i)=S               
 !                     u(i)=P(23)     

                      
                      w(i)=P(23)-k(i)*(rhosw0-1)/ph(i)
!                       w(i)=u(i)-k(i)*(rhosw0-1)/ph(i)*(1+dpdx*P(33)/(ph(i)-P(34)))
  !                    
                              
                   
 !!!!!!!!!!!!!!!!!!!!                     
                      Pe1=w(i)*P(16)/(2*dca(i))
                      Pe2=w(i)*P(16)/(2*dco(i))
                      Pe3=w(i)*P(16)/(2*P(35))
                      if(Pe1.gt.0.) then
                         sigca(i)=(1.+dexp(-2*Pe1))/(1.-dexp(-2*Pe1))-1./Pe1
                      else
                         sigca(i)=(1.+dexp(2*Pe1))/(-1.+dexp(2*Pe1))-1./Pe1 
                      endif      
                      if(Pe2.gt.0.) then
                         sigco(i)=(1.+dexp(-2*Pe2))/(1.-dexp(-2*Pe2))-1./Pe2
                      else
                         sigco(i)=(1.+dexp(2*Pe2))/(-1.+dexp(2*Pe2))-1./Pe2 
                      endif 
                       if(Pe3.gt.0.) then
                         sigpo(i)=(1.+dexp(-2*Pe3))/(1.-dexp(-2*Pe3))-1./Pe3
                      else
                         sigpo(i)=(1.+dexp(2*Pe3))/(-1.+dexp(2*Pe3))-1./Pe3 
                      endif               
             
              endif          
            endif
 !  other reaction rates
             if (1-ARA(i)-CAL(i).lt.1.d-70)then
               RAR(i)=-P(4)*P(26)*(1-CAL(i))-P(26)*OA(i)*(1-CAL(i))*CAL(i)-P(26)*P(3)*OC(i)*(1-CAL(i))*CAL(i)
               RCAL(i)=P(4)*P(26)*(1-CAL(i))+P(3)*P(26)*OC(i)*CAL(i)*(1-CAL(i))+P(26)*OA(i)*CAL(i)*(1-CAL(i))
               Rp(i)=P(26)*((1-CAL(i))*OA(i)-P(3)*CAL(i)*OC(i))*(1-ph(i))
             else
               RAR(i)=-P(4)*P(26)*ARA(i)-P(26)*OA(i)*ARA(i)*(1-ARA(i))-P(26)*P(3)*OC(i)*ARA(i)*CAL(i)
               RCAL(i)=P(4)*P(26)*ARA(i)+P(3)*P(26)*OC(i)*CAL(i)*(1-CAL(i))+P(26)*OA(i)*CAL(i)*ARA(i)
               Rp(i)=P(26)*(ARA(i)*OA(i)-P(3)*CAL(i)*OC(i))*(1-ph(i))
            endif 
! debug          
!            u(i)=P(23)
            
             
            S(i)=u(i)/dabs(u(i))                        
            
                   
82         continue
!           write (13,111) t,ph(N/4),k(N/4),u(N/4),w(N/4)
!111        format (5(d12.5,1x))
           return
           end
!

           subroutine projectX(n,ARA,CAL,U,S,RAR,RCAL,ARAhalf,CALhalf)
! This routine calculates the projected solid phase fields at half-time in order to treat the non-linearities
           real*8 U(0:1000),ARA(0:1000),CAL(0:1000),ARAhalf(0:1000),CALhalf(0:1000)
           REAL*8 S(0:1000)
           REAL*8 dt,dx,P(35),a
           REAL*8 RAR(0:1000),RCAL(0:1000)
        
               integer N,i
           COMMON/par/P
           dx=P(16)
           dt=P(15)
           
           a=dt/(2*dx)
        
           
       
             ARAhalf(0)=ARA(0)
             CALhalf(0)=CAL(0)
          
                      

           do 30 i=1,n-1
         
           
             ARAhalf(i)=ARA(i)-a*u(i)*(ARA(i)*S(i)-ARA(i-1)*0.5*(S(i)+1.)+ARA(i+1)*0.5*(1.-S(i)))&
 &              +0.5*dt*RAR(i)
             CALhalf(i)=CAL(i)-a*u(i)*(CAL(i)*S(i)-CAL(i-1)*0.5*(S(i)+1.)+CAL(i+1)*0.5*(1.-S(i)))&
 &               +0.5*dt*RCAL(i)
             if(1-ARAhalf(i)-CALhalf(i).lt.1.d-70) ARAhalf(i)=1-CALhalf(i)
             if(1-CALhalf(i).lt.1.d-10) CALhalf(i)=1.
             if(1-ARAhalf(i).lt.1.d-10) ARAhalf(i)=1.
             if(ARAhalf(i).lt.1.d-70) ARAhalf(i)=0.
             if(CALhalf(i).lt.1.d-70) CALhalf(i)=0.
        
30         continue
 
             ARAhalf(n)=ARA(n)-a*u(n)*S(n)*(ARA(n)-ARA(n-1))+0.5*dt*RAR(n)
             CALhalf(n)=CAL(n)-a*u(n)*S(n)*(CAL(n)-CAL(n-1))+0.5*dt*RCAL(n)
            
             if(1-ARAhalf(n)-CALhalf(n).lt.1.d-70) ARAhalf(n)=1-CALhalf(n)
             if(1-CALhalf(n).lt.1.d-10) CALhalf(n)=1.
             if(1-ARAhalf(n).lt.1.d-10) ARAhalf(n)=1.
             if(ARAhalf(n).lt.1.d-70) ARAhalf(n)=0.
             if(CALhalf(n).lt.1.d-70) CALhalf(n)=0.
         
           return
           end
!
!
           SUBROUTINE  projectY(n,ph,ca,co,W,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,phalf,cahalf,cohalf)
! This routine calculates the projected aqueous phase fields at half-time in order to treat the non-linearities
           real*8 ph(0:1000),W(0:1000),cahalf(0:1000),cohalf(0:1000),phalf(0:1000),sigpo(0:1000)
           REAL*8 ca(0:1000),co(0:1000),dca(0:1000),dco(0:1000),sigca(0:1000),sigco(0:1000),Rca(0:1000),Rco(0:1000)
           REAL*8 dt,dx,P(35),a,b,eps,difpor,Rp(0:1000)
           integer N,i
           COMMON/par/P
           dx=P(16)
           dt=P(15)
           eps=P(29)
           a=dt/(4*dx)
           b=dt/(4*dx*dx)
           difpor=P(35)
          
           cahalf(0)=ca(0)
           cohalf(0)=co(0)
           phalf(0)=ph(0)
           do 30 i=1,n-1
             phalf(i)=ph(i)-a*((1-sigpo(i))*ph(i+1)*w(i+1)+2*sigpo(i)*ph(i)*w(i)-&
&          (1+sigpo(i))*ph(i-1)*w(i-1))&
&             +2*b*difpor*(ph(i-1)-2*ph(i)+ph(i+1))+0.5*dt*Rp(i)
             if(phalf(i).lt.eps) phalf(i)=eps
             if(1-phalf(i).lt.eps) phalf(i)=1.-eps    
!               
             if(ph(i).le.eps) then
                cahalf(i)=ca(i)-a*w(i)*((1-sigca(i))*ca(i+1)+2*sigca(i)*ca(i)-(1+sigca(i))*ca(i-1)) &
&             +0.5*dt*Rca(i)
                cohalf(i)=co(i)-a*w(i)*((1-sigco(i))*co(i+1)+2*sigco(i)*co(i)-(1+sigco(i))*co(i-1)) &
&             +0.5*dt*Rco(i)
           else
            cahalf(i)=ca(i)-a*w(i)*((1-sigca(i))*ca(i+1)+2*sigca(i)*ca(i)-(1+sigca(i))*ca(i-1)) &
&             +b*(ph(i+1)*dca(i+1)+ph(i)*dca(i))*(ca(i+1)-ca(i))/ph(i) &
&             -b*(ph(i-1)*dca(i-1)+ph(i)*dca(i))*(ca(i)-ca(i-1))/ph(i)+0.5*dt*Rca(i)
            cohalf(i)=co(i)-a*w(i)*((1-sigco(i))*co(i+1)+2*sigco(i)*co(i)-(1+sigco(i))*co(i-1)) &
&             +b*(ph(i+1)*dco(i+1)+ph(i)*dco(i))*(co(i+1)-co(i))/ph(i) &
&             -b*(ph(i-1)*dco(i-1)+ph(i)*dco(i))*(co(i)-co(i-1))/ph(i)+0.5*dt*Rco(i)
            endif
             if(cahalf(i).lt.1.d-15) cahalf(i)=0.
             if(cohalf(i).lt.1.d-15) cohalf(i)=0.
             
          
30         continue
              phalf(n)=ph(n)+2*a*(sigpo(n)*ph(n-1)*w(n-1)-sigpo(n)*ph(n)*w(n)) &
&             +4*b*difpor*(ph(n-1)-ph(n))+0.5*dt*Rp(n)
             if(phalf(n).lt.eps) phalf(n)=eps
             if(1-phalf(n).lt.eps) phalf(n)=1.-eps    
!               
             if (ph(n).le.eps) then
!          
                cahalf(n)=ca(n)-2*a*w(n)*sigca(n)*(ca(n)-ca(n-1)) +0.5*dt*Rca(n)
           cohalf(n)=co(n)-2*a*w(n)*sigco(n)*(co(n)-co(n-1)) +0.5*dt*Rco(n)
             else  
           cahalf(n)=ca(n)-2*a*w(n)*sigca(n)*(ca(n)-ca(n-1)) &
&             +2*b*(ph(n-1)*dca(n-1)+ph(n)*dca(n))*(ca(n-1)-ca(n))/ph(n)+0.5*dt*Rca(n)
           cohalf(n)=co(n)-2*a*w(n)*sigco(n)*(co(n)-co(n-1)) &
&             +2*b*(ph(n-1)*dco(n-1)+ph(n)*dco(n))*(co(n-1)-co(n))/ph(n)+0.5*dt*Rco(n)
             endif
             if(cahalf(n).lt.1.d-15) cahalf(n)=0.
             if(cohalf(n).lt.1.d-15) cohalf(n)=0.
 
           return
           end
!

           subroutine upwind(n,ARA,CAL,U,S,RAR,RCAL)
! Implementation of the upwind scheme for the advection equations
           real*8 U(0:1000),ARA(0:1000),CAL(0:1000)
           REAL*8 RAR(0:1000),RCAL(0:1000),ARAnew(0:1000),CALnew(0:1000)
           REAL*8 S(0:1000)
           REAL*8 dt,dx,P(35),a
          
               integer N,i
           COMMON/par/P
           dx=P(16)
           dt=P(15)
           
           a=dt/dx
         
              
           ARAnew(0)=ARA(0)
           CALnew(0)=CAL(0)

           do 30 i=1,n-1
         
         
             ARANEW(i)=ARA(i)-a*u(i)*(ARA(i)*S(i)-ARA(i-1)*0.5*(S(i)+1.)+ARA(i+1)*0.5*(1.-S(i)))&
 &              +dt*RAR(i)
             CALNEW(i)=CAL(i)-a*u(i)*(CAL(i)*S(i)-CAL(i-1)*0.5*(S(i)+1.)+CAL(i+1)*0.5*(1.-S(i)))&
 &               +dt*RCAL(i)
             if(1-ARAnew(i)-CALnew(i).lt.1.d-70) ARAnew(i)=1-CALnew(i)
             if(1-CALnew(i).lt.1.d-10) CALnew(i)=1.
             if(1-ARAnew(i).lt.1.d-10) ARAnew(i)=1.
             if(ARAnew(i).lt.1.d-70) ARAnew(i)=0.
             if(CALnew(i).lt.1.d-70) CALnew(i)=0.
                   
30         continue

             ARAnew(n)=ARA(n)-a*u(n)*S(n)*(ARA(n)-ARA(n-1))+dt*RAR(n)
             CALnew(n)=CAL(n)-a*u(n)*S(n)*(CAL(n)-CAL(n-1))+dt*RCAL(n)
             
             
             if(1-ARAnew(n)-CALnew(n).lt.1.d-70) ARAnew(n)=1-CALnew(n)
             if(1-CALnew(n).lt.1.d-10) CALnew(n)=1.
             if(1-ARAnew(n).lt.1.d-10) ARAnew(n)=1.
             if(ARAnew(n).lt.1.d-70) ARAnew(n)=0.
             if(CALnew(n).lt.1.d-70) CALnew(n)=0.
          
            
!    update
             do 44 i=0,n
                
                ARA(i)=ARAnew(i)
                CAL(i)=CALnew(i)
44           continue                                 

            
           return
           end
!
!
           subroutine CRANK(t,n,ph,ca,co,Whalf,dca,dco,sigpo,sigca,sigco,Rp,Rca,Rco,phalf)
! implementation of the Crank-Nicholson scheme for the advection-diffusion equations: X(n+1)=(A^-1).B.X(n)
           real*8 ph(0:1000),ca(0:1000),co(0:1000),Whalf(0:1000),sigpo(0:1000),phalf(0:1000)
           REAL*8 dt,dx,P(35),sigca(0:1000),sigco(0:1000),dca(0:1000),dco(0:1000),Rca(0:1000),Rco(0:1000)
           REAL*8 aA1(0:1000),bA1(0:1000),cA1(0:1000),aA2(0:1000),bA2(0:1000),cA2(0:1000),Byca(0:1000),Byco(0:1000)
           REAL*8 aA3(0:1000),bA3(0:1000),cA3(0:1000),Bypo(0:1000),Rp(0:1000)
           REAL*8 a,b,t,eps
           integer N,i
           COMMON/par/P
           dx=P(16)
           dt=P(15)
           eps=P(29)
           a=dt/(4*dx)
         
           b=dt/(4*dx*dx)
           call matA(n,a,b,phalf,whalf,sigpo,sigca,sigco,dca,dco,aA1,bA1,cA1,aA2,bA2,cA2,aA3,bA3,cA3)
           call matBy(n,dt,a,b,ph,ca,co,phalf,whalf,sigpo,sigca,sigco,dca,dco,Rp,Rca,Rco,Bypo,Byca,Byco)
 !          if (t.gt.320.2) then
 !             WRITE(13,*) t,aA1(25),aA1(26),bA1(25),bA1(26),cA1(25),cA1(26)
  !            WRITE(13,*) t,aA2(25),aA2(26),bA2(25),bA2(26),cA2(25),cA2(26)
  !         endif
!  find solution
          call tridag(t,n,aA1,bA1,cA1,Byca,ca)
          call tridag(t,n,aA2,bA2,cA2,Byco,co)
      call tridag(t,n,aA3,bA3,cA3,Bypo,ph)
          do 30 i=0,n
              IF(ca(i).lt.1.d-15) ca(i) = 0.
              IF(co(i).lt.1.d-15) co(i) = 0.
              if(ph(i).lt.eps) ph(i)=eps
             if(1-ph(i).lt.eps) ph(i)=1.-eps    
               
30        continue
          return
          end
  !
  !
  !
          subroutine matA(n,a,b,phalf,whalf,sigpo,sigca,sigco,dca,dco,aA1,bA1,cA1,aA2,bA2,cA2,aA3,bA3,cA3)
 ! Determine the matrix elements of the tridiagonal matrix A in the Crank-Nicholson scheme
          real*8 Whalf(0:1000),phalf(0:1000)
           REAL*8 sigpo(0:1000),sigca(0:1000),sigco(0:1000),dca(0:1000),dco(0:1000)
           REAL*8 aA1(0:1000),bA1(0:1000),cA1(0:1000),aA2(0:1000),bA2(0:1000),cA2(0:1000)
           Real*8 aA3(0:1000),bA3(0:1000),cA3(0:1000)
           REAL*8 a,b,eps,P(35),difpor
           integer N,i
           common/par/P
           eps=P(29)
           difpor=P(35)
           aA1(0)=0.
           bA1(0)=1.
           cA1(0)=0.
           aA2(0)=0.
           bA2(0)=1.
           cA2(0)=0.
           aA3(0)=0.
           bA3(0)=1.
           cA3(0)=0.
          do 10 i=1,n-1
           aA3(i)=-a*whalf(i-1)*(1+sigpo(i-1))-2*difpor*b
           bA3(i)=1.+2*a*whalf(i)*sigpo(i)+4*difpor*b
           cA3(i)=a*whalf(i+1)*(1-sigpo(i+1))-2*difpor*b
 !     
           if (phalf(i).le.eps) then
           aA1(i)= -a*whalf(i)*(1+sigca(i))
           bA1(i)=1.+a*2*whalf(i)*sigca(i)
           cA1(i)=a*whalf(i)*(1-sigca(i))
           aA2(i)= -a*whalf(i)*(1+sigco(i))
           bA2(i)=1.+a*2*whalf(i)*sigco(i)
           cA2(i)= a*whalf(i)*(1-sigco(i))

           else  
               aA1(i)= -a*whalf(i)*(1+sigca(i))-b*(dca(i-1)*phalf(i-1)+dca(i)*phalf(i))/phalf(i)
           bA1(i)=1.+a*2*whalf(i)*sigca(i)+b*(2*phalf(i)*dca(i)+phalf(i+1)*dca(i+1)+phalf(i-1)*dca(i-1))/phalf(i)
           cA1(i)=a*whalf(i)*(1-sigca(i))-b*(dca(i+1)*phalf(i+1)+dca(i)*phalf(i))/phalf(i)
           aA2(i)= -a*whalf(i)*(1+sigco(i))-b*(dco(i-1)*phalf(i-1)+dco(i)*phalf(i))/phalf(i)
           bA2(i)=1.+a*2*whalf(i)*sigco(i)+b*(2*phalf(i)*dco(i)+phalf(i+1)*dco(i+1)+phalf(i-1)*dco(i-1))/phalf(i)
           cA2(i)= a*whalf(i)*(1-sigco(i))-b*(dco(i+1)*phalf(i+1)+dco(i)*phalf(i))/phalf(i)
           endif
10       continue
           aA3(n)=-2*a*whalf(n-1)*sigpo(n-1)-4*difpor*b
           bA3(n)=1.+2*a*whalf(n)*sigpo(n)+4*difpor*b
           cA3(n)=1.
 !          
           if(phalf(n).le.eps) then
           aA1(n)= -2*a*whalf(n)*sigca(n)
           bA1(n)=1.+2*a*whalf(n)*sigca(n)
           cA1(n)=1.
           aA2(n)= -2*a*whalf(n)*sigco(n)
           bA2(n)=1.+2*a*whalf(n)*sigco(n)
           cA2(n)=1.
           else
           aA1(n)= -2*a*whalf(n)*sigca(n)-b*(dca(n-1)*phalf(n-1)+3*dca(n)*phalf(n))/phalf(n)
           bA1(n)=1.+2*a*whalf(n)*sigca(n)+b*(3*dca(n)*phalf(n)+phalf(n-1)*dca(n-1))/phalf(n)
           cA1(n)=1.
           aA2(n)= -2*a*whalf(n)*sigco(n)-b*(dco(n-1)*phalf(n-1)+3*dco(n)*phalf(n))/phalf(n)
           bA2(n)=1.+2*a*whalf(n)*sigco(n)+b*(3*dco(n)*phalf(n)+phalf(n-1)*dco(n-1))/phalf(n)
           cA2(n)=1.
           endif
          return
          end
!
!
          subroutine matby(n,dt,a,b,ph,ca,co,phalf,whalf,sigpo,sigca,sigco,dca,dco,Rp,Rca,Rco,Bypo,Byca,Byco)
! Determine the matrix elements of the tridiagonal matrix B in the Crank-Nicholson scheme
           REAL*8 phalf(0:1000),whalf(0:1000),sigpo(0:1000), sigca(0:1000),sigco(0:1000),dca(0:1000)
       REAL*8 dco(0:1000),Rca(0:1000),Rco(0:1000),Rp(0:1000)
           REAL*8 Byca(0:1000),Byco(0:1000),Bypo(0:1000),co(0:1000),ca(0:1000),ph(0:1000)
           REAL*8 a,b,dt,P(35),eps,difpor
           integer N,i
           common/par/P
          eps=P(29)
          difpor=P(35)
          Byca(0)=ca(0)
          Byco(0)=co(0)
          Bypo(0)=ph(0)


          do 10 i=1,n-1
       Bypo(i)=(a*whalf(i-1)*(1+sigpo(i))+2*b*difpor)*ph(i-1) &
 &             +(1.-a*2*whalf(i)*sigpo(i)-4*b*difpor)*ph(i) &
 &             +  (-a*whalf(i+1)*(1-sigpo(i))+2*b*difpor)*ph(i+1) + dt*Rp(i)
      
 !
      if(phalf(i).le.eps) then
        Byca(i)= (a*whalf(i)*(1+sigca(i)))*ca(i-1) &
 &             +(1.-a*2*whalf(i)*sigca(i))*ca(i) &
 &             +  (-a*whalf(i)*(1-sigca(i)))*ca(i+1) + dt*Rca(i)
           Byco(i)=(a*whalf(i)*(1+sigco(i)))*co(i-1) &
 &             +(1.-a*2*whalf(i)*sigco(i))*co(i) &
 &             +  (-a*whalf(i)*(1-sigco(i)))*co(i+1) + dt*Rco(i)
      else 
           Byca(i)= (a*whalf(i)*(1+sigca(i))+b*(dca(i-1)*phalf(i-1)+dca(i)*phalf(i))/phalf(i))*ca(i-1) &
 &             +(1.-a*2*whalf(i)*sigca(i)-b*(2*dca(i)*phalf(i)+dca(i+1)*phalf(i+1)+dca(i-1)*phalf(i-1))/phalf(i))*ca(i) &
 &             +  (-a*whalf(i)*(1-sigca(i))+b*(dca(i+1)*phalf(i+1)+dca(i)*phalf(i))/phalf(i))*ca(i+1) + dt*Rca(i)
           Byco(i)=(a*whalf(i)*(1+sigco(i))+b*(dco(i-1)*phalf(i-1)+dco(i)*phalf(i))/phalf(i))*co(i-1) &
 &             +(1.-a*2*whalf(i)*sigco(i)-b*(2*dco(i)*phalf(i)+dco(i+1)*phalf(i+1)+dco(i-1)*phalf(i-1))/phalf(i))*co(i) &
 &             +  (-a*whalf(i)*(1-sigco(i))+b*(dco(i+1)*phalf(i+1)+dco(i)*phalf(i))/phalf(i))*co(i+1) + dt*Rco(i)
         endif
10       continue
         Bypo(n)=(2*a*whalf(n-1)*sigpo(n-1)+4*b*difpor)*ph(n-1) &
 &             +(1.-a*2*whalf(n)*sigpo(n)-4*b*difpor)*ph(n) + dt*Rp(n)
!
         if(phalf(n).le.eps) then
                    Byca(n)= 2*a*whalf(n)*sigca(n)*ca(n-1) &
 &             +(1.-a*2*whalf(n)*sigca(n))*ca(n)+dt*Rca(n)
            Byco(n)=2*a*whalf(n)*sigco(n)*co(n-1) &
 &             +(1.-a*2*whalf(n)*sigco(n))*co(n)+dt*Rco(n)
         else            
           Byca(n)= (2*a*whalf(n)*sigca(n)+b*(dca(n-1)*phalf(n-1)+3*dca(n)*phalf(n))/phalf(n))*ca(n-1) &
 &             +(1.-a*2*whalf(n)*sigca(n)-b*(3*dca(n)*phalf(n)+dca(n-1)*phalf(n-1))/phalf(n))*ca(n)+dt*Rca(n)
            Byco(n)=(2*a*whalf(n)*sigco(n)+b*(dco(n-1)*phalf(n-1)+3*dco(n)*phalf(n))/phalf(n))*co(n-1) &
 &             +(1.-a*2*whalf(n)*sigco(n)-b*(3*dco(n)*phalf(n)+dco(n-1)*phalf(n-1))/phalf(n))*co(n)+dt*Rco(n)
         endif
         return
          end
!
!
        subroutine tridag(t,m,a,b,c,r,sol)
! Find the inverse of the tridiagonal matrix A by Thomas algorithm
          implicit none
          real*8 a(0:1000),b(0:1000),c(0:1000),r(0:1000),sol(0:1000)
          real*8 h(0:1000),g(0:1000),denom,t
          INTEGER m,i
          g(m-1)=-a(m)/b(m)
          h(m-1)=r(m)/b(m)
          do 10 i=m-2,0,-1
             denom=b(i+1)+c(i+1)*g(i+1)
                 if (dabs(denom).le.1.d-70) then
                    write (6,*)  'tridag a echoue a i et t=',i,t
                    WRITE(6,*) 'denom,b,c,g=',denom,b(i+2),c(i+2),g(i+2),b(i+1),c(i+1),g(i+1)
                        !pause
                 endif        
                 g(i)=-a(i+1)/denom
                 h(i)=(r(i+1)-c(i+1)*h(i+1))/denom
10     continue
        sol(0)=(r(0)-c(0)*h(0))/(b(0)+c(0)*g(0))
          do 11 i=1,m
11       sol(i)=sol(i-1)*g(i-1)+h(i-1)
       return
          end


!
!
       subroutine init
! This routine defines all parameter values and initial conditions, to be passed in vector P.
           implicit none
           real*8 pho(0:1000),cao(0:1000),coo(0:1000),ARAo(0:1000),CALo(0:1000)
           REAL*8 dt,dx, P(35),Th,mua,rhoa,rhoc,rhot,rhow,D0ca,D0co3,Ka,Kc,beta,k1,k2,k3,length,xdis
           REAL*8 m,nn,S,phi0,ca0,co30,ccal0,cara0,a,b,Vscale,rhos0,Xs,Ts,eps,xcem,xcemf
           real*8 phi00,ca00,co300,ccal00,cara00,phiinf,cAthy
!           REAL*8 arg,tgh
           integer tmax,N,outt,i,outx
           COMMON/general/ pho,cao,coo, ARAo,CALo,tmax,outx,outt, N
           COMMON/par/P
!   DIMENSIONLESS TIMES, tmax=max time index, outx = time index for graphic x output
!   outt= time index for graphic t outputs at four points
!  dt, dx = discrete steps (dimensionless) ; N +1 = total number of nodes;
! N must be a multiple of 4
!  Th = thickness of active dissolution zone (cm); length = system size in cm, xdis=position of the beginning of the active dissolution zone (cm)
! mua in g/mol, densities in g/cm3, diff (in pure water) in cm2/a, Ka, Kc in M2
!
! beta in cm/a, k1, k2,k3 in 1/a, m, n dimensionless
! S in cm/a, phio,cal0,cara0 dimensionless, ca0, co30 in M, a and b (tangh profile for initial concentrations) in cm
! cAthy Athys constant in cm-1; phiinf=porosity at the bottom of the compressed system
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           LAX=0
           ! dt=1.d-3
           dt=5.d-6
        
           
 !          length=2000.
 
           xdis=50.
           xcem=-100.
           xcemf=1000.
           length=500.
 
           Th=100.
           
           eps=1.d-2
       
           tmax=200000
           outt=1000
           
           
           
          
           outx=tmax/4
           N=200
           mua=100.09
           rhoa=2.95
           rhoc=2.71
           rhot=2.8
           rhow=1.023
           D0ca=131.9
           D0co3=272.6
           Ka=10.**(-6.19)
           Kc=10.**(-6.37)
           

           beta=0.01
           k1=0.0
           k2=1.d-2
           
           k3=1.d-3
         
           m=2.48
           nn=2.8
 !          m=1.
  !         nn=1.
          
           
             
           S=0.01
          
           cAthy=0.1
         
!           phiinf=0.05
           phiinf=0.7
           phiinf=eps
!  new incoming sediment          
            phi0=0.7 
           ca0=0.5*dsqrt(Kc)
           co30=0.5*dsqrt(Kc)
           ccal0=0.3
           cara0=0.6
!  old uniform sediment
           phi00=0.7
           
           
           ca00=0.5*dsqrt(Kc)
           co300=0.5*dsqrt(Kc)
           ccal00=0.3
           cara00=0.6           
           a=0.7 *length
           b=length/10.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !          K0=beta*phi0**3/(1-phi0)**2
 ! Scale time with beta
           Vscale=beta
 ! Scale time with S   
           Vscale=S
                  
           rhos0=rhot/(1-cara0*(1-rhot/rhoa)-ccal0*(1-rhot/rhoc))
           Xs=D0ca/Vscale
           Ts=D0ca/Vscale/Vscale
           dx=(length/Xs)/N
!         
           P(1)= m
           P(2)= nn
           P(3)= k3/k2
           P(4)= k1/k2
           P(5)=cara0/ccal0
           P(6)= rhos0*ccal0/(mua*dsqrt(Kc))
           P(7)= rhos0*ccal0/rhoc
           P(8)=phi0
           P(9)= 1.-rhot/rhoc
           P(10)= 1.-rhot/rhoa
           P(11)= rhot/rhow
           P(12)=rhoc/rhow
           P(13)=rhos0/rhow-1.
           P(14)=Th/Xs
           P(15)=dt
           P(16)=dx
           P(17)=1-cara0-ccal0
           P(18)=cara0
           P(19)=ccal0
           P(20)=rhos0/rhow
           P(21)=a/Xs
           P(22)=xdis/Xs
           P(23)=S/Vscale
           P(24)=Kc/Ka
           P(25)= D0co3/D0ca
           P(26)=k2*Ts
           P(27)=rhoc/rhoa
           P(28)=b/Xs
           P(29)=eps
           P(30)=xcem/Xs         
           P(31)=xcemf/Xs
           P(32)=beta/Vscale
           P(33)=1./(cAthy*Xs)
           P(34)=phiinf
           P(35)=beta*phi0**3/((phi0-phiinf)*cAthy)/D0ca
!  case where all densities are equal to rhos
           P(7)=ccal0
           P(9)=0.
           P(10)=0.
           P(11)=rhos0/rhow
           P(12)=P(11)
           P(27)=1.  
    
! Initial conditions (dimensionless)

           pho(0)=phi0
           cao(0)=ca0/dsqrt(Kc)
           coo(0)=co30/dsqrt(Kc)
           ARAo(0)=cara0
           CALo(0)=ccal0
          do 10 i=1,N
           pho(i)=phi00
!           arg=(i*dx-P(21))/P(28)
!               IF(arg.ge.0.) tgh=(1-dexp(-2*arg))/(1+dexp(-2*arg))
!               IF(arg.lt.0.) tgh=(dexp(2*arg)-1)/(dexp(2*arg)+1)
!           cao(i)=ca0/dsqrt(Kc)*0.5*(1-tgh)
!           coo(i)=co30/dsqrt(Kc) *0.5*(1-tgh)
           cao(i)=ca00/dsqrt(Kc)
           coo(i)=co300/dsqrt(Kc)
           ARAo(i)=cara00
           CALo(i)=ccal00

 10     continue
!
           write (6,*) 'Xs (cm), Ts (a)', Xs,Ts
           write (6,*) 'dt,dx,dtS/dx =', dt,dx,P(15)*(P(23)/P(16))
           WRITE(6,*) 'dx^2/2d=', P(16)**2/(2*P(25))
           IF(P(15)*P(23)/P(16).gt.1./5.) then
              WRITE(6,*)' problem: possible instability'
              !pause
           endif
           write (6,*) 'scale for MA, MC =',rhos0*cara0*(1-phi0), rhos0*ccal0*(1-phi0)
           write (6,*) 'scale for c=',dsqrt(Kc)
           write (6,*) 'Damkohler number Da=, scaled sedimentation rate,rhos0/rhow=',P(26),P(23),rhos0/rhow
           WRITE(6,*) 'scaled length, position of dissolution zone=', length/Xs, xdis/Xs,(xdis+Th)/Xs
           WRITE(6,*) 'a/Xs,b/Xs, rhosw-1=',P(21),P(28),P(13)
           write(6,*) '1/cXs,Dpor/Dca=',P(33),P(35)



           !PAUSE
       return
       end

