      subroutine op_fermion(y,t3,q,cc,ia,iz,ip,im,ia2,iz2,iw2,
     & cew,leg,ch)
      !------------------------------------------------------------
      ! the subroutine returns the operator values as follows:
      ! y -- hypercharge; q -- electric charge; t3 -- 3rd isospin
      ! cc -- SU(2) casimir; ia, iz, ip, im -- gauge couplings
      ! ia2, iz2, iw2 -- square of ia, iz, ip/m
      ! cew -- electroweak casimir operator
      ! leg -- external legs; ch -- chirality of the leg
      ! the legs use PDG code; ch = -1 right; ch = +1 left
      ! the matricies ip(2,2) and im(2,2) have the components which
      ! are off-diagonal, exchanging isospins. The order of the 
      ! indices are (u,d) and (nu,l) for fermion external leg.
      !------------------------------------------------------------
      implicit none
      include 'types.f'
      include 'ewcouple.f'
      include 'cplx.h'
      real(dp):: ia,iz,ip,im,ia2,iz2,iw2,cew,q,t3,cc,y
      integer leg,ch
      real(dp):: sw,cw,chir,dleg
      
      sw=sqrt(xw)
      cw=sqrt(1._dp-xw)
      chir=real(ch,dp)
      dleg=real(leg,dp)
!     ip=0._dp
!     im=0._dp
      
!     if(leg==0) then
!     write("particle # input error...please use the right corresponding #s")
!     stop
!     end if
      
      if (abs(leg) .gt. 0) then
c---  quarks and leptons
      
         if(abs(leg) >= 11 .and. abs(leg) <= 16) then
            if(mod(leg,2) == 0) then
               chir=1._dp
               q=0._dp
               t3=sign(1._dp/2._dp,dleg)
               y=2._dp*(q-t3)
            else
               q=-sign(1._dp,dleg)
               t3=-sign(1._dp/2._dp,dleg)*(chir+1._dp)/2._dp
               y=2._dp*(q-t3)
            end if
            ia=-q
            iz=(t3-xw*q)/sw/cw
         else if(abs(leg) >= 1 .and. abs(leg) <= 6) then
            q=sign(1._dp/3._dp,dleg)*sign(1._dp,real(1-2*mod(abs(leg),2),dp))
     &           *(2._dp-mod(abs(leg),2))
            t3=sign(1._dp/2._dp,dleg)*(chir+1._dp)/2._dp
     &           *sign(1._dp,real(1-2*mod(abs(leg),2),dp))
            y=2._dp*(q-t3)
            ia=-q
            iz=(t3-xw*q)/sw/cw
         end if
         
      
         if((mod(leg,2) == 0 .and. leg > 0).or.(mod(leg,2) == -1)) then
            ip=1._dp/sqrt(2._dp)/sw*sign(1._dp,dleg)*(chir+1._dp)/2._dp
            im=0._dp
         else if(leg /= 0) then
            im=1._dp/sqrt(2._dp)/sw*sign(1._dp,dleg)*(chir+1._dp)/2._dp
            ip=0._dp
         else
            ip=0._dp
            im=0._dp
         end if
        
         cc=3._dp/4._dp*(chir+1._dp)/2._dp
         ia2=ia**2
         iz2=iz**2
         iw2=(cc-t3**2)/xw
         cew=ia2+iz2+iw2
      
      else
c-- gluons      
      
         ia=0._dp
         ia2=0._dp
         iz=0._dp
         iz2=0._dp
         ip=0._dp
         im=0._dp
         iw2=0._dp
         cew=0._dp
         cc=0._dp
         y=0._dp
         t3=0._dp
         q=0._dp
      
      endif
      
      end subroutine op_fermion
      
            
      
      subroutine op_scalar(ia,iz,ip,im,ia2,iz2,iw2,cew,leg)

C --- the # correspondence of the scalar fields
C --- phi^+ ---> 27 -> 1, i.e., read in # 27 and reassigned to 1 in subroutine
C --- phi^- ---> 28 -> 2
C --- H -------> 25 -> 3
C --- chi -----> 26 -> 4
      implicit none
      include 'types.f'
      include 'ewcouple.f'
      include 'cplx.h'
      complex(dp):: ia,iz,ip(2),im(2)
      real(dp):: ia2,iz2,iw2,cew,sw,cw
      integer leg,rleg

C-------------------------------------------------
C --- ia ---> photonic coupling 
C --- iz ---> Z boson clupling
C --- ip ---> W^+ coupling
C --- im ---> W^- coupling
C --- where ip and im have 2 components, which is
C --- necessary for phi^(+/-) coupled to two neutral
C --- scalar fields (H or chi) with different couplings;
C --- (e.g., I^+_(phi^+ --> H) and I^+_(phi^+ --> chi)
C --- have different values)
C --- while such couplings associated with scalar fields 
C --- only have one non-zero value.
C --- (e.g., I^+_(H --> phi^-))
C --- In the representation of charged scalar fields:
C --- (Ip(m)(1) --> H, Ip(m)(2) --> chi)
C --- for neutral fields:
C --- (Ip(m)(1) --> phi^+; Ip(m)(2) --> phi^-)
C ------------------------------------------------


      if(leg == 27) then
         rleg = 1
      elseif(leg == 28) then
         rleg = 2
      elseif(leg == 25) then
         rleg = 3
      elseif(leg == 26) then
         rleg = 4
      else
         write(*,*) "incorrect particle code for the scalar fields"
         write(*,*) "please input the following #: ",
     .        "26 ---> phi^+",
     .        "27 ---> phi^-",
     .        "25 ---> H",
     .        "26 ---> chi"
      endif

      sw=sqrt(xw)
      cw=sqrt(1._dp-xw)

      if(rleg == 1 .or. rleg == 2) then
         ia = (-1._dp,0._dp)
         iz = cplx1((1._dp-2._dp*xw)/sw/cw/2._dp)
         ip(1) = cplx2(+1._dp/2._dp/sw,0._dp)
         ip(2) = cplx2(0._dp,+1._dp/2._dp/sw)
         im(1) = cplx2(-1._dp/2._dp/sw,0._dp)
         im(2) = cplx2(0._dp,+1._dp/2._dp/sw)
         ia = ia*real(3-2*rleg,dp)
         iz = iz*real(3-2*rleg,dp)
         ip = ip*real(2-rleg,dp)
         im = im*real(rleg-1,dp)
      elseif(rleg == 3 .or. rleg == 4) then
         ia = (0._dp,0._dp)
         iz = cplx2(0._dp,-1._dp/2._dp/sw/cw)
         ip(1) = (0._dp,0._dp)
         ip(2) = + cplx2(-1._dp/2._dp/sw,0._dp)*real(4-rleg,dp)
     .           + cplx2(0._dp,-1._dp/2._dp/sw)*real(rleg-3,dp)
         im(1) = + cplx2(+1._dp/2._dp/sw,0._dp)*real(4-rleg,dp)
     .           + cplx2(0._dp,-1._dp/2._dp/sw)*real(rleg-3,dp)
         im(2) = (0._dp,0._dp)
         iz = iz*real(7-2*rleg,dp)

      endif
         
C --- iw2 is diagonal, more specifically is proportional to I_(2X2), 
C --- so the following expression w.r.t. iw2 = |I^+|^2+|I^-|^2
      ia2 = abs(ia)**2
      iz2 = abs(iz)**2
      iw2 = + abs(ip(1))**2 + abs(ip(2))**2
     .      + abs(im(1))**2 + abs(im(2))**2

      cew=ia2+iz2+iw2

C     return

      end subroutine op_scalar
      


      subroutine op_vector(ia,iz,ip,im,ia2,iz2,iw2,cew,leg)
C --- the subroutine returns the values of charge operators for external 
C --- gauge bosons which are transversely polarized
C --- There are two sub groups 1 = (W^+,W^-), 2 = (A,Z)
C --- particle code
C --- W^+/W^- ---> 24/-24
C --- A       ---> 22
C --- Z       ---> 23 
      implicit none
      include 'types.f'
      include 'ewcouple.f'
      real(dp):: ia,iz,ip(2),im(2),ia2,iz2,iw2(2),cew(2),sw,cw
      integer leg

      sw=sqrt(xw)
      cw=sqrt(1._dp-xw)

      if(abs(leg) == 24) then
         ia = -sign(1._dp,real(leg,dp))
         iz = +cw/sw*sign(1._dp,real(leg,dp))
         ip(1) = +1._dp*real(24+leg,dp)/4.8d1
         ip(2) = -cw/sw*real(24+leg,dp)/4.8d1
         im(1) = -1._dp*real(24-leg,dp)/4.8d1
         im(2) = +cw/sw*real(24-leg,dp)/4.8d1
         iw2(1) = 1._dp/xw
         iw2(2) = 0._dp
         cew(1) = 2._dp/xw
         cew(2) = 0._dp
      elseif(leg == 22) then
         ia = 0._dp
         iz = 0._dp
         ip(1) = 0._dp
         ip(2) = -1._dp
         im(1) = +1._dp
         im(2) = 0._dp
         iw2(1) = 2._dp
         iw2(2) = -2._dp*cw/sw
         cew = iw2
      elseif(leg == 23) then
         ia = 0._dp
         iz = 0._dp
         ip(1) = 0._dp
         ip(2) = +cw/sw
         im(1) = -cw/sw
         im(2) = 0._dp
         iw2(1) = -2._dp*cw/sw
         iw2(2) = 2._dp*cw**2/xw
         cew = iw2
      endif
         
      ia2 = ia**2
      iz2 = iz**2

      end subroutine op_vector
