      subroutine sudakov_DY(dd,born,legs,chs,rr)
! so far the considered process is 2 -> 2 NC, n=4
! the output is dd(5) with corresponding returns:
! dd(1) -- double logs from SU(2) symmetry
! dd(2) -- double logs from photon contribution in SU(2)
! dd(3) -- single logs from SU(2) symmetry
! dd(4) -- single logs from photon contribution in SU(2)
! dd(5) -- single logs from lz
! dd(6) -- single logs from lc+lyuk+lpr, proportional to log(musq/Mw**2)
! the scale musq is upto user's option, i.e., musq=Mw**2 makes these
! contributions vanish 
! dd(7) -- single logs from lcna, i.e., collinear w.r.t. photon
! born -- Yi*Yj/4/cw2+T3i*T3j/xw
      implicit none
      include 'types.f'
      include 'masses.f'
      include 'ewcouple.f'
      integer n
      parameter(n=4)
      real(dp):: dd(7),born
      integer legs(n),chs(n)
      real(dp):: rr(n,n)!rr(n*(n-1)/2)
      real(dp):: y(4),t3(4),q(4),cc(4),ia(4),iz(4),ip(4),im(4),
     & ia2(4),iz2(4),iw2(4),cew(4)
      real(dp):: cw2,dlm,iiz2,lssca,lsscz,lsscw,lc,lcna,lyuk,
     & lpr,mag(n,n)
      real(dp):: baz,baa,Rff,Dff
      integer i,j
      
      cw2=1._dp-xw
      dlm=log(zmass**2/wmass**2)
      baz=-(1.9d1+2.2d1*xw)/6._dp/sqrt(xw*cw2)
      baa=-11._dp/3._dp
      
      dd=0._dp
      iiz2=0._dp
      Rff=0._dp
      Dff=0._dp
      lssca=0._dp
      lsscz=0._dp
      lsscw=0._dp
      lc=0._dp
      lcna=0._dp
      lyuk=0._dp
      lpr=0._dp
      mag=0._dp


! calculate the quantum numbers of the doublet partner
      if(legs(2) == 6 .or. legs(2) == 16) then
         call op_fermion(y(2),t3(2),q(2),cc(2),ia(2),iz(2),ip(2),
     &        im(2),ia2(2),iz2(2),iw2(2),cew(2),legs(2)-1,chs(2))
      else
         call op_fermion(y(2),t3(2),q(2),cc(2),ia(2),iz(2),ip(2),
     &        im(2),ia2(2),iz2(2),iw2(2),cew(2),legs(2)+1,chs(2))
      endif

      call op_fermion(y(3),t3(3),q(3),cc(3),ia(3),iz(3),ip(3),im(3),
     &     ia2(3),iz2(3),iw2(3),cew(3),legs(3),chs(3))

      mag(1,2)=ia(2)*ia(3)+iz(2)*iz(3)

      call op_fermion(y(2),t3(2),q(2),cc(2),ia(2),iz(2),ip(2),im(2),
     &     ia2(2),iz2(2),iw2(2),cew(2),legs(2),chs(2))

      if(legs(3) == 6 .or. legs(3) == 16) then
         call op_fermion(y(3),t3(3),q(3),cc(3),ia(3),iz(3),ip(3),
     &        im(3),ia2(3),iz2(3),iw2(3),cew(3),legs(3)-1,chs(3))
      else
         call op_fermion(y(3),t3(3),q(3),cc(3),ia(3),iz(3),ip(3),
     &        im(3),ia2(3),iz2(3),iw2(3),cew(3),legs(3)+1,chs(3))
      endif

      mag(3,4)=ia(2)*ia(3)+iz(2)*iz(3)

!  mag(1,2)=-1._dp
!  mag(3,4)=-1._dp
! generally speaking, the values of mag(1,2) and mag(3,4) do not affect
! the overall result when one considers the s-channel scattering, 
! because log(r12/s)=0

      mag(1,3)=1._dp/2._dp/xw
      mag(2,4)=mag(1,3)
      mag(2,3)=mag(1,3)
      mag(1,4)=mag(1,3)

      do i=1,n
         call op_fermion(y(i),t3(i),q(i),cc(i),ia(i),iz(i),ip(i),im(i),
     &        ia2(i),iz2(i),iw2(i),cew(i),legs(i),chs(i))
         dd(1)=dd(1)-cew(i)/2._dp
         dd(2)=dd(2)-q(i)**2/2._dp
         iiz2=iiz2+iz2(i)
         lc=lc+3._dp/2._dp*cew(i)
         lcna=lcna+q(i)**2*3._dp/2._dp
         if(abs(legs(i)) == 6) then
            lyuk=lyuk-1._dp/8._dp/xw*(3._dp/2._dp-chs(i)/2._dp)
     &           *Mt**2/wmass**2
         else if(abs(legs(i)) == 5 .and. chs(i) == 1) then
            lyuk=lyuk-1._dp/8._dp/xw*Mt**2/wmass**2
         endif
      enddo

      do i=1,n-1
         do j=i+1,n
            lssca=lssca+2._dp*ia(i)*ia(j)*rr(i,j)
            lsscz=lsscz+2._dp*iz(i)*iz(j)*rr(i,j)
            lsscw=lsscw+2._dp*(ip(i)*im(j)+im(i)*ip(j))*rr(i,j)*mag(i,j)
         enddo
      enddo
  
      Rff=ia(2)*ia(3)+iz(2)*iz(3)
      Dff=-1._dp/4._dp/cw2*y(2)*y(3)+cw2/xw**2*t3(2)*t3(3)
      Dff=Dff/Rff
      lsscw=lsscw/Rff           ! in order to be proportional to NC born

      lpr=sqrt(xw/cw2)*baz*Dff-baa

      dd(3)=lssca+lsscz+lsscw
!     dd(4)=lcna+lssca
      dd(4)=lssca
      dd(5)=iiz2*dlm
      dd(6)=lc+lyuk+lpr
      dd(7)=lcna

      born=Rff**2
      
      !------------------------------------------------------------
      ! checking with numerical example in ref, eqs.(6.15) - (6.17)
      !  write(*,*) "DL: ", dd(1)
      !  write(*,*) "lz: ", iiz2*dlm
      !  write(*,*) "lssc: ", lssca+lsscz+lsscw
      !  write(*,*) "lc: ", lc
      !  write(*,*) "lyuk: ", lyuk
      !  write(*,*) "lpr: ", lpr
      !-------------------------------------------------------------
            
      end subroutine sudakov_DY
      
      
      
      
      
      subroutine mandelstam_var(rr,s,t,u,s12)
      implicit none
      include 'types.f'
      real(dp):: rr(4,4)
      real(dp):: s,t,u,s12
      
      rr=0._dp
      
      rr(1,2)=log(abs(s/s12))
      rr(3,4)=rr(1,2)
      
      rr(1,3)=log(abs(t/s12))
      rr(2,4)=rr(1,3)
      
      rr(2,3)=log(abs(u/s12))
      rr(1,4)=rr(2,3)
      
      end subroutine mandelstam_var



      subroutine sudakov_ttb(dd,legs,chs,rr)
C --- sudakov logarithms to ttb production
      implicit none
      include 'types.f'
      include 'masses.f'
      include 'ewcouple.f'
      integer, parameter:: n=4
      real(dp):: dd(7)
      integer legs(n),chs(n)
      real(dp):: rr(n,n)
      real(dp):: y(4),t3(4),q(4),cc(4),ia(4),iz(4),ip(4),im(4),
     & ia2(4),iz2(4),iw2(4),cew(4)
      real(dp):: cw2,dlm,iiz2,lssca,lsscz,lc,lcna,lyuk,
     & lpr,mag(n,n)
      integer i,j
      
      cw2=1._dp-xw
      dlm=log(zmass**2/wmass**2)
      
      dd=0._dp
      iiz2=0._dp
      lssca=0._dp
      lsscz=0._dp
      lc=0._dp
      lcna=0._dp
      lyuk=0._dp

      do i=1,n
         call op_fermion(y(i),t3(i),q(i),cc(i),ia(i),iz(i),ip(i),im(i),
     &        ia2(i),iz2(i),iw2(i),cew(i),legs(i),chs(i))
         dd(1)=dd(1)-cew(i)/2._dp
         dd(2)=dd(2)-q(i)**2/2._dp
         iiz2=iiz2+iz2(i)
         lc=lc+3._dp/2._dp*cew(i)
         lcna=lcna+q(i)**2*3._dp/2._dp
         if(abs(legs(i)) == 6) then
            lyuk=lyuk-1._dp/8._dp/xw*(3._dp/2._dp-chs(i)/2._dp)
     &           *Mt**2/wmass**2
         else if(abs(legs(i)) == 5 .and. chs(i) == 1) then
            lyuk=lyuk-1._dp/8._dp/xw*Mt**2/wmass**2
         endif
      enddo

      do i=1,n-1
         do j=i+1,n
            lssca=lssca+2._dp*ia(i)*ia(j)*rr(i,j)
            lsscz=lsscz+2._dp*iz(i)*iz(j)*rr(i,j)
         enddo
      enddo

      dd(3)=lssca+lsscz
!     dd(4)=lcna+lssca
      dd(4)=lssca
      dd(5)=iiz2*dlm
      dd(6)=lc+lyuk
      dd(7)=lcna
      
      end subroutine sudakov_ttb
      
