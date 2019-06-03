c
c ======================================================================
c
      subroutine iovxxx(fi,fo,vc,g , vertex)
c
c this subroutine computes an amplitude of the fermion-fermion-vector   
c coupling.                                                             
c                                                                       
c input:                                                                
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v 
c       real    g(2)           : coupling constants                  gvf
c                                                                       
c output:                                                               
c       complex vertex         : amplitude                     <fo|v|fi>
c
      complex*16 fi(6),fo(6),vc(6),vertex
      real*8    g(2)
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
      complex*16 c_imag
      parameter( c_imag=(0d0,1d0) )
c

      vertex =  g(1)*( (fo(3)*fi(1)+fo(4)*fi(2))*vc(1)
     &                +(fo(3)*fi(2)+fo(4)*fi(1))*vc(2)
     &                -(fo(3)*fi(2)-fo(4)*fi(1))*vc(3)*c_imag
     &                +(fo(3)*fi(1)-fo(4)*fi(2))*vc(4)        )
c
      if ( g(2) .ne. r_zero ) then
         vertex = vertex
     &          + g(2)*( (fo(1)*fi(3)+fo(2)*fi(4))*vc(1)
     &                  -(fo(1)*fi(4)+fo(2)*fi(3))*vc(2)
     &                  +(fo(1)*fi(4)-fo(2)*fi(3))*vc(3)*c_imag
     &                  -(fo(1)*fi(3)-fo(2)*fi(4))*vc(4)        )
      end if
c
      return
      end
