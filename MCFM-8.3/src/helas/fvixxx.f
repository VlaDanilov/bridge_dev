c
c ----------------------------------------------------------------------
c
      subroutine fvixxx(fi,vc,g,fmass,fwidth , fvi)
c
c this subroutine computes an off-shell fermion wavefunction from a     
c flowing-in external fermion and a vector boson.                       
c                                                                       
c input:                                                                
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex vc(6)          : input    vector                      v 
c       real    g(2)           : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'             
c       real    fwidth         : width of output fermion f'             
c                                                                       
c output:                                                               
c       complex fvi(6)         : off-shell fermion             |f',v,fi>
c
      complex*16 fi(6),vc(6),fvi(6),sl1,sl2,sr1,sr2,d
      real*8    g(2),pf(0:3),fmass,fwidth,pf2
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
      complex*16 c_imag
      parameter( c_imag=(0d0,1d0) )
c
      fvi(5) = fi(5)-vc(5)
      fvi(6) = fi(6)-vc(6)
c
      pf(0)=dble( fvi(5))
      pf(1)=dble( fvi(6))
      pf(2)=dimag(fvi(6))
      pf(3)=dimag(fvi(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)
c
      d=-r_one/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),r_zero))
      sl1= (vc(1)+       vc(4))*fi(1)
     &    +(vc(2)-c_imag*vc(3))*fi(2)
      sl2= (vc(2)+c_imag*vc(3))*fi(1)
     &    +(vc(1)-       vc(4))*fi(2)
c
      if ( g(2) .ne. r_zero ) then
         sr1= (vc(1)-       vc(4))*fi(3)
     &       -(vc(2)-c_imag*vc(3))*fi(4)
         sr2=-(vc(2)+c_imag*vc(3))*fi(3)
     &       +(vc(1)+       vc(4))*fi(4)
c
         fvi(1) = ( g(1)*((pf(0)-pf(3))*sl1 -dconjg(fvi(6))*sl2)
     &             +g(2)*fmass*sr1)*d
         fvi(2) = ( g(1)*(      -fvi(6)*sl1 +(pf(0)+pf(3))*sl2)
     &             +g(2)*fmass*sr2)*d
         fvi(3) = ( g(2)*((pf(0)+pf(3))*sr1 +dconjg(fvi(6))*sr2)
     &             +g(1)*fmass*sl1)*d
         fvi(4) = ( g(2)*(       fvi(6)*sr1 +(pf(0)-pf(3))*sr2)
     &             +g(1)*fmass*sl2)*d
c
      else          
         fvi(1) = g(1)*((pf(0)-pf(3))*sl1 -dconjg(fvi(6))*sl2)*d
         fvi(2) = g(1)*(      -fvi(6)*sl1 +(pf(0)+pf(3))*sl2)*d
         fvi(3) = g(1)*fmass*sl1*d
         fvi(4) = g(1)*fmass*sl2*d
      end if
c
      return          
      end
