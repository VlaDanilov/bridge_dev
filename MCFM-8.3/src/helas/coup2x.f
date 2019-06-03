c
c ----------------------------------------------------------------------
c
      subroutine coup2x(sw2 , gal,gau,gad,gwf,gzn,gzl,gzu,gzd,g1)
c
c this subroutine sets up the coupling constants for the fermion-       
c fermion-vector vertices in the standard model.  the array of the      
c couplings specifies the chirality of the flowing-in fermion.  g??(1)  
c denotes a left-handed coupling, and g??(2) a right-handed coupling.   
c                                                                       
c input:                                                                
c       real    sw2            : square of sine of the weak angle       
c                                                                       
c output:                                                               
c       real    gal(2)         : coupling with a of charged leptons     
c       real    gau(2)         : coupling with a of up-type quarks      
c       real    gad(2)         : coupling with a of down-type quarks    
c       real    gwf(2)         : coupling with w-,w+ of fermions        
c       real    gzn(2)         : coupling with z of neutrinos           
c       real    gzl(2)         : coupling with z of charged leptons     
c       real    gzu(2)         : coupling with z of up-type quarks      
c       real    gzd(2)         : coupling with z of down-type quarks    
c       real    g1(2)          : unit coupling of fermions              
c
      real*8 gal(2),gau(2),gad(2),gwf(2),gzn(2),gzl(2),gzu(2),gzd(2),
     &     g1(2),sw2,alpha,fourpi,ee,sw,cw,ez,ey
c
      real*8 r_zero, r_half, r_one, r_two, r_three, r_four, r_ote
      real*8 r_pi, r_ialph
      parameter( r_zero=0.0d0, r_half=0.5d0, r_one=1.0d0, r_two=2.0d0,
     $     r_three=3.0d0 )
      parameter( r_four=4.0d0, r_ote=128.878d0 )
      parameter( r_pi=3.14159265358979323846d0, r_ialph=137.0359895d0 )
c
      alpha = r_one / r_ote
c      alpha = r_one / r_ialph
      fourpi = r_four * r_pi
      ee=sqrt( alpha * fourpi )
      sw=sqrt( sw2 )
      cw=sqrt( r_one - sw2 )
      ez=ee/(sw*cw)
      ey=ee*(sw/cw)
c
      gal(1) =  ee
      gal(2) =  ee
      gau(1) = -ee*r_two/r_three
      gau(2) = -ee*r_two/r_three
      gad(1) =  ee   /r_three
      gad(2) =  ee   /r_three
      gwf(1) = -ee/sqrt(r_two*sw2)
      gwf(2) =  r_zero
      gzn(1) = -ez*  r_half
      gzn(2) =  r_zero
      gzl(1) = -ez*(-r_half+sw2)
      gzl(2) = -ey
      gzu(1) = -ez*( r_half-sw2*r_two/r_three)
      gzu(2) =  ey*          r_two/r_three
      gzd(1) = -ez*(-r_half+sw2   /r_three)
      gzd(2) = -ey             /r_three
      g1(1)  =  r_one
      g1(2)  =  r_one
c
      return
      end
