c
c **********************************************************************
c
      subroutine coup1x(sw2, gw,gwwa,gwwz)
c
c this subroutine sets up the coupling constants of the gauge bosons in 
c the standard model.                                                   
c                                                                       
c input:                                                                
c       real    sw2            : square of sine of the weak angle       
c                                                                       
c output:                                                               
c       real    gw             : weak coupling constant                 
c       real    gwwa           : dimensionless coupling of w-,w+,a      
c       real    gwwz           : dimensionless coupling of w-,w+,z      
c
      real*8    sw2,gw,gwwa,gwwz,alpha,fourpi,ee,sw,cw
c
      real*8 r_one, r_four, r_ote, r_pi, r_ialph
      parameter( r_one=1.0d0, r_four=4.0d0, r_ote=128.878d0 )
      parameter( r_pi=3.14159265358979323846d0, r_ialph=137.0359895d0 )
c
      alpha = r_one / r_ote
c      alpha = r_one / r_ialph
      fourpi = r_four * r_pi
      ee=sqrt( alpha * fourpi )
      sw=sqrt( sw2 )
      cw=sqrt( r_one - sw2 )
c
      gw    =  ee/sw
      gwwa  =  ee
      gwwz  =  ee*cw/sw
c
      return
      end
