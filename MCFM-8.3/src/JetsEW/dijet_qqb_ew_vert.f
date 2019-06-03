      subroutine dijet_vt1(vrt,mv,ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: vrt,mv,ss,tt,uu,sxs,f1,vos,cs

      vos = mv**2/ss
C --- where 'cs' is color structure
      cs = 2._dp
      vrt = cs*f1(vos)*sxs(ss,tt,uu)

      end subroutine dijet_vt1



C --- returns gluon vertex correction
c      subroutine dijet_vt2(vrt,ss,tt,uu,ep)

c      implicit none


c      real(dp):: vrt,mv,ss,tt,uu,fg,cs1,cs2,sxt
c      integer ep


C --- the color structure is according to the order of diagrams, which is 
C --- strong first and weak second
c      cs1 = -2._dp/3._dp
c      cs2 = 16._dp/3._dp
c      vrt = (cs1*fg(ss,ep) + cs2*fg(tt,ep))*sxt(ss,tt,uu)


c      end subroutine dijet_vt2



c      real(dp):: function fg(x,ep)


c      implicit none

c      include 'constants.f'
c      include 'scale.f'


c      real(dp):: x,xI2,xI3
c      integer ep


c      fg = - 3._dp*xI2(x,0._dp,0._dp,musq,ep) 
c     .     - 2._dp*x*xI3(0._dp,x,0._dp,0._dp,0._dp,0._dp,musq,ep)


c      write(6,*) "fg --- 1:", fg, ep, musq, x
c      fg = - 3._dp*real(ep+2,dp)*(-real(ep,dp) +
c     .     real(ep+1,dp)*(log(musq/abs(x))/2._dp + 1._dp))
c     .     - 4._dp*(real(ep+1,dp)*real(ep,dp)/2._dp 
c     .     - log(musq/abs(x))*real(ep+2,dp)*real(ep,dp)
c     .     + real(ep+2,dp)*real(ep+1,dp)*(1._dp/4._dp*log(musq/abs(x))**2 
c     .     - pi**2/4._dp))
c      write(6,*) "fg --- 2:", fg, ep, musq, x

c      pause


c      end function fg



       subroutine dijet_vt2(vrt,ss,tt,uu)
C --- if it is still strong amp first and weak second
C --- include charge renormalization (MS-bar scheme)
       implicit none
       include 'types.f'
       include 'nf.f'
       include 'constants.f'
       include 'epinv.f'
       include 'scale.f'
       include 'masses.f'
       real(dp):: vrt,ss,tt,uu,frm1(2),frm2(2),sig,cs_slf,sxt,
     .      ren_g,duv,cs1,cs2

C --- ren_g is the strong coupling renormalization (the part propt to gs_bare)
C --- top is decoupled from the strong running coupling, nf = 5
c       duv = (real(nf,dp)/3._dp - 11._dp/2._dp)*(epinv + log(musq/mt**2))
c       ren_g = duv/(8._dp*pi**2) + epinv/(24._dp*pi**2)
       duv = (real(nf,dp)/3._dp - 11._dp/2._dp)*epinv
       ren_g = 2._dp*(duv + (epinv + log(musq/mt**2))/3._dp)
c       ren_g = (2._dp*real(nf,dp)/3._dp - 13._dp/3._dp - 2._dp)*epinv

! for checking the single poles
c       ren_g = (2._dp/3._dp*real(nf+1,dp) - 11._dp)*epinv

       call qcd_vert(frm1,ss)
       call qcd_vert(frm2,tt)
       call qcd_slf(sig,ss)

c       cs1 = -2._dp/3._dp
c       cs2 = 16._dp/3._dp
c       write(6,*) "vert: ", frm1(1) + frm2(2), frm1(2) + frm2(1), 
c     .      2._dp*(cs1*log(abs(ss)/musq) + cs2*log(abs(tt)/musq))
c       write(6,*) "self: ", (ren_g + sig)/1d3
c       pause


C --- cs_slf is color structure to self-energy contribution
C --- cs_slf = Tr(T_a*T_b)*delta(a,b)
       cs_slf = 4._dp
       vrt = (2._dp*frm1(1) + 2._dp*frm2(2) + cs_slf*sig
     .      + cs_slf*ren_g)*sxt(ss,tt,uu)
! for checking the double poles
c       vrt = (2._dp*frm1(1) + 2._dp*frm2(2))*sxt(ss,tt,uu)

c       write(6,*) "frm1(1), frm2(2): ", frm1(1), frm2(2)

       vrt = 2._dp*vrt

       end subroutine dijet_vt2




C --- QCD correction to vertex qqg and qqV(V = W,Z/gamma); 
C --- form factor w/ full color structure as that in matrix element
C --- the order of the leading amp is amp(strong) X amp(weak)
      subroutine qcd_vert(frm,x)
C --- x-channel; returns form factor 'frm' which has two components
C --- frm(1) is that to vertex qqg
C --- frm(2) is that to vertex qqv
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'cplx.h'
      real(dp):: frm(2),x,lam1,lam2,bb0,cc0,cs1,cs2

C --- bubble and triangle int in terms of logarithms (real)
C --- arguments: bb0(x;0,0), cc0(0,x,0;0,0,0)
      bb0 = epinv + log(musq/abs(x)) + 2._dp
C --- cc0 = x*cc0
      cc0 = epinv*epinv2 + epinv*real(log(cplx1(-musq/x)),dp) 
     .     + real(log(cplx1(-musq/x))**2/2._dp,dp)
c     .     + log(musq/abs(x))**2/2._dp
c      if(x > 0) then
c         cc0 = cc0 - pi**2/2._dp
c      end if

! for checking the single poles
c      bb0 = 1._dp*0
c      cc0 = log(musq/abs(x))

c      epinv = 0._dp
c      epinv2 = 0._dp

! for checking poles
c      bb0 = epinv
c      cc0 = epinv2**2 + epinv*log(musq/abs(x))
c      cc0 = 0._dp

C --- form factors lam1 and lam2 to the vertex qqg
      lam1 = -3._dp*bb0 - 2._dp*cc0
C --- finite part from eps*(1/eps)
     .     - 2._dp
      lam2 = -bb0
! for checking the double poles
c      lam1 = -2._dp*epinv2**2
c      lam2 = 0._dp

      cs1 = -2._dp/3._dp
      cs2 = 16._dp/3._dp

      frm(1) = cs1*lam1 + (cs2-cs1)*lam2
      frm(2) = cs2*lam1

      end subroutine qcd_vert



C --- transverse scalar self-energy, including gluon, ghost and fermion loops
      subroutine qcd_slf(sig,ss)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: sig,ss,xI1,xI2

c      epinv = 0._dp

      sig = 1._dp/3._dp*(
c     .     epinv*(
c     .     + 4._dp*xI1(mt**2,musq,-1) 
c     .     + (15._dp - 2._dp*real(nf,dp))*ss*xI2(ss,0._dp,0._dp,musq,-1)
c     .     - 2._dp*(2._dp*mt**2 + ss)*xI2(ss,mt**2,mt**2,musq,-1))
     .     + 4._dp*xI1(mt**2,musq,0) 
     .     + (15._dp - 2._dp*real(nf,dp))*ss*xI2(ss,0._dp,0._dp,musq,0)
     .     - 2._dp*(2._dp*mt**2 + ss)*xI2(ss,mt**2,mt**2,musq,0)
     .     )

      sig = sig + (13._dp -2._dp*real(nf,dp))*ss/3._dp*epinv
C --- finite part from ep*(1/ep)
     .     + (12._dp*mt**2 - (23._dp + 2._dp*real(nf,dp))*ss)/18._dp

      if(.false.) then
      write(6,*) "compact form for the single pole in vertex"
      write(6,*) "comparing w/ the poles returned by QCDLoop"
      write(6,*) (+ 4._dp*xI1(mt**2,musq,-1) 
     .     + (15._dp - 2._dp*real(nf,dp))*ss*xI2(ss,0._dp,0._dp,musq,-1)
     .     - 2._dp*(2._dp*mt**2 + ss)*xI2(ss,mt**2,mt**2,musq,-1))/
     .     ((13._dp -2._dp*real(nf,dp))*ss)
      pause
      end if

! for checking the single poles
c      sig = (13._dp -2._dp*real(nf,dp))*ss/3._dp*epinv

      sig = sig/ss

      end subroutine qcd_slf



      subroutine dijet_vt3(vrt,mv,ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: vrt,mv,ss,tt,uu,vos,vot,f1,cs,sxt

      vos = mv**2/ss
      vot = mv**2/tt
      cs = -2._dp/3._dp
      vrt = cs*(f1(vos) + f1(vot))*sxt(ss,tt,uu)
c      vrt = cs*f1(vos)*sxt(ss,tt,uu)

      end subroutine dijet_vt3



C --- vertex corrections to 
C --- (1) u_i+ub_i -> u_j + ub_j; (2) d_i + db_i -> d_j + db_j
      subroutine vertex1(vrt,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
C --- vrt(i,j): indices 'i', 'j' describe initial and final states (quarks)
C --- anti-quarks are adjoint to corresponding quarks
      real(dp):: vrt(nf,nf),vv(2),ss,tt,uu,mz,mw,cplz(2),cplw

      mz = zmass
      mw = wmass

c      cplz(1:2) = l(1:2)**2 + r(1:2)**2
c      cplw = gwsq/4._dp
      cplz(1:2) = (l(1:2)**2 + r(1:2)**2)/2._dp
      cplw = 1._dp/8._dp/xw

      call dijet_vt1(vv(1),mz,ss,tt,uu)
      call dijet_vt1(vv(2),mw,ss,tt,uu)

      vrt = 0._dp
C --- including initial and final vertex corrections (first 2._dp for)
      vrt(1,3) = -2._dp*(cplz(1)*vv(1) + 2._dp*cplw*vv(2))
      vrt(2,4) = -2._dp*(cplz(2)*vv(1) + 2._dp*cplw*vv(2))
      vrt(3,1) = vrt(1,3)
      vrt(4,2) = vrt(2,4)

C --- '2' from 2*Re(dM*M)
c      vrt = 2._dp*vrt/(4._dp*pi)**2
      vrt = vrt/(4._dp*pi)**2
      vrt = gsq**2*esq*vrt*aveqq

      end subroutine vertex1



C --- vertex corrections to 
C --- (1) u_i+ub_i -> d_j + db_j; (2) d_i + db_i -> u_j + ub_j
      subroutine vertex2(vrt,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
!      include 'epinv.f'
!      include 'epinv2.f'
      real(dp):: vrt(nf,nf),ss,tt,uu,v1(nf,nf),vv1(2),
     .     vv2,propw,mz,mw,ccz,ccw
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

!      epinv=0._dp ! DEBUG
      mz = zmass
      mw = wmass
      propw = tt*(tt-mw**2)/((tt-mw**2)**2 + wwidth**2*mw**2)

      call dijet_vt1(vv1(1),mz,ss,tt,uu)
      call dijet_vt1(vv1(2),mw,ss,tt,uu)
      call dijet_vt2(vv2,ss,tt,uu)

      v1 = 0._dp
      ccz = (l(1)**2 + r(1)**2 + l(2)**2 + r(2)**2)/2._dp
      ccw = 1._dp/8._dp/xw
c      v1(1,2) = -2._dp*(ccz*vv1(1) + 2._dp*ccw*vv1(2))
      v1(1,2) = -(ccz*vv1(1) + 4._dp*ccw*vv1(2))
      v1(1,4) = v1(1,2)
      v1(2,1) = v1(1,2)
      v1(2,3) = v1(1,2)
      v1(3,2) = v1(1,2)
      v1(3,4) = v1(1,2)
      v1(4,1) = v1(1,2)
      v1(4,3) = v1(1,2)

c      vv2 = v2(0) + v2(-1)*epinv + v2(-2)*epinv2**2

c      vrt = 0._dp
c      vrt(1,2) = -2._dp*vv2*Vud**2*propw
c      vrt(1,4) = -2._dp*vv2*Vcd**2*propw
c      vrt(2,1) = vrt(1,2)
c      vrt(2,3) = -2._dp*vv2*Vus**2*propw
c      vrt(3,2) = vrt(2,3)
c      vrt(3,4) = -2._dp*vv2*Vcs**2*propw
c      vrt(4,1) = vrt(4,1)
c      vrt(4,3) = vrt(3,4)

      vrt = 0._dp
      vrt(1,2) = +vv2*Vud**2*propw
      vrt(1,4) = +vv2*Vcd**2*propw
      vrt(2,1) = vrt(1,2)
      vrt(2,3) = +vv2*Vus**2*propw
      vrt(3,2) = vrt(2,3)
      vrt(3,4) = +vv2*Vcs**2*propw
      vrt(4,1) = vrt(1,4)
      vrt(4,3) = vrt(3,4)
      vrt = vrt*2._dp*ccw

C --- '2' from 2*Re(dM*M)
      vrt = v1 + vrt
c      vrt = 2._dp*vrt/(4._dp*pi)**2
      vrt = vrt/(4._dp*pi)**2
      vrt = gsq**2*esq*vrt*aveqq

      end subroutine vertex2



C --- vertex corrections to 
C --- (1) u_i+ub_i -> u_i + ub_i; (2) d_i + db_i -> d_i + db_i
      subroutine vertex3(vrt,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
!      include 'epinv.f'
!      include 'epinv2.f'
      real(dp):: vrt(nf,nf),v1(nf,nf),vgg(nf,nf),vmix(nf,nf),
     .     vg1,vg2,vvgz,vvzg,propz_s,propz_t,mz,mw,ss,tt,uu,
     .     ccz(2),ccw,vv1(2),vv2(2),Qf(nf)
      logical flg/.false./

!      epinv=0._dp ! DEBUG
c      cplz(1:2) = (l(1:2)**2 + r(1:2)**2)/2._dp
c      cplw = gwsq/4._dp
c      cplw = 1._dp/8._dp/xw
      mz = zmass
      mw = wmass
      propz_s = ss*(ss-mz**2)/((ss-mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt-mz**2)/((tt-mz**2)**2 + zwidth**2*mz**2)

      call dijet_vt1(vv1(1),mz,ss,tt,uu)
      call dijet_vt1(vv1(2),mw,ss,tt,uu)
      call dijet_vt1(vv2(1),mz,tt,ss,uu)
      call dijet_vt1(vv2(2),mw,tt,ss,uu)
      call dijet_vt3(vg1,mz,ss,tt,uu)
      call dijet_vt3(vg2,mw,ss,tt,uu)

      ccz(1:2) = (l(1:2)**2 + r(1:2)**2)/2._dp
      ccw = 1._dp/8._dp/xw

C--- flag determines whether photon mediated amplitude is included

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)

      v1(1,1) = -2._dp*(ccz(1)*(vv1(1) + vv2(1)) 
     .     + 2._dp*ccw*(vv1(2) + vv2(2)))
      v1(2,2) = -2._dp*(ccz(2)*(vv1(1) + vv2(1)) 
     .     + 2._dp*ccw*(vv1(2) + vv2(2)))
      v1(3,3) = v1(1,1)
      v1(4,4) = v1(2,2)

!      v1 = 2._dp*v1/(4._dp*pi)**2
      v1 = v1/(4._dp*pi)**2

      vgg = 0._dp
C --- the last 2._dp due to the crossing ss <-> tt
      vgg(1,1) = -2._dp*(ccz(1)*vg1 + 2._dp*ccw*vg2)!*2._dp
      vgg(2,2) = -2._dp*(ccz(2)*vg1 + 2._dp*ccw*vg2)!*2._dp
      vgg(3,3) = vgg(1,1)
      vgg(4,4) = vgg(2,2)

c      vgg = 2._dp*vgg/(4._dp*pi)**2
      vgg = vgg/(4._dp*pi)**2

c      call dijet_vt2(vgz(0),ss,tt,uu,0)
c      call dijet_vt2(vgz(-1),ss,tt,uu,-1)
c      call dijet_vt2(vgz(-2),ss,tt,uu,-2)
c      call dijet_vt2(vzg(0),tt,ss,uu,0)
c      call dijet_vt2(vzg(-1),tt,ss,uu,-1)
c      call dijet_vt2(vzg(-2),tt,ss,uu,-2)

c      vvgz = vgz(0) + vgz(-1)*epinv + vgz(-2)*epinv2**2
c      vvzg = vzg(0) + vzg(-1)*epinv + vzg(-2)*epinv2**2

C --- sxs channel vanishes due to Tr(t_a) = 0
      call dijet_vt2(vvgz,ss,tt,uu)
      call dijet_vt2(vvzg,tt,ss,uu)

      vmix = 0._dp
c      vmix(1,1) = +2._dp*((vvgz*propz_t + vvzg*propz_s)*ccz(1) 
c     .     + Qf(1)**2*(vvgz + vvzg))
c      vmix(2,2) = +2._dp*((vvgz*propz_t + vvzg*propz_s)*ccz(2)
c     .     + Qf(2)**2*(vvgz + vvzg))
      vmix(1,1) = +((vvgz*propz_t + vvzg*propz_s)*ccz(1) 
     .     + Qf(1)**2*(vvgz + vvzg))
      vmix(2,2) = +((vvgz*propz_t + vvzg*propz_s)*ccz(2)
     .     + Qf(2)**2*(vvgz + vvzg))
      vmix(3,3) = vmix(1,1)
      vmix(4,4) = vmix(2,2)

c      vmix = 2._dp*vmix/(4._dp*pi)**2
C --- no factor of 2 (included in routine dijet_vt2 already)
      vmix = vmix/(4._dp*pi)**2

      vrt = v1 + vgg + vmix
      vrt = gsq**2*esq*vrt*aveqq

      end subroutine vertex3

      
