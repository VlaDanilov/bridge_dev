      subroutine qqb_z_ew_exact(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'first.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'scale.f'
      include 'zprods_decl.f'
      include 'cplx.h'
      include 'delr.f'
      real(dp)::msq(fn:nf,fn:nf),p(mxpart,4),ss,tt,uu,bx(nf,8,2),
     .     bxtot(fn:nf,fn:nf),born(fn:nf,fn:nf),vert(fn:nf,fn:nf),
     .     gvq(nf),gaq(nf),f1,fac,bfac,propa,mz,mw,gvl,gal,xI1,xI2,xI3,
     .     fbw,aemmz,self(fn:nf,fn:nf),bbrn(fn:nf,fn:nf),ggw,sw2,
     .     cw2,nng(2),nnz(2),llg(2),llz(2),ulg(2),ulz(2),dlg(2),dlz(2),
     .     ctzz,ctaa,ctaz,ctww,p31(2),p32(2),p34(2),p35(2),
     .     drho,xt,gmu,da,cpl(2,2),ccpl(2,0:1),bbx(8),rpropz
      complex(dp)::qqb,qbq,BLL(nf,4),BRR(nf,4),BLR(nf,4),BRL(nf,4),
     .     sfLL(fn:nf,fn:nf),sfRR(fn:nf,fn:nf),sfLR(fn:nf,fn:nf),
     .     sfRL(fn:nf,fn:nf),propz,ffgLL(fn:nf,fn:nf,2),
     .     ffgRR(fn:nf,fn:nf,2),ffgLR(fn:nf,fn:nf,2),
     .     ffgRL(fn:nf,fn:nf,2),ffzLL(fn:nf,fn:nf,2),
     .     ffzRR(fn:nf,fn:nf,2),ffzLR(fn:nf,fn:nf,2),cmwos,cmzos,
     .     ffzRL(fn:nf,fn:nf,2),lam2,lam3,fwgq(nf),fwgl,fwzq(nf),fwzl,
     .     slfzz,slfaa,slfaz,slfww,s1,s2,ds1,ds2,qlI2,
     .     cpropa,cpropz,cpropaz,cpropza,sigmz,dsigmz
      integer j,k,ep,ier
      common/em/aemmz

      ep = 0
      fac = aemmz/four/pi
      mz = zmass
      mw = wmass
      sw2 = xw
      cw2 = one - sw2

c      call getdoreenPSpoint(p)

      call spinoru(4,p,za,zb)

      qqb = za(2,3)*zb(4,1)
      qbq = za(1,3)*zb(4,2)

      call dotem(4,p,s)
      ss = s(1,2)
      tt = s(1,3)
      uu = s(2,3)


      bfac = 4._dp*esq**2*xn*aveqq

      propz = 1._dp/cplx2((ss-zmass**2),zmass*zwidth)
      propa = 1._dp/ss
      rpropz = 1._dp/(ss-mz**2)
      if (first) then
        first=.false. 
        call delta_r
!        write(6,*) 'delta r is: ', dr
      endif


!      fbw = -1._dp + (2._dp*mt**2 - 2._dp*mw**2 - 3._dp*ss)/ss*
!     .     ((xI1(mw**2, musq, ep) - xI1(mt**2, musq, ep))/(mt**2 
!     .     - mw**2) + xI2(ss, mt**2, mt**2, musq, ep)) 
!     .     - 2._dp/ss*((mt**2 - mw**2)**2 + ss*(ss - mt**2 
!     .     + 2._dp*mw**2))*xI3(0._dp, 0._dp, ss, mt**2, mw**2, mt**2, musq, 
!     .     ep) - (mt**2 - mw**2)*db0(0._dp, mt**2, mw**2)

!      fbw = -4._dp*ggw**2*fbw
!      fbw = 0._dp

      cmwos = cplx2(mw**2/ss,-1.e-12_dp)
      cmzos = cplx2(mz**2/ss,-1.e-12_dp)
!      cmwos = mw**2/ss
!      cmzos = mz**2/ss

C -- couplings related to A/Z->WW->ffb, so far I haven't generalized the 
C -- relation between quantum numbers. It must be connected to the respective 
C -- quantum numbers somehow. But for now, I write the couplings in terms of 
C -- the mixing angles. Let's call the coupling ll, nn for lepton and neutrino 
C -- and ul, dl for up and down-quarks. Each coupling has two components, 1st 
C -- couples to lam2 and 2nd to lam3
C -- Reference: Fortschritte der Physik 38 (1990) 3, 165-260
C -- the renormalized vertex expression can be found in Appendix C

      nng = 0._dp
      llg(1) = 0._dp
      llg(2) = - 3._dp/4._dp/sw2
      nnz(1) = + (2._dp*sw2 - 1._dp)/8._dp/sw2/sqrt(sw2*cw2)
      nnz(2) = + 3._dp/4._dp*cw2/sw2/sqrt(sw2*cw2)
      llz(1) = + 1._dp/8._dp/sw2/sqrt(sw2*cw2)
      llz(2) = - 3._dp/4._dp/sw2*sqrt(cw2/sw2)
C -- for process 31 - 35, no 33
      p31(1) = 1._dp + q1
      p31(2) = le - l1
      p32(1) = -q1
      p32(2) = ln*sqrt(3._dp) - l1
      p34(1) = Q(1)*sqrt(3._dp*xn) - q1
      p34(2) = l(1)*sqrt(3._dp*xn) - l1
      p35(1) = Q(2)*sqrt(2._dp*xn) - q1
      p35(2) = l(2)*sqrt(2._dp*xn) - l1

!      llg = -q1*llg
!      llz = (ln*sqrt(3._dp) - l1)*llz + (le - l1)*nnz*sqrt(3._dp)
!      llz = llz/(ln*sqrt(3._dp) + le - 2._dp*l1)

      ulg(1) = - 1._dp/12._dp/sw2
      ulg(2) = + 3._dp/4._dp/sw2
      dlg(1) = + 1._dp/6._dp/sw2
      dlg(2) = - 3._dp/4._dp/sw2
      ulz(1) = - (1._dp - 2._dp/3._dp*sw2)/8._dp/sw2/sqrt(sw2*cw2)
      ulz(2) = + 3._dp/4._dp/sw2*sqrt(cw2/sw2)
      dlz(1) = + (1._dp - 4._dp/3._dp*sw2)/8._dp/sw2/sqrt(sw2*cw2)
      dlz(2) = - 3._dp/4._dp/sw2*sqrt(cw2/sw2)

C -- select the wanted process for the final states
      llg = + p32(1)*p34(1)*p35(1)*llg 
     .     + p31(1)*p32(1)*p35(1)*ulg*sqrt(2._dp*xn)
     .     + p31(1)*p32(1)*p34(1)*dlg*sqrt(3._dp*xn)
      llg = llg/(p32(1)*p34(1)*p35(1) + p31(1)*p32(1)*p35(1) 
     .     + p31(1)*p32(1)*p34(1))
      llz =  + p32(2)*p34(2)*p35(2)*llz + p31(2)*p34(2)*p35(2)*nnz
     .     + p31(2)*p32(2)*p35(2)*ulg*sqrt(2._dp*xn)
     .     + p31(2)*p32(2)*p34(2)*dlg*sqrt(3._dp*xn)
      llz = llz/(p32(2)*p34(2)*p35(2) + p31(2)*p34(2)*p35(2) 
     .     + p31(2)*p32(2)*p35(2) + p31(2)*p32(2)*p34(2))
C -- for process 34 & 35
!      llg = (Q(1)*sqrt(3._dp*xn) - q1)*ulg*sqrt(2._dp*xn)
!     .     + (Q(2)*sqrt(2._dp*xn) - q1)*dlg*sqrt(3._dp*xn)
!      llg = llg/(Q(1)*sqrt(3._dp*xn) + Q(2)*sqrt(2._dp*xn) - 2._dp*q1)
!      llz = (l(1)*sqrt(3._dp*xn) - l1)*ulz*sqrt(2._dp*xn) 
!     .     + (l(2)*sqrt(2._dp*xn) - l1)*dlz*sqrt(3._dp*xn)
!      llz = llz/(l(1)*sqrt(3._dp*xn) + l(2)*sqrt(2._dp*xn) - 2._dp*l1)

C -- process 33 is postponed for no b quark is considered for vertex for now

      fwgq(1) = dlg(1)*lam2(cmwos) + dlg(2)*lam3(cmwos)
      fwgq(2) = ulg(1)*lam2(cmwos) + ulg(2)*lam3(cmwos)
      fwgq(3) = fwgq(1)
      fwgq(4) = fwgq(2)
      fwgq(5) = 0._dp
      fwgq(5) = fwgq(1)
      fwzq(1) = dlz(1)*lam2(cmwos) + dlz(2)*lam3(cmwos)
      fwzq(2) = ulz(1)*lam2(cmwos) + ulz(2)*lam3(cmwos)
      fwzq(3) = fwzq(1)
      fwzq(4) = fwzq(2)
      fwzq(5) = 0._dp
      fwzq(5) = fwzq(1)
      fwgl = llg(1)*lam2(cmwos) + llg(2)*lam3(cmwos)
      fwzl = llz(1)*lam2(cmwos) + llz(2)*lam3(cmwos)

c      ss = 90164.6136_dp
c      tt = -10550.1699_dp
c      uu = -79614.4437_dp
      call z_ew_box(bx(:,:,1),ss,tt,uu)
      call z_ew_box(bx(:,:,2),ss,uu,tt)

c      write(6,*) bx(2,4,1) + bx(2,8,1)
c      write(6,*) bx(2,4,2) + bx(2,8,2)
c      write(6,*) bx(1,3,1) + bx(1,7,1)
c      write(6,*) bx(1,3,2) + bx(1,7,2)
c      stop

      bx(2,3,:) = 0._dp
      bx(4,3,:) = 0._dp
      bx(1,4,:) = 0._dp
      bx(3,4,:) = 0._dp
      bx(5,4,:) = 0._dp
      bx(2,7,:) = 0._dp
      bx(4,7,:) = 0._dp
      bx(1,8,:) = 0._dp
      bx(3,8,:) = 0._dp
      bx(5,8,:) = 0._dp


      bxtot = 0._dp
      do j = 1,8
         do k = 1,nf
            bxtot(-k,k) = bxtot(-k,k) + bx(k,j,2)
            bxtot(k,-k) = bxtot(k,-k) + bx(k,j,1)
         end do
      end do


c      ss = 73.9645607_dp**2
      call znen(slfaa,aemmz,sw2,ss,ep,'AA')
      call znen(slfaz,aemmz,sw2,ss,ep,'AZ')
      call znen(slfzz,aemmz,sw2,ss,ep,'ZZ')
      call znenteghaon(ctaa,aemmz,sw2,ss,ep,'AA')
      call znenteghaon(ctaz,aemmz,sw2,ss,ep,'AZ')
      call znenteghaon(ctzz,aemmz,sw2,ss,ep,'ZZ')
      call dznenzz(dsigmz,aemmz,sw2,mz**2,ep)
      call znen(sigmz,aemmz,sw2,mz**2,ep,'ZZ')

      if(abs(ss-mz**2) < 1.e-3_dp) then
         slfzz = dsigmz*(ss-mz**2) + sigmz
      end if




C -- BLL(:,4) etc. denotes the born amplitude w/o mediated gauge boson 
C -- propagator, may add the prop later according to the types of couplings
C -- BXX(:,1) is the amplitude with photon-photon coupling
C -- BXX(:,2) w/ Z-Z coupling
C -- BXX(:,3) w/ A-Z coupling
C -- BXX(:,4) w/ Z-A coupling
C -- ffgLL(nf,:) etc. have 2 components, (1) denotes the initial vertex and 
C -- (2) the final
   
      msq = 0._dp

      do j = fn, nf
         k = -j
         if(j == 0) then
            msq(j,k) = 0._dp
         else if(j > 0) then

            BLL(j,1) = Q(j)*q1*qqb*propa
            BRR(j,1) = Q(j)*q1*qqb*propa
            BLR(j,1) = Q(j)*q1*qbq*propa
            BRL(j,1) = Q(j)*q1*qbq*propa
            BLL(j,2) = L(j)*l1*qqb*propz
            BRR(j,2) = R(j)*r1*qqb*propz
            BLR(j,2) = L(j)*r1*qbq*propz
            BRL(j,2) = R(j)*l1*qbq*propz
            BLL(j,3) = Q(j)*l1*qqb
            BRR(j,3) = Q(j)*r1*qqb
            BLR(j,3) = Q(j)*r1*qbq
            BRL(j,3) = Q(j)*l1*qbq
            BLL(j,4) = L(j)*q1*qqb
            BRR(j,4) = R(j)*q1*qqb
            BLR(j,4) = L(j)*q1*qbq
            BRL(j,4) = R(j)*q1*qbq

            sfLL(j,k) = - real(slfaa + ctaa,dp)*BLL(j,1)*propa 
     .           - real(slfzz + ctzz,dp)*BLL(j,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BLL(j,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BLL(j,4)*propa*propz
            sfRR(j,k) = - real(slfaa + ctaa,dp)*BRR(j,1)*propa 
     .           - real(slfzz + ctzz,dp)*BRR(j,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BRR(j,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BRR(j,4)*propa*propz
            sfLR(j,k) = - real(slfaa + ctaa,dp)*BLR(j,1)*propa 
     .           - real(slfzz + ctzz,dp)*BLR(j,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BLR(j,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BLR(j,4)*propa*propz
            sfRL(j,k) = - real(slfaa + ctaa,dp)*BRL(j,1)*propa 
     .           - real(slfzz + ctzz,dp)*BRL(j,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BRL(j,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BRL(j,4)*propa*propz


            ffgLL(j,k,1) = +fac*real(Q(j)*L(j)**2*lam2(cmzos) 
     .           + 2._dp*fwgq(j),dp)*q1*qqb*propa
            ffgLL(j,k,2) = +fac*real(q1*l1**2*lam2(cmzos)
     .           + 2._dp*fwgl,dp)*Q(j)*qqb*propa
            ffgRR(j,k,1) = +fac*real(Q(j)*R(j)**2*lam2(cmzos),dp)*
     .           q1*qqb*propa
            ffgRR(j,k,2) = +fac*real(q1*r1**2*lam2(cmzos),dp)*Q(j)*
     .           qqb*propa
            ffgLR(j,k,1) = +fac*real(Q(j)*L(j)**2*lam2(cmzos) 
     .           + 2._dp*fwgq(j),dp)*q1*qbq*propa
            ffgLR(j,k,2) = +fac*real(q1*r1**2*lam2(cmzos),dp)*Q(j)*
     .           qbq*propa
            ffgRL(j,k,1) = +fac*real(Q(j)*R(j)**2*lam2(cmzos),dp)*
     .           q1*qbq*propa
            ffgRL(j,k,2) = +fac*real(q1*l1**2*lam2(cmzos)
     .           +2._dp*fwgl,dp)*Q(j)*qbq*propa


            ffzLL(j,k,1) = +fac*real(L(j)**3*lam2(cmzos) 
     .           + 2._dp*fwzq(j),dp)*l1*qqb*propz
            ffzLL(j,k,2) = +fac*real(l1**3*lam2(cmzos)
     .           + 2._dp*fwzl,dp)*L(j)*qqb*propz
            ffzRR(j,k,1) = +fac*real(R(j)**3*lam2(cmzos),dp)*r1*qqb*propz
            ffzRR(j,k,2) = +fac*real(r1**3*lam2(cmzos),dp)*R(j)*qqb*propz
            ffzLR(j,k,1) = +fac*real(L(j)**3*lam2(cmzos) 
     .           + 2._dp*fwzq(j),dp)*r1*qbq*propz
            ffzLR(j,k,2) = +fac*real(r1**3*lam2(cmzos),dp)*L(j)*qbq*propz
            ffzRL(j,k,1) = +fac*real(R(j)**3*lam2(cmzos),dp)*l1*qbq*propz
            ffzRL(j,k,2) = +fac*real(l1**3*lam2(cmzos)
     .           + 2._dp*fwzl,dp)*R(j)*qbq*propz


            self(j,k) = + real(sfLL(j,k)*conjg(BLL(j,1) + BLL(j,2)),dp) 
     .           + real(sfRR(j,k)*conjg(BRR(j,1) + BRR(j,2)),dp)
     .           + real(sfLR(j,k)*conjg(BLR(j,1) + BLR(j,2)),dp)
     .           + real(sfRL(j,k)*conjg(BRL(j,1) + BRL(j,2)),dp)

            self(j,k) = 2._dp*bfac*self(j,k)

!            born(j,k) = + cdabs(BLL(j,1) + BLL(j,2))**2
!     .           + cdabs(BRR(j,1) + BRR(j,2))**2
!     .           + cdabs(BLR(j,1) + BLR(j,2))**2
!     .           + cdabs(BRL(j,1) + BRL(j,2))**2
!            born(j,k) = bfac*born(j,k)

C -- initially set k < nf for excluding the initial b-quark
            if(j <= nf) then
!               vert(j,k) = -fac*born(j,k)*(
!     .              + (gvq(j)**2+gaq(j)**2+gvl**2+gal**2)*f1(mz**2/ss)
!     .              + 2._dp*ggw**2*2._dp*f1(mw**2/ss))
               vert(j,k) = + real((ffgLL(j,k,1) + ffgLL(j,k,2) 
     .              + ffzLL(j,k,1) + ffzLL(j,k,2))*conjg(BLL(j,1) 
     .              + BLL(j,2)),dp)
     .              + real((ffgRR(j,k,1) + ffgRR(j,k,2) + ffzRR(j,k,1) 
     .              + ffzRR(j,k,2))*conjg(BRR(j,1) + BRR(j,2)),dp)
     .              + real((ffgLR(j,k,1) + ffgLR(j,k,2) + ffzLR(j,k,1) 
     .              + ffzLR(j,k,2))*conjg(BLR(j,1) + BLR(j,2)),dp)
     .              + real((ffgRL(j,k,1) + ffgRL(j,k,2) + ffzRL(j,k,1) 
     .              + ffzRL(j,k,2))*conjg(BRL(j,1) + BRL(j,2)),dp)
               vert(j,k) = 2._dp*bfac*vert(j,k)
            else
!               vert(j,k) = -fac*born(j,k)*(
!     .              + (gvq(j)**2+gaq(j)**2+gvl**2+gal**2)*f1(mz**2/ss)
!     .              + 2._dp*ggw**2*f1(mw**2/ss) + fbw)
               vert(j,k) = 0._dp
            end if

            msq(j,k) = vert(j,k) + self(j,k) + bxtot(j,k)
         else if(j < 0) then

            BLL(k,1) = Q(k)*q1*qbq*propa
            BRR(k,1) = Q(k)*q1*qbq*propa
            BLR(k,1) = Q(k)*q1*qqb*propa
            BRL(k,1) = Q(k)*q1*qqb*propa
            BLL(k,2) = L(k)*l1*qbq*propz
            BRR(k,2) = R(k)*r1*qbq*propz
            BLR(k,2) = L(k)*r1*qqb*propz
            BRL(k,2) = R(k)*l1*qqb*propz
            BLL(k,3) = Q(k)*l1*qbq
            BRR(k,3) = Q(k)*r1*qbq
            BLR(k,3) = Q(k)*r1*qqb
            BRL(k,3) = Q(k)*l1*qqb
            BLL(k,4) = L(k)*q1*qbq
            BRR(k,4) = R(k)*q1*qbq
            BLR(k,4) = L(k)*q1*qqb
            BRL(k,4) = R(k)*q1*qqb

            sfLL(j,k) = - real(slfaa + ctaa,dp)*BLL(k,1)*propa 
     .           - real(slfzz + ctzz,dp)*BLL(k,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BLL(k,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BLL(k,4)*propa*propz
            sfRR(j,k) = - real(slfaa + ctaa,dp)*BRR(k,1)*propa 
     .           - real(slfzz + ctzz,dp)*BRR(k,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BRR(k,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BRR(k,4)*propa*propz
            sfLR(j,k) = - real(slfaa + ctaa,dp)*BLR(k,1)*propa 
     .           - real(slfzz + ctzz,dp)*BLR(k,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BLR(k,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BLR(k,4)*propa*propz
            sfRL(j,k) = - real(slfaa + ctaa,dp)*BRL(k,1)*propa 
     .           - real(slfzz + ctzz,dp)*BRL(k,2)*rpropz 
     .           + real(slfaz + ctaz,dp)*BRL(k,3)*propa*propz 
     .           + real(slfaz + ctaz,dp)*BRL(k,4)*propa*propz



            ffgLL(j,k,1) = +fac*real(Q(k)*L(k)**2*lam2(cmzos) 
     .           + 2._dp*fwgq(k),dp)*q1*qbq*propa
            ffgLL(j,k,2) = +fac*real(q1*l1**2*lam2(cmzos)
     .           +2._dp*fwgl,dp)*Q(k)*qbq*propa
            ffgRR(j,k,1) = +fac*real(Q(k)*R(k)**2*lam2(cmzos),dp)*
     .           q1*qbq*propa
            ffgRR(j,k,2) = +fac*real(q1*r1**2*lam2(cmzos),dp)*Q(k)*
     .           qbq*propa
            ffgLR(j,k,1) = +fac*real(Q(k)*L(k)**2*lam2(cmzos) 
     .           + 2._dp*fwgq(k),dp)*q1*qqb*propa
            ffgLR(j,k,2) = +fac*real(q1*r1**2*lam2(cmzos),dp)*Q(k)*
     .           qqb*propa
            ffgRL(j,k,1) = +fac*real(Q(k)*R(k)**2*lam2(cmzos),dp)*
     .           q1*qqb*propa
            ffgRL(j,k,2) = +fac*real(q1*l1**2*lam2(cmzos)
     .           + 2._dp*fwgl,dp)*Q(k)*qqb*propa


            ffzLL(j,k,1) = +fac*real(L(k)**3*lam2(cmzos) 
     .           + 2._dp*fwzq(k),dp)*l1*qbq*propz
            ffzLL(j,k,2) = +fac*real(l1**3*lam2(cmzos)
     .           + 2._dp*fwzl,dp)*L(k)*qbq*propz
            ffzRR(j,k,1) = +fac*real(R(k)**3*lam2(cmzos),dp)*r1*qbq*propz
            ffzRR(j,k,2) = +fac*real(r1**3*lam2(cmzos),dp)*R(k)*qbq*propz
            ffzLR(j,k,1) = +fac*real(L(k)**3*lam2(cmzos) 
     .           + 2._dp*fwzq(k),dp)*r1*qqb*propz
            ffzLR(j,k,2) = +fac*real(r1**3*lam2(cmzos),dp)*L(k)*qqb*propz
            ffzRL(j,k,1) = +fac*real(R(k)**3*lam2(cmzos),dp)*l1*qqb*propz
            ffzRL(j,k,2) = +fac*real(l1**3*lam2(cmzos)
     .           + 2._dp*fwzl,dp)*R(k)*qqb*propz


            self(j,k) = + real(sfLL(j,k)*conjg(BLL(k,1) + BLL(k,2)),dp) 
     .           + real(sfRR(j,k)*conjg(BRR(k,1) + BRR(k,2)),dp)
     .           + real(sfLR(j,k)*conjg(BLR(k,1) + BLR(k,2)),dp)
     .           + real(sfRL(j,k)*conjg(BRL(k,1) + BRL(k,2)),dp)

            self(j,k) = 2._dp*bfac*self(j,k)


!            born(j,k) = + cdabs(BLL(k,1) + BLL(k,2))**2
!     .           + cdabs(BRR(k,1) + BRR(k,2))**2
!     .           + cdabs(BLR(k,1) + BLR(k,2))**2
!     .           + cdabs(BRL(k,1) + BRL(k,2))**2

!            born(j,k) = bfac*born(j,k)

C -- initially set k < nf for excluding the initial b-quark
            if(k <= nf) then
!               vert(j,k) = -fac*born(j,k)*(
!     .              + (gvq(k)**2+gaq(k)**2+gvl**2+gal**2)*f1(mz**2/ss)
!     .              + 2._dp*ggw**2*2._dp*f1(mw**2/ss))
               vert(j,k) = + real((ffgLL(j,k,1) + ffgLL(j,k,2) 
     .              + ffzLL(j,k,1) + ffzLL(j,k,2))*conjg(BLL(k,1) 
     .              + BLL(k,2)),dp)
     .              + real((ffgRR(j,k,1) + ffgRR(j,k,2) + ffzRR(j,k,1) 
     .              + ffzRR(j,k,2))*conjg(BRR(k,1) + BRR(k,2)),dp)
     .              + real((ffgLR(j,k,1) + ffgLR(j,k,2) + ffzLR(j,k,1) 
     .              + ffzLR(j,k,2))*conjg(BLR(k,1) + BLR(k,2)),dp)
     .              + real((ffgRL(j,k,1) + ffgRL(j,k,2) + ffzRL(j,k,1) 
     .              + ffzRL(j,k,2))*conjg(BRL(k,1) + BRL(k,2)),dp)
               vert(j,k) = 2._dp*bfac*vert(j,k)
            else
!               vert(j,k) = -fac*born(j,k)*(
!     .              + (gvq(k)**2+gaq(k)**2+gvl**2+gal**2)*f1(mz**2/ss)
!     .              + 2._dp*ggw**2*f1(mw**2/ss) + fbw)
               vert(j,k) = 0._dp
            end if

            msq(j,k) = vert(j,k) + self(j,k) + bxtot(j,k)
         end if
      end do



! my dr = 3.02905165411810687d-2; sigw0 = 81.812060683631486, ct = 69.633881468816412
      call qqb_z(p,bbrn)
c      dr = 0.0304911655_dp
      msq = msq - 2._dp*bbrn*dr

      end subroutine





      subroutine delta_r

      include 'types.f'
      include 'delr.f'
      include 'masses.f'
      include 'constants.f'
      include 'ewcouple.f'

      real(dp) :: ctww,sw2,cw2,mw,aemmz
      complex(dp) :: sigw
      integer :: ep
      common/em/aemmz

      ep = 0
      sw2 = xw
      cw2 = 1._dp - sw2
      mw = wmass

      call sigw0(sigw,aemmz,sw2,ep)
      call znenteghaon(ctww,aemmz,sw2,0.d0,ep,"WW")
c!      sigw = sigw + ctww
      dr = real(sigw + ctww)/mw**2 + aemmz/4.d0/pi/sw2*(6.d0
     .     + (7.d0 - 4.d0*sw2)/2.d0/sw2*log(cw2))
c      print*, 'sigw: ', sigw, ctww
c      print*, aemmz,sw2
c      call znen(sigw,aemmz,sw2,mw**2,ep,'WW')
c      print*, sigw
c      print*, dr
c      stop
      return

      end subroutine delta_r





C -- function lam2 and lam3 are extracted from literature Nucl. Phys. B365, 24
C -- (1991) by W. Beenakker et al. See in (3.13). They are used in Z/gamma-ffb 
C -- vertex when f is light quark(no bottom) or lepton.
      function lam2(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: lam2,x,cli2
!      real(dp)::x,ddilog


      lam2 = -7._dp/2._dp - 2._dp*x - (2._dp*x + 3._dp)*log(-x) 
     .     + 2._dp*(1._dp + x)**2*(cli2(1._dp + 1._dp/x) - pisqo6)

!      lam2 = cmplx(-7._dp/2._dp - 2._dp*x - (2._dp*x + 3._dp)*log(x)
!     .     + 2._dp*(1._dp + x)**2*(log(x)*log((1._dp + x)/x) -
!     .     ddilog(-1._dp/x)), -pi*(3._dp + 2._dp*x 
!     .     - 2._dp*(x + 1._dp)**2*log((1._dp + x)/x)))

      end function lam2




      function lam3(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: lam3,x,y
!      real(dp)::x,y

!      if(x > 0.25_dp) then
!         y = sqrt(4._dp*x - 1._dp)

!         lam3 = cmplx(5._dp/6._dp - 2._dp*x/3._dp + 2._dp/3._dp*(2._dp*x 
!     .        + 1._dp)*y*datan(1._dp/y) 
!     .        - 8._dp/3._dp*x*(x + 2._dp)*(datan(1._dp/y))**2)
!      else
!         lam3 = 0._dp
!      end if

      y = (sqrt(1._dp - 4._dp*x) - 1._dp)/(sqrt(1._dp - 4._dp*x) 
     .     + 1._dp)
      lam3 = +5._dp/6._dp - 2._dp*x/3._dp -(2._dp*x 
     .     + 1._dp)/3._dp*sqrt(1._dp - 4._dp*x)*log(y) 
     .     + 2._dp/3._dp*x*(x + 2._dp)*log(y)**2

      end function lam3





      subroutine getdoreenPSpoint(p)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4)
      
C -- cross checked with Doreen at one phase space point
      p(1,4) = -0.8660254037844385e+03_dp
      p(1,1) = -0.0000000000000000e+00_dp
      p(1,2) = -0.0000000000000000e+00_dp
      p(1,3) = -0.8660254037844387e+03_dp

      p(2,4) = -0.8660254037844387e+03_dp
      p(2,1) = 0.0000000000000000e+00_dp
      p(2,2) = 0.0000000000000000e+00_dp
      p(2,3) = 0.8660254037844387e+03_dp


      p(3,4) = 0.8660254037844385e+03_dp
      p(3,1) = 0.2569122186398718e+03_dp
      p(3,2) =-0.7906945058071030e+03_dp
      p(3,3) = 0.24248362913813273e+03_dp

      p(4,4) = 0.8660254037844387e+03_dp
      p(4,1) =-0.2569122186398718e+03_dp
      p(4,2) = 0.7906945058071030e+03_dp
      p(4,3) = -0.24248362913813273e+03_dp

      return
      end
      





c      subroutine sigwa(sig,alpha,psq,ep)

c      implicit none
c      include 'types.f'
c      include 'constants.f'
c      include 'masses.f'
c      include 'scale.f'

c      complex(dp)::sig,qlI2
c      real(dp)::alpha,psq
c      integer::ep

      
c      sig = (two*wmass**2+five*psq)*qlI2(psq,wmass**2,zip,musq,ep)
c     .     - two*wmass**2*qlI2(zip,wmass**2,wmass**2,musq,ep)
c     .     + psq/three

c      if(psq /= zip) then
c         sig = sig 
c     .        - wmass**4/psq*(qlI2(psq,wmass**2,zip,musq,ep) 
c     .        - qlI2(zip,wmass**2,zip,musq,ep))
c      else
c         sig = sig - wmass**2/two
c      end if

c      sig = -alpha/fourpi*two/three*sig

c      return

c      end subroutine sigwa

      














!      bx(1,3,2) = 0._dp
!      bx(3,3,2) = 0._dp
!      bx(5,3,2) = 0._dp
!      bx(2,4,2) = 0._dp
!      bx(4,4,2) = 0._dp
!      bx(1,7,2) = 0._dp
!      bx(3,7,2) = 0._dp
!      bx(5,7,2) = 0._dp
!      bx(2,8,2) = 0._dp
!      bx(4,8,2) = 0._dp

!      box = 0._dp
!      write(6,*) "uubar: "
!      do j = 1,8
!         write(6,*) "uubar: ", j, bx(2,j,1)
!         write(6,*) "ubaru: ", j, bx(2,j,2)
!      end do
!      write(6,*)
!      write(6,*) "ddbar: "
!      do j = 1,8
!         write(6,*) "ddbar: ", j, bx(1,j,1)
!         write(6,*) "dbard: ", j, bx(1,j,2)
!      end do
!      write(6,*)
!      write(6,*) "bb~: "
!      do j = 1,8
!         write(6,*) "bbbar: ", j, bx(5,j,1)
!         write(6,*) "bbarb: ", j, bx(5,j,2)
!      end do
!      stop     











!      if(.false.) then
!      ep = -1
!      call self_ZZ(slfzz,ss,ep)
!      call self_ct_ZZ(ctzz,ss,ep)
!      call self_AA(slfaa,ss,ep)
!      call self_ct_AA(ctaa,ss,ep)
!      call self_AZ(slfaz,ss,ep)
!      call self_ct_AZ(ctaz,ss,ep)
!
!      print*, "slfzz: ", slfzz
!      print*, "ctzz:  ", ctzz
!      print*, "slfaa: ", slfaa
!      print*, "ctaa:  ", ctaa
!      print*, "slfaz: ", slfaz
!      print*, "ctaz:  ", ctaz
!
!      call self_VV(slfzz,ctzz,ss,ep,'ZZ')
!      call self_VV(slfaa,ctaa,ss,ep,'AA')
!      call self_VV(slfaz,ctaz,ss,ep,'AZ')
!      print*, "slfzz: ", slfzz
!      print*, "ctzz:  ", ctzz
!      print*, "slfaa: ", slfaa
!      print*, "ctaa:  ", ctaa
!      print*, "slfaz: ", slfaz
!      print*, "ctaz:  ", ctaz
!
!      print*, 'pole 1: ', qlI2(ss,mw**2,mw**2,musq,ep)
!      print*, 'pole 2: ', qlI2(0._dp,mw**2,mw**2,musq,ep)
!
!      call self_ZZ(s1,mz**2,ep)
!      call self_ZZ(s2,mz**2+1._dp,ep)
!      call dself_ZZ(ds1,mz**2,ep)
!      call dself_ZZ(ds2,mz**2+1._dp,ep)
!      print*, 's1: ', s1
!      print*, 's2: ', s2
!      print*, 'ds1: ', ds1
!      print*, 'ds2: ', ds2

!      call self_AA(s1,mz**2,ep)
!      call self_AA(s2,mz**2+1._dp,ep)
!      call dself_AA(ds1,mz**2,ep)
!      call dself_AA(ds2,mz**2+1._dp,ep)
!      print*, 's1: ', s1
!      print*, 's2: ', s2
!      print*, 'ds1: ', ds1
!      print*, 'ds2: ', ds2

!      call self_AZ(s1,mz**2,ep)
!      call self_AZ(s2,mz**2+1._dp,ep)
!      call dself_AZ(ds1,mz**2,ep)
!      call dself_AZ(ds2,mz**2+1._dp,ep)
!      print*, 's1: ', s1
!      print*, 's2: ', s2
!      print*, 'ds1: ', ds1
!      print*, 'ds2: ', ds2

!      call self_VV(slfzz,ctzz,ss,-1,'ZZ')
!      call self_VV(slfaa,ctaa,ss,-1,'AA')
!      call self_VV(slfaz,ctaz,ss,-1,'AZ')

!      print*, "slfzz: ", slfzz
!      print*, "ctzz:  ", ctzz
!      print*, "slfaa: ", slfaa
!      print*, "ctaa:  ", ctaa
!      print*, "slfaz: ", slfaz
!      print*, "ctaz:  ", ctaz

!      call self_VV(s1,ds1,mz**2,ep,"ZZ")
!      call self_VV(s2,ds2,mz**2+1._dp,ep,"ZZ")
!      print*, 's1: ', s1
!      print*, 's2: ', s2
!      print*, 'ds1: ', ds1
!      print*, 'ds2: ', ds2
!      call self_VV(s1,ds1,mz**2,ep,"AA")
!      call self_VV(s2,ds2,mz**2+1._dp,ep,"AA")
!      print*, 's1: ', s1
!      print*, 's2: ', s2
!      print*, 'ds1: ', ds1
!      print*, 'ds2: ', ds2
!      call self_VV(s1,ds1,mz**2,ep,"AZ")
!      call self_VV(s2,ds2,mz**2+1._dp,ep,"AZ")
!      print*, 's1: ', s1
!      print*, 's2: ', s2
!      print*, 'ds1: ', ds1
!      print*, 'ds2: ', ds2

!      call self_ZZ(slfzz,ss,ep)
!      call self_ct_ZZ(ctzz,ss,ep)
!      call self_AA(slfaa,ss,ep)
!      call self_ct_AA(ctaa,ss,ep)
!      call self_AZ(slfaz,ss,ep)
!      call self_ct_AZ(ctaz,ss,ep)
!      end if

!      ep = -1
!      call self_VV(slfzz,ctzz,ss,ep,'ZZ')
!      call self_VV(slfaa,ctaa,ss,ep,'AA')
!      call self_VV(slfaz,ctaz,ss,ep,'AZ')
!      print*, "slfzz: ", slfzz
!      print*, "ctzz:  ", ctzz
!      print*, "slfaa: ", slfaa
!      print*, "ctaa:  ", ctaa
!      print*, "slfaz: ", slfaz
!      print*, "ctaz:  ", ctaz
!      print*, "ss: ", ss
!      call znen(slfzz,2._dp*mz**2,ep,'ZZ')
!      call znenteghaon(ctzz,2._dp*mz**2,ep,'ZZ')
!      call znen(slfww,2._dp*mw**2,ep,'WW')
!      call znenteghaon(ctww,2._dp*mw**2,ep,'WW')
!      print*, "renormalized: "
!      print*, "slfzz: ", slfzz
!      print*, "ctzz:  ", ctzz
!      print*, "slfww: ", slfww
!      print*, "ctww:  ", ctww



!      call damz(da,aemmz,ep)
!      print*, "slfzz: ", slfzz
!      print*, "ctzz:  ", ctzz
!      print*, "slfaa: ", slfaa
!      print*, "ctaa:  ", ctaa
!      print*, "slfaz: ", slfaz
!      print*, "ctaz:  ", ctaz
!      stop

!      cpropa = cmplx(1._dp/(ss + slfaa + ctaa 
!     .     - (slfaz + ctaz)**2/(ss - mz**2 + slfaa + ctaa)))
!      cpropz = cmplx(1._dp/(ss - mz**2 + slfzz + ctzz 
!     .     - (slfaz + ctaz)/(ss + slfaa + ctaa)))
!      cpropaz = -(slfaz + ctaz)/(ss + slfaa + ctaa)*cpropz
!      cpropza = -(slfaz + ctaz)/(ss - mz**2 + slfzz + ctzz)*cpropa

!      print*, 'propa, cpropa: ', propa, cpropa
!      print*, 'propz, cpropz: ', propz, cpropz
!      print*, 'propaz, cpropaz: ', propa*propz, -cpropaz/(slfaz+ctaz)
!      print*, (slfaz + ctaz)**2/(mz**2 + slfaa + ctaa)


!      stop

!      slfzz = 0._dp
!      ctzz = 0._dp



!      call sigw0(sigw,aemmz,sw2,ep)
!      call znenteghaon(ctww,aemmz,sw2,0.d0,ep,"WW")
!      sigw = sigw + ctww
!      dr = dreal(sigw + ctww)/mw**2 + aemmz/4.d0/pi/sw2*(6.d0 
!     .     + (7.d0 - 4.d0*sw2)/2.d0/sw2*dlog(cw2))
!      print*, 'sigw: ', sigw, ctww
!      call znen(sigw,aemmz,sw2,mw**2,ep,'WW')
!      print*, sigw
!      print*, dr
!      stop

!      write(6,*) 'self: ', self(-1,1),self(1,-1),self(-2,2),self(2,-2)
!      write(6,*) 'self+vert: ', 
!     .     self(2,-2)+vert(2,-2),self(-2,2)+vert(-2,2),
!     .     self(1,-1)+vert(1,-1),self(-1,1)+vert(-1,1)
!      write(6,*) 'vert: ', vert(-1,1),vert(1,-1),vert(-2,2),vert(2,-2)
!      write(6,*) 'box: ', bxtot(2,-2),bxtot(-2,2),bxtot(1,-1),
!     .     bxtot(-1,1)
!      write(6,*) 'tota: ', msq(2,-2), msq(-2,2), msq(1,-1), msq(-1,1)
!      write(6,*) "born1: ", born(-1,1),born(1,-1),born(-2,2),born(2,-2)
!      write(6,*) "born2: ", bbrn(-1,1),bbrn(1,-1),bbrn(-2,2),bbrn(2,-2)
!      stop
!      write(6,*) "xw: ", xw, 1._dp - mw**2/mz**2
!      pause
    


!      slfzz = (0._dp,0._dp)
!      ctzz = 0._dp
!      slfaa = (0._dp,0._dp)
!      ctaa = 0._dp
!      slfzz = (0._dp,0._dp)
