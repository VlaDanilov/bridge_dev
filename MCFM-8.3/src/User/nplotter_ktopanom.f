      subroutine nplotter_ktopanom(pjet,wt,wt2,nd)
          use ieee_arithmetic
          use types
          use jettagging
          use bbfrac_m, only: bbfrac
          use qsort_m
          use mathfun
          implicit none

          include 'first.f'
          include 'nproc.f'
          include 'jetlabel.f' ! jets
          include 'nplot.f' ! nextnplot
          include 'masses.f'
          include 'plabel.f'

          real(dp), intent(in) :: pjet(mxpart,4)
          real(dp), intent(in) :: wt,wt2
          integer, intent(in) :: nd

          real(dp) :: mainwgt, bbwgt

          integer, parameter:: tagbook=1, tagplot=2
          integer :: tag

          integer :: jetindices(3)

          integer :: bjetindices(3)
          integer :: nonbjetindices(3)

          integer :: bjetindices_bbfrac(3)
          integer :: nonbjetindices_bbfrac(3)

          real(dp) :: ptpure, etarappure, puremass, rpure
          real(dp) :: costheta, costheta_bbfrac
          integer :: i,n

          integer :: imax, imaxb

          logical :: failed, failed_bbfrac

          real(dp) :: ptop(4), ptopalt(4), ptoprec(4), ptoprest(4)
          real(dp) :: plep(4), plight(4), pnu(4)
          real(dp) :: pleprest(4), plightrest(4)
          real(dp) :: mtcand, mtcandalt, mtrec

          real(dp) :: pzn(2), pnu_cand1(4), pnu_cand2(4)

          real(dp) :: pw(4),pwrest(4),pwrest_bbfrac(4)
          real(dp) :: pwrest_w(4)
          real(dp) :: pspecrest(4), pspecrest_bbfrac(4)

          real(dp) :: plepWrest(4)
          real(dp) :: plepWrestproj(3)
          real(dp) :: coslstar, coslN, coslT
          real(dp) :: cosphiTstar, cosphiNstar

          real(dp) :: plepWrestproj_bbfrac(3)
          real(dp) :: coslstar_bbfrac, coslN_bbfrac, coslT_bbfrac
          real(dp) :: cosphiTstar_bbfrac, cosphiNstar_bbfrac

          real(dp) :: x(3),y(3),z(3)
          real(dp) :: st(3), st_bbfrac(3)
          real(dp) :: x_bbfrac(3), y_bbfrac(3), z_bbfrac(3)

          real(dp) :: xprod(3), yprod(3), zprod(3)
          real(dp) :: xprod_bbfrac(3), yprod_bbfrac(3), zprod_bbfrac(3)
          real(dp) :: beam(4), beamrest(4), beamrest_bbfrac(4)
          real(dp) :: direction 

          real(dp) :: cosxprod, cosyprod, coszprod
          real(dp) :: cosxprod_bbfrac, cosyprod_bbfrac, coszprod_bbfrac

          real(dp) :: ptop_bbfrac(4), ptopalt_bbfrac(4), ptoprec_bbfrac(4), ptoprest_bbfrac(4)
          real(dp) :: plight_bbfrac(4)
          real(dp) :: pleprest_bbfrac(4), plightrest_bbfrac(4)
          real(dp) :: mtcand_bbfrac, mtcandalt_bbfrac, mtrec_bbfrac

          integer :: nplotmax
          common/nplotmax/nplotmax

          if (first) then
              tag = tagbook
          else
              tag = tagplot
          endif

          jetindices(:) = [5,6,7]

c         jetpt = 0._dp
c         etarap = 0._dp
c         do i=1,jets
c           jetpt(i) = ptpure(pjet(4+i,:))
c           etarap(i) = etarappure(pjet(4+i,:))
c         enddo

c         ! it looks like it's already pt ordered, so the following is not necessary
c         call qsort_dp_with_jetindex(jetpt, jetindices)
c         write (*,*) jetpt, jetindices

          ! READ THIS IF YOU IMPLEMENT OR MODIFY THINGS:
          ! _bbfrac pieces contain two b-quarks in the final state
          ! and have to be handled separately due to the way things are set up
          ! they are associatiated with the bbwgt weight

          ! the other pieces contain only one b-quark in the final state
          ! and are associated with the mainwgt weight

          bjetindices = 0
          nonbjetindices = 0

          ! for main part with at most one b
          imax = 1
          imaxb = 1
          do i=1,jets
              if (iand(bitflag(5),jetcontent(jetindices(i))) > 0) then
                  bjetindices(imaxb) = jetindices(i)
                  imaxb = imaxb + 1
              else
                  nonbjetindices(imax) = jetindices(i)
                  imax = imax + 1
              endif
          enddo

          bjetindices_bbfrac = 0
          nonbjetindices_bbfrac = 0

          imax = 1
          imaxb = 1
          do i=1,jets
              if ((iand(bitflag(5),jetcontent(jetindices(i))) > 0) .or. (iand(bitflag(7),jetcontent(jetindices(i))) > 0)) then
                  bjetindices_bbfrac(imaxb) = jetindices(i)
                  imaxb = imaxb + 1
              else
                  nonbjetindices_bbfrac(imax) = jetindices(i)
                  imax = imax + 1
              endif
          enddo

          mainwgt = wt*(1-bbfrac(nd))
          bbwgt = wt*bbfrac(nd)


ccc       ! we already check for exactly two jets in gencuts.f for ktopanom

c         if ((.not. first) .and. jets /= 2) then
c             return
c         endif

ccc       !!! these conditions check for at least one b-jet and one light jet

          ! return if both main and bbfrac contributions are zero due to
          ! jet constraint
          if ( ((bjetindices(1) == 0 .or. nonbjetindices(1) == 0) .and.
     &         (bjetindices_bbfrac(1) == 0 .or. nonbjetindices_bbfrac(1) == 0)) .and.
     &          (.not. first) ) then
              return 
          endif

          ! otherwise only one contribution might vanish, so do not return

          ! set failed to .true. not when first == .true. !!
          failed = .false.
          failed_bbfrac = .false.

          if ((bjetindices(1) == 0 .or. nonbjetindices(1) == 0)
     &              .and. (.not. first)) then
              failed = .true.
              mainwgt = 0._dp
          endif

          if ((bjetindices_bbfrac(1) == 0 .or. nonbjetindices_bbfrac(1) == 0)
     &              .and. (.not. first)) then
              failed_bbfrac = .true.
              bbwgt = 0._dp
          endif

          !! WARNING: some observables below depend on the presence of
          !! a light quark jet. Additional checks need to be made
          !! when there is no such jet, for example in events with
          !! two jets, both b-tagged.

          if (.not. first) then

          ! W mass reconstruction
          if (plabel(3) == 'el') then
              plep(:) = pjet(3,:)
              pnu(:) = pjet(4,:)
          else
              plep(:) = pjet(4,:)
              pnu(:) = pjet(3,:)
          endif


CCC COMMENT THIS BLOCK OUT TO HAVE EXACT W RECONSTRUCTION
CCC OTHERWISE RECONSTRUCT NEUTRINO Z COMPONENT
CCC BY REQUIRING ON-SHELL W; SELECT SOLUTION WITH CLOSEST ON-SHELL TOP MASS
c         pzn(1) = (wmass**2*plep(3) + 2*plep(1)*plep(3)*pnu(1) + 2*plep(2)*plep(3)*pnu(2) - 
c    -    Sqrt((plep(1)**2 + plep(2)**2 + plep(3)**2)*
c    -      (wmass**4 - 4*(plep(2)*pnu(1) - plep(1)*pnu(2))**2 + 4*wmass**2*(plep(1)*pnu(1) + plep(2)*pnu(2)))))/
c    -  (2._dp*(plep(1)**2 + plep(2)**2))

c         pzn(2) = (wmass**2*plep(3) + 2*plep(1)*plep(3)*pnu(1) + 2*plep(2)*plep(3)*pnu(2) + 
c    -    Sqrt((plep(1)**2 + plep(2)**2 + plep(3)**2)*
c    -      (wmass**4 - 4*(plep(2)*pnu(1) - plep(1)*pnu(2))**2 + 4*wmass**2*(plep(1)*pnu(1) + plep(2)*pnu(2)))))/
c    -  (2._dp*(plep(1)**2 + plep(2)**2))

c         pnu_cand1(:) = pnu(:)
c         pnu_cand1(3) = pzn(1)
c         pnu_cand1(4) = sqrt(pnu(1)**2 + pnu(2)**2 + pzn(1)**2)

c         pnu_cand2(:) = pnu(:)
c         pnu_cand2(3) = pzn(2)
c         pnu_cand2(4) = sqrt(pnu(1)**2 + pnu(2)**2 + pzn(2)**2)

c         if (ieee_is_nan(pzn(1)) .or. ieee_is_nan(pzn(2))) then
c            return
c         else
c            if (  abs(puremass(pnu_cand1(:) + plep(:) + pjet(bjetindices(1),:))-mt) < 
c    &             abs(puremass(pnu_cand2(:) + plep(:) + pjet(bjetindices(1),:))-mt) ) then
c               pnu(:) = pnu_cand1(:)
c           else
c               pnu(:) = pnu_cand2(:)
c           endif
c         endif
CCC END W RECONSTRUCTION, NEUTRINO Z COMPONENT

          ! alternative smearing
          call wrecons(plep, pnu, pnu_cand1)
          pnu = pnu_cand1

          !! begin: reconstruction of top
          ptop(:) = 0._dp
          ptop_bbfrac(:) = 0._dp
          ptop(:) = pnu(:) + plep(:) + pjet(bjetindices(1),:)
          ptop_bbfrac(:) = pnu(:) + plep(:) + pjet(bjetindices_bbfrac(1),:)

          ptopalt(:) = 0._dp
          ptopalt_bbfrac(:) = 0._dp

          pw(:) = 0._dp
          pwrest(:) = 0._dp
          pwrest_bbfrac(:) = 0._dp

          ! when there is a second b-jet, we have another candidate for the top mass reconstruction
          ! although this will be overridden by default with the 2-jet requirement above
          if (bjetindices(2) > 0) then
              ptopalt(:) = pnu(:) + plep(:) + pjet(bjetindices(2),:)
              if ( abs(puremass(ptop)-mt) < abs(puremass(ptopalt)-mt) ) then
                mtrec = puremass(ptop)
                ptoprec = ptop
              else
                mtrec = puremass(ptopalt)
                ptoprec = ptopalt
              endif
          else
              mtrec = puremass(ptop)
              ptoprec = ptop
          endif

          if (bjetindices_bbfrac(2) > 0) then
              ptopalt_bbfrac(:) = pnu(:) + plep(:) + pjet(bjetindices_bbfrac(2),:)
              if ( abs(puremass(ptop_bbfrac)-mt) < abs(puremass(ptopalt_bbfrac)-mt) ) then
                mtrec_bbfrac = puremass(ptop_bbfrac)
                ptoprec_bbfrac = ptop_bbfrac
              else
                mtrec_bbfrac = puremass(ptopalt_bbfrac)
                ptoprec_bbfrac = ptopalt_bbfrac
              endif
          else
              mtrec_bbfrac = puremass(ptop_bbfrac)
              ptoprec_bbfrac = ptop_bbfrac
          endif

c         ! if top reconstruction fails, don't record
c         if (failed) then
c             mainwgt = 0._dp
c         endif

c         if (failed_bbfrac) then
c             bbwgt = 0._dp
c         endif

c         ! additional reconstructed top mass range

c         if (mtrec > 200._dp .or. mtrec < 130._dp) then
c           mainwgt = 0._dp
c         endif

c         if (mtrec_bbfrac > 200._dp .or. mtrec_bbfrac < 130._dp) then
c           bbwgt = 0._dp
c         endif

c         ! additional rapidity cuts on b jet 
c         if (bjetindices(1) > 0) then
c             if (abs(etarappure(pjet(bjetindices(1),:))) > 2.5_dp) then
c                 mainwgt = 0._dp
c             endif
c         endif

c         if (bjetindices_bbfrac(1) > 0) then
c             if (abs(etarappure(pjet(bjetindices_bbfrac(1),:))) > 2.5_dp) then
c                 bbwgt = 0._dp
c             endif
c         endif

c         ! additional rapidity cuts on light jet
c         if (nonbjetindices(1) > 0) then
c             if (abs(etarappure(pjet(nonbjetindices(1),:))) < 2.0_dp) then
c                 mainwgt = 0._dp
c             endif
c         endif

c         if (nonbjetindices_bbfrac(1) > 0) then
c             if (abs(etarappure(pjet(nonbjetindices_bbfrac(1),:))) < 2.0_dp) then
c                 bbwgt = 0._dp
c             endif
c         endif

          endif


          ! leading b jet
          n = nextnplot
          if (bjetindices(1) > 0) then
            call bookplot(n,tag,'ptb1',ptpure(pjet(bjetindices(1),:)),mainwgt,mainwgt**2,0d0,300d0,5d0,'lin')
          endif
          if (bjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'ptb1',ptpure(pjet(bjetindices_bbfrac(1),:)),bbwgt,bbwgt**2,0d0,300d0,5d0,'lin')
          endif

          n = n+1
          if (bjetindices(1) > 0) then
            call bookplot(n,tag,'etab1',etarappure(pjet(bjetindices(1),:)),mainwgt,mainwgt**2,-4.0d0,4.0d0,0.1d0,'lin')
          endif
          if (bjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'etab1',etarappure(pjet(bjetindices_bbfrac(1),:)),bbwgt,bbwgt**2,-4.0d0,4.0d0,0.1d0,'lin')
          endif

c         ! this is some debugging histogram to check for rapidity spectrum asymmetries
c         ! things have been debugged, so this serves as a check if uncertainties are correctly estimated:
c         ! every bin should be consistent with zero
          n = n+1
          if (bjetindices(1) > 0) then
            call bookplot(n,tag,'etab1asy',abs(etarappure(pjet(bjetindices(1),:))),
     &              sign(1d0,etarappure(pjet(bjetindices(1),:)))*mainwgt,mainwgt**2,0d0,4.0d0,0.25d0,'lin')
          endif
          if (bjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'etab1asy',abs(etarappure(pjet(bjetindices_bbfrac(1),:))),
     &              sign(1d0,etarappure(pjet(bjetindices_bbfrac(1),:)))*bbwgt,bbwgt**2,0d0,4.0d0,0.25d0,'lin')
          endif
c         ! end debugging histogram

          ! subleading b jet
          n = n+1
          if (bjetindices_bbfrac(2) > 0) then
            call bookplot(n,tag,'ptb2',ptpure(pjet(bjetindices_bbfrac(2),:)),bbwgt,bbwgt**2,0d0,300d0,5d0,'lin')
          endif

          n = n+1
          if (bjetindices(1) > 0) then
            call bookplot(n,tag,'mbl',puremass(pjet(bjetindices(1),:) + plep(:)), mainwgt, mainwgt**2, 0d0, 200d0, 5d0, 'lin')
          endif
          if (bjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'mbl',puremass(pjet(bjetindices_bbfrac(1),:) + plep(:)), bbwgt, bbwgt**2, 0d0, 200d0, 5d0, 'lin')
          endif

          n = n+1
          if (bjetindices(1) > 0) then
            call bookplot(n,tag,'ptbl',ptpure(pjet(bjetindices(1),:) + plep(:)), mainwgt, mainwgt**2, 0d0, 300d0, 5d0, 'lin')
          endif
          if (bjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'ptbl',ptpure(pjet(bjetindices_bbfrac(1),:) + plep(:)), bbwgt, bbwgt**2, 0d0, 300d0, 5d0, 'lin')
          endif


          ! angular correlations, observables in top rest frame

          if (.not. first) then

          ptoprest(:) = 0._dp
          ptoprest(4) = mtrec

          ptoprest_bbfrac(:) = 0._dp
          ptoprest_bbfrac(4) = mtrec_bbfrac

          plight(:) = pjet(nonbjetindices(1),:)
          plight_bbfrac(:) = pjet(nonbjetindices_bbfrac(1),:)

          call boostx(plep,ptoprec,ptoprest,pleprest)
          call boostx(plight,ptoprec,ptoprest,plightrest)

          costheta = sum(pleprest(1:3)*plightrest(1:3))
     &          /sqrt(sum(pleprest(1:3)**2))/sqrt(sum(plightrest(1:3)**2))


          call boostx(plep,ptoprec_bbfrac,ptoprest_bbfrac,pleprest_bbfrac)
          call boostx(plight_bbfrac,ptoprec_bbfrac,ptoprest_bbfrac,plightrest_bbfrac)

          costheta_bbfrac = sum(pleprest_bbfrac(1:3)*plightrest_bbfrac(1:3))
     &          /sqrt(sum(pleprest_bbfrac(1:3)**2))/sqrt(sum(plightrest_bbfrac(1:3)**2))

          ! construction of coordinate system in top rest frame
          ! the following observables probe the decay part (1005.5382)

          pw = pnu(:) + plep(:)
          pwrest_W(:) = 0._dp
          pwrest_W(4) = puremass(pw)

          ! pwrest is W vector in top rest frame
          ! pwrest_W is W vector in W rest fram

          ! pwrest defines z axis
          call boostx(pw,ptoprec,ptoprest,pwrest)
          call boostx(pw,ptoprec_bbfrac,ptoprest_bbfrac,pwrest_bbfrac)

          ! normalize axis
          z(1:3) = pwrest(1:3)
          z(:) = z(:) / euclideanNorm(z)

          z_bbfrac(1:3) = pwrest_bbfrac(1:3)
          z_bbfrac(:) = z_bbfrac(:) / euclideanNorm(z_bbfrac)

          ! pspecrest defines \hat{s}_t
          call boostx(pjet(nonbjetindices(1),:), ptoprec, ptoprest, pspecrest)
          call boostx(pjet(nonbjetindices_bbfrac(1),:), ptoprec_bbfrac, ptoprest_bbfrac, pspecrest_bbfrac)

          st(1:3) = pspecrest(1:3)
          st(:) = st(:) / euclideanNorm(st)

          st_bbfrac(1:3) = pspecrest_bbfrac(1:3)
          st_bbfrac(:) = st_bbfrac(:) / euclideanNorm(st_bbfrac)

          y(:) = crossprod(z, st)
          y(:) = y(:) / euclideanNorm(y)
          y_bbfrac(:) = crossprod(z_bbfrac, st_bbfrac)
          y_bbfrac(:) = y_bbfrac(:) / euclideanNorm(y_bbfrac)

          x(:) = crossprod(y,z)
          x(:) = x(:) / euclideanNorm(x)
          x_bbfrac(:) = crossprod(y_bbfrac, z_bbfrac)
          x_bbfrac(:) = x_bbfrac(:) / euclideanNorm(x_bbfrac)

          ! lepton in W rest frame
          call boostx(plep, pw, pwrest_W, plepWrest)

          coslstar = sum(z(:)*plepWrest(1:3))/ euclideanNorm(plepWrest(1:3))
          coslstar_bbfrac = sum(z_bbfrac(1:3)*plepWrest(1:3))/ euclideanNorm(plepWrest(1:3))

          coslN = sum(-y(:)*plepWrest(1:3)) / euclideanNorm(plepWrest(1:3))
          coslN_bbfrac = sum(-y_bbfrac(:)*plepWrest(1:3)) / euclideanNorm(plepWrest(1:3))

          coslT = sum(x(:)*plepWrest(1:3)) / euclideanNorm(plepWrest(1:3))
          coslT_bbfrac = sum(x_bbfrac(:)*plepWrest(1:3)) / euclideanNorm(plepwRest(1:3))

          plepWrestproj(:) = plepWrest(1:3) - sum(plepWrest(1:3)*z(1:3)) * z(1:3)
          plepWrestproj_bbfrac(:) = plepWrest(1:3) - sum(plepWrest(1:3)*z(1:3)) * z_bbfrac(1:3)

          ! angle between N (-y) and the projection
          cosphiNstar = sum(-y(:)*plepWrestproj(:)) / euclideanNorm(plepWrestproj)
          cosphiNstar_bbfrac = sum(-y_bbfrac(:)*plepWrestproj_bbfrac(:)) / euclideanNorm(plepWrestproj_bbfrac)

          ! angle between T (x) and the projection
          cosphiTstar = sum(x(:)*plepWrestproj(:)) / euclideanNorm(plepWrestproj)
          cosphiTstar_bbfrac = sum(x_bbfrac(:)*plepWrestproj_bbfrac(:)) / euclideanNorm(plepWrestproj_bbfrac)


          ! construction of the coordinate system that probes the production part
          ! 1404.1585

          zprod(:) = st(:)
          zprod_bbfrac(:) = st_bbfrac(:)

          ! opposite sign for p1, since MCFM uses negative incoming momenta
          if ( (-pjet(1,3)) * pjet(nonbjetindices(1),3) > 0d0 ) then
              beam(:) = pjet(1,:)
          else
              beam(:) = pjet(2,:)
          endif

          call boostx(beam, ptoprec, ptoprest, beamrest)
          call boostx(beam, ptoprec_bbfrac, ptoprest_bbfrac, beamrest_bbfrac)

          yprod(:) = crossprod(zprod, beamrest(1:3))
          yprod(:) = yprod(:) / euclideanNorm(yprod)

          yprod_bbfrac(:) = crossprod(zprod_bbfrac, beamrest_bbfrac(1:3))
          yprod_bbfrac(:) = yprod_bbfrac(:) / euclideanNorm(yprod_bbfrac)

          xprod(:) = crossprod(yprod, zprod)
          xprod_bbfrac(:) = crossprod(yprod_bbfrac, zprod_bbfrac)

          ! angles of lepton w.r.t xprod,yprod,zprod

          cosxprod = sum(xprod(:)*pleprest(1:3)) / euclideanNorm(pleprest(1:3))
          cosyprod = sum(yprod(:)*pleprest(1:3)) / euclideanNorm(pleprest(1:3))
          coszprod = sum(zprod(:)*pleprest(1:3)) / euclideanNorm(pleprest(1:3))

          cosxprod_bbfrac = sum(xprod_bbfrac(:)*pleprest_bbfrac(1:3)) / euclideanNorm(pleprest_bbfrac(1:3))
          cosyprod_bbfrac = sum(yprod_bbfrac(:)*pleprest_bbfrac(1:3)) / euclideanNorm(pleprest_bbfrac(1:3))
          coszprod_bbfrac = sum(zprod_bbfrac(:)*pleprest_bbfrac(1:3)) / euclideanNorm(pleprest_bbfrac(1:3))

          endif


          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'azfb', 0.5d0, sign(1d0,coszprod)*mainwgt, mainwgt**2, 0d0,1d0,1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'azfb', 0.5d0, sign(1d0,coszprod_bbfrac)*bbwgt, bbwgt**2, 0d0,1d0,1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'axfb', 0.5d0, sign(1d0,cosxprod)*mainwgt, mainwgt**2, 0d0,1d0,1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'axfb', 0.5d0, sign(1d0,cosxprod_bbfrac)*bbwgt, bbwgt**2, 0d0,1d0,1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'ayfb', 0.5d0, sign(1d0,cosyprod)*mainwgt, mainwgt**2, 0d0,1d0,1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'ayfb', 0.5d0, sign(1d0,cosyprod_bbfrac)*bbwgt, bbwgt**2, 0d0,1d0,1d0,'lin')
          endif


          n=n+1
          if (.not.failed) then
              if (ieee_is_nan(cosxprod)) then
                  call abort
              endif
              call bookplot(n,tag,'cosxprod', cosxprod, mainwgt, mainwgt**2, -1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              if (ieee_is_nan(cosxprod_bbfrac)) then
                  call abort
              endif
              call bookplot(n,tag,'cosxprod', cosxprod_bbfrac, bbwgt, bbwgt**2, -1d0,1d0,0.1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'cosyprod', cosyprod, mainwgt, mainwgt**2, -1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'cosyprod', cosyprod_bbfrac, bbwgt, bbwgt**2, -1d0,1d0,0.1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'coszprod', coszprod, mainwgt, mainwgt**2, -1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'coszprod', coszprod_bbfrac, bbwgt, bbwgt**2, -1d0,1d0,0.1d0,'lin')
          endif

          ! depends on top reconstruction
          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'coslstar', coslstar, mainwgt, mainwgt**2,-1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'coslstar', coslstar_bbfrac, bbwgt, bbwgt**2,-1d0,1d0,0.1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'coslN', coslN, mainwgt, mainwgt**2,-1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'coslN', coslN_bbfrac, bbwgt, bbwgt**2,-1d0,1d0,0.1d0,'lin')
          endif


          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'coslT', coslT, mainwgt, mainwgt**2,-1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'coslT', coslT_bbfrac, bbwgt, bbwgt**2,-1d0,1d0,0.1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'acoslstarfb', 0.5d0, sign(1d0,coslstar)*mainwgt, mainwgt**2, 0d0,1d0,1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'acoslstarfb', 0.5d0, sign(1d0,coslstar_bbfrac)*bbwgt, bbwgt**2, 0d0,1d0,1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'acoslnfb', 0.5d0, sign(1d0,cosln)*mainwgt, mainwgt**2, 0d0,1d0,1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'acoslnfb', 0.5d0, sign(1d0,cosln_bbfrac)*bbwgt, bbwgt**2, 0d0,1d0,1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'acosltfb', 0.5d0, sign(1d0,coslt)*mainwgt, mainwgt**2, 0d0,1d0,1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'acosltfb', 0.5d0, sign(1d0,coslt_bbfrac)*bbwgt, bbwgt**2, 0d0,1d0,1d0,'lin')
          endif


          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'cosphiNstar', cosphiNstar, mainwgt, mainwgt**2,-1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'cosphiNstar', cosphiNstar_bbfrac, bbwgt, bbwgt**2,-1d0,1d0,0.1d0,'lin')
          endif

          n=n+1
          if (.not.failed) then
              call bookplot(n,tag,'cosphiTstar', cosphiTstar, mainwgt, mainwgt**2,-1d0,1d0,0.1d0,'lin')
          endif
          if (.not.failed_bbfrac) then
              call bookplot(n,tag,'cosphiTstar', cosphiTstar_bbfrac, bbwgt, bbwgt**2,-1d0,1d0,0.1d0,'lin')
          endif


          n=n+1
          if (.not. failed) then
              call bookplot(n,tag,'costheta',costheta,mainwgt,mainwgt**2,-1d0,1d0,0.1d0,'lin')
          endif
          if (.not. failed_bbfrac) then
              call bookplot(n,tag,'costheta',costheta_bbfrac,bbwgt,bbwgt**2,-1d0,1d0,0.1d0,'lin')
          endif

          n=n+1
          if (.not. failed) then
            call bookplot(n,tag,'mt',mtrec,mainwgt,mainwgt**2,120.5d0,220.5d0,1.0d0,'lin')
          endif
          if (.not. failed_bbfrac) then
            call bookplot(n,tag,'mt',mtrec_bbfrac,bbwgt,bbwgt**2,120.5d0,220.5d0,1.0d0,'lin')
          endif

          n=n+1
          if (.not. failed) then
            call bookplot(n,tag,'mtlarge',mtrec,mainwgt,mainwgt**2,0d0,1000d0,20d0,'lin')
          endif
          if (.not. failed_bbfrac) then
            call bookplot(n,tag,'mtlarge',mtrec_bbfrac,bbwgt,bbwgt**2,0d0,1000d0,20d0,'lin')
          endif

          n=n+1
          if (.not. failed) then
            if (nonbjetindices(1) > 0) then
                call bookplot(n,tag,'m_tj',puremass(ptop(:)+pjet(nonbjetindices(1),:)),
     &                                  mainwgt,mainwgt**2, 180.5d0,280.5d0,1.0d0,'lin')
            endif
          endif
          if (.not. failed_bbfrac) then
            if (nonbjetindices_bbfrac(1) > 0) then
                call bookplot(n,tag,'m_tj',puremass(ptop_bbfrac(:)+pjet(nonbjetindices_bbfrac(1),:)),
     &                                  mainwgt,mainwgt**2, 180.5d0,280.5d0,1.0d0,'lin')
            endif
          endif

          n=n+1
          call bookplot(n,tag,'mtransW', sqrt(pw(4)**2 - pw(3)**2), mainwgt, mainwgt**2, 0d0,300d0,5d0,'lin')
          call bookplot(n,tag,'mtransW', sqrt(pw(4)**2 - pw(3)**2), bbwgt, bbwgt**2, 0d0,300d0,5d0,'lin')

          !! end: top reconstruction and angular correlation

          n=n+1
          if (bjetindices(1) > 0 .and. nonbjetindices(1) > 0) then
              call bookplot(n,tag,'delR_jb',rpure(pjet(bjetindices(1),:), pjet(nonbjetindices(1),:)),
     &                  mainwgt, mainwgt**2, 0d0, 5d0, 0.1d0, 'lin')
          endif
          if (bjetindices_bbfrac(1) > 0 .and. nonbjetindices_bbfrac(1) > 0) then
              call bookplot(n,tag,'delR_jb',rpure(pjet(bjetindices_bbfrac(1),:), pjet(nonbjetindices_bbfrac(1),:)),
     &                  bbwgt, bbwgt**2, 0d0, 5d0, 0.1d0, 'lin')
          endif

          n=n+1
          call bookplot(n,tag,'ptW', ptpure(pnu(:) + plep(:)), mainwgt, mainwgt**2,0d0,300d0,5d0,'lin')
          call bookplot(n,tag,'ptW', ptpure(pnu(:) + plep(:)), bbwgt, bbwgt**2,0d0,300d0,5d0,'lin')

          n=n+1
          if (.not. failed) then
              call bookplot(n,tag,'pttoprec', ptpure(ptoprec), mainwgt, mainwgt**2, 0d0, 200d0, 5d0, 'lin')
          endif
          if (.not. failed_bbfrac) then
              call bookplot(n,tag,'pttoprec', ptpure(ptoprec_bbfrac), bbwgt, bbwgt**2, 0d0, 200d0, 5d0, 'lin')
          endif

          n=n+1
          if (.not. failed) then
              call bookplot(n,tag,'etatoprec', etarappure(ptoprec), mainwgt, mainwgt**2, -5d0, 5d0, 0.1d0, 'lin')
          endif
          if (.not. failed_bbfrac) then
              call bookplot(n,tag,'etatoprec', etarappure(ptoprec_bbfrac), bbwgt, bbwgt**2, -5d0, 5d0, 0.1d0, 'lin')
          endif

          ! leading non b jet
          n = n+1
          if (nonbjetindices(1) > 0) then
            call bookplot(n,tag,'ptj1',ptpure(pjet(nonbjetindices(1),:)),mainwgt,mainwgt**2,0d0,300d0,5d0,'lin')
          endif
          if (nonbjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'ptj1',ptpure(pjet(nonbjetindices_bbfrac(1),:)),bbwgt,bbwgt**2,0d0,300d0,5d0,'lin')
          endif

          n = n+1
          if (nonbjetindices(1) > 0) then
            call bookplot(n,tag,'etaj1',etarappure(pjet(nonbjetindices(1),:)),mainwgt,mainwgt**2,-4.0d0,4.0d0,0.1d0,'lin')
          endif
          if (nonbjetindices_bbfrac(1) > 0) then
            call bookplot(n,tag,'etaj1',etarappure(pjet(nonbjetindices_bbfrac(1),:)),bbwgt,bbwgt**2,-4.0d0,4.0d0,0.1d0,'lin')
          endif

          ! subleading non b jet
          n = n+1
          if (nonbjetindices(2) > 0) then
            call bookplot(n,tag,'ptj2',ptpure(pjet(nonbjetindices(2),:)),mainwgt,mainwgt**2,0d0,300d0,5d0,'lin')
          endif
          if (nonbjetindices_bbfrac(2) > 0) then
            call bookplot(n,tag,'ptj2',ptpure(pjet(nonbjetindices_bbfrac(2),:)),bbwgt,bbwgt**2,0d0,300d0,5d0,'lin')
          endif

          ! electron pt and rapidity
          n = n + 1
          call bookplot(n,tag,'pte',ptpure(plep(:)),mainwgt,mainwgt**2,0d0,300d0,5d0,'lin')
          call bookplot(n,tag,'pte',ptpure(plep(:)),bbwgt,bbwgt**2,0d0,300d0,5d0,'lin')

          n = n + 1
          call bookplot(n,tag,'etae',etarappure(plep(:)),mainwgt,mainwgt**2,-4.0d0,4.0d0,0.1d0,'lin')
          call bookplot(n,tag,'etae',etarappure(plep(:)),bbwgt,bbwgt**2,-4.0d0,4.0d0,0.1d0,'lin')

          ! jet binning
          n = n + 1
          call bookplot(n,tag,'jets',real(jets,dp),mainwgt,mainwgt**2,-0.5d0,3.5d0,1d0,'lin')
          call bookplot(n,tag,'jets',real(jets,dp),bbwgt,bbwgt**2,-0.5d0,3.5d0,1d0,'lin')

          n = n + 1
          call bookplot(n,tag,'totalcross', 0.5d0,mainwgt,mainwgt**2,0d0,1d0,1d0,'lin')
          call bookplot(n,tag,'totalcross', 0.5d0,bbwgt,bbwgt**2,0d0,1d0,1d0,'lin')


          if (first) then
              first = .false.
              nplotmax = n
          endif

      end

