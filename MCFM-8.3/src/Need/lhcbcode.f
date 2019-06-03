c--- 9/2016:  minimally-adapted by JC from the original code
c---          written by Philip Ilten <philten@cern.ch>
c
c     -------------------------------------------------------------------------
c     Perform the LHCb specific forward-backward cuts.
c
c     prts:  four-vectors of the particles in the event.
c     njets: number of jets in the event.
c     dir:   direction to perform the cuts, 1 is forward and 2 is backward.
c     -------------------------------------------------------------------------
      function lhcb_cuts(prts, njets, dir)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'plabel.f'
      include 'is_functions_com.f'
      include 'leptcuts.f'
      include 'nproc.f'
      include 'lhcb.f'
      include 'first.f'
      include 'mpicommon.f'
      logical lhcb_cuts
c     Arguments.
      real(dp):: prts(mxpart, 4)
      integer njets, dir
c     Functions.
      real(dp):: etarap, pt, R
      logical is_lepton, is_hadronic, is_heavy, is_photon, leptmask(mxpart)
c     Local variables.
      real(dp):: etan, ptn
      integer idx, idx1, idx2, idxj, idxb, idxg, nl, nj, nb
c     Saved variables.
      integer nprts
      logical heavy, onecl, twocl, gamma, hardest
      save heavy, onecl, twocl, gamma, nprts
!      data first/-1/
!      save first, heavy, onecl, twocl, gamma
!$omp threadprivate(heavy,onecl,twocl,gamma,nprts)

c     First pass for process.
      if (first) then
         first = .false.
         nprts = 2
         heavy = .false.
         twocl = .false.
         onecl = .false.
         gamma = .false.
         if (pt1_min > 0._dp) then
            onecl = .true.
         endif
         if (pt2_min > 0._dp) then
            twocl = .true.
         endif
         if (ptg_min > 0._dp) then
            gamma = .true.
         endif

c     Handle the masking and un-masking of charged leptons.
         leptmask(:)=.true.
         if (nproc == 71) then
            leptmask(4) = (cut_mode .ne. 2)
            leptmask(5) = (cut_mode .ne. 3)
            leptmask(6) = (cut_mode .ne. 3)
         else if (nproc == 76) then
            leptmask(3) = (cut_mode .ne. 2)
            leptmask(5) = (cut_mode .ne. 3)
            leptmask(6) = (cut_mode .ne. 3)
         else if (nproc == 181) then
            leptmask(4) = (cut_mode .ne. 2)
         else if (nproc == 186) then
            leptmask(3) = (cut_mode .ne. 2)
         endif

c     Print the particles.
!$omp master
         if (rank == 0) then
         write(6,*)
         write(6,*)  '************* LHCb code: particles *****************'
         endif
!$omp end master
         do idx = 3, mxpart
            if (is_heavy(idx)) then
               heavy = .true.
            endif
            if (plabel(idx) .ne. '') then
               nprts = nprts + 1
            else
               exit
            endif
!$omp master
            if (rank == 0) then
            write(6,*) '* particle ', idx, plabel(idx), 
     &           '                         *'
            endif
!$omp end master
             
         enddo
!$omp master
      if (rank == 0) then
         write(6,*)  '**************************************************
     &**'

c     Print the cuts.
         write(6,*)  
         write(6,*)  '******************  LHCb cuts  *********************'
         write(6,*)  '*        direction       =', dir_mode, 
     &        '            *'
         write(6,*)  '*        heavy hardest   =', hard_bjet,
     &        '                      *'
         write(6,99) '*        pt(lept)        >   ', pt1_min,
     &        ' GeV        *'
         write(6,99) '*        pt(lept)        <   ', pt1_max,
     &        ' GeV        *'
         write(6,99) '*       eta(lept)        >   ', eta1_min(dir),
     &        '            *'
         write(6,99) '*       eta(lept)        <   ', eta1_max(dir),
     &        '            *'
         if (twocl) then
            write(6,99) '*     pt(2nd lept)       >   ', pt2_min,
     &           ' GeV        *'
            write(6,99) '*     pt(2nd lept)       <   ', pt2_max,
     &           ' GeV        *'
            write(6,99) '*    eta(2nd lept)       >   ', eta2_min(dir),
     &           '            *'
            write(6,99) '*    eta(2nd lept)       <   ', eta2_max(dir),
     &           '            *'
         endif
         write(6,99) '*        pt(jet)         >   ', ptj_min,
     &        ' GeV        *'
         write(6,99) '*        pt(jet)         <   ', ptj_max,
     &        ' GeV        *'
         write(6,99) '*       eta(jet)         >   ', etaj_min(dir),
     &        '            *'
         write(6,99) '*       eta(jet)         <   ', etaj_max(dir),
     &        '            *'
         if (heavy) then
            write(6,99) '*        pt(heavy jet)   >   ', ptb_min,
     &           ' GeV        *'
            write(6,99) '*        pt(heavy jet)   <   ', ptb_max,
     &           ' GeV        *'
            write(6,99) '*       eta(heavy jet)   >   ', etab_min(dir),
     &           '            *'
            write(6,99) '*       eta(heavy jet)   <   ', etab_max(dir),
     &           '            *'
         endif
         if (gamma) then
            write(6,99) '*        pt(gamma)       >   ', ptg_min,
     &           ' GeV        *'
            write(6,99) '*        pt(gamma)       <   ', ptg_max,
     &           ' GeV        *'
            write(6,99) '*       eta(gamma)       >   ', etag_min(dir),
     &           '            *'
            write(6,99) '*       eta(gamma)       <   ', etag_max(dir),
     &           '            *'
         endif
         write(6,99) '*        pt(lept + jet)  >   ', pt1j_min,
     &        ' GeV        *'
         write(6,99) '*        pt(lept + jet)  <   ', pt1j_max,
     &        ' GeV        *'
         write(6,99) '*        R(lept, jet)    >   ', Rjlmin,
     &        '            *'
         write(6,*)  '****************************************************'
         call flush(6)
         endif
!$omp end master
 99      format(1x,a29,f10.3,a13)
      endif

c     Determine the highest pT lepton(s) and jet.
      lhcb_cuts = .true.
      hardest   = .false.
      idx1 = -1; idx2 = -1; idxj = -1; idxb = -1; idxg = -1
      nl = 0; nj = 0; nb = 0
      pt1(dir)  = 0; pt2(dir)  = 0; ptj(dir)  = 0; ptb(dir)  = 0;
      ptg(dir)  = 0
      eta1(dir) = 0; eta2(dir) = 0; etaj(dir) = 0; etab(dir) = 0;
      etag(dir) = 0
      do idx = 3, nprts
         ptn  = pt(idx, prts) 
         etan = etarap(idx, prts)
         if (is_lepton(idx) .and. leptmask(idx)) then
            if ((ptn > pt1(dir)) .and. 
     &           (etan .gt. eta1_min(dir)) .and. 
     &           (etan < eta1_max(dir))) then
               if ((pt1(dir) > pt2(dir)) .and.
     &              (eta1(dir) .gt. eta2_min(dir)) .and.
     &              (eta1(dir) < eta2_max(dir))) then
                  pt2(dir) = pt1(dir); idx2 = idx1
                  eta2(dir) = eta1(dir)
               endif
               pt1(dir) = ptn; idx1 = idx; eta1(dir) = etan
            else if ((ptn > pt2(dir)) .and.
     &              (etan .gt. eta2_min(dir)) .and.
     &              (etan < eta2_max(dir))) then
               pt2(dir) = ptn; idx2 = idx; eta2(dir) = etan
            endif
            if ((ptn > pt1_min) .and. 
     &           (etan .gt. eta1_min(dir)) .and. 
     &           (etan < eta1_max(dir))) then
               nl = nl + 1
            endif
         else if (is_hadronic(idx)) then
            if ((ptn > ptj(dir)) .and.
     &           (etan .gt. etaj_min(dir)) .and.
     &           (etan < etaj_max(dir))) then
               ptj(dir) = ptn; idxj = idx; etaj(dir) = etan
               hardest = is_heavy(idx)
            endif
            if (is_heavy(idx) .and. (ptn > ptb(dir)) .and.
     &           (etan .gt. etab_min(dir)) .and.
     &           (etan < etab_max(dir))) then
               ptb(dir) = ptn; idxb = idx; etab(dir) = etan
            endif
            if ((ptn > ptj_min) .and. 
     &           (etan .gt. etaj_min(dir)) .and. 
     &           (etan < etaj_max(dir))) then
               if (is_heavy(idx)) then
                  nb = nb + 1
               else
                  nj = nj + 1
               endif
            endif
         else if (is_photon(idx)) then
            if ((ptn > ptg(dir)) .and.
     &           (etan .gt. etag_min(dir)) .and.
     &           (etan < etag_max(dir))) then
               ptg(dir) = ptn; idxg = idx; etag(dir) = etan
            endif
         endif
      enddo

c     Cut on number of final state particles.
      if ((nl_min > 0) .and. (nl < nl_min)) then
!         write(6,*) 'failed nl_min'
         return
      else if ((nj_min > 0) .and. (nj < nj_min)) then
!         write(6,*) 'failed nj_min'
         return
      else if ((nb_min > 0) .and. (nb < nb_min)) then
!         write(6,*) 'failed nb_min'
         return
      endif

c     Cut on leptons.
      if (onecl) then
         if ((idx1 == -1) .or. (pt1(dir) < pt1_min) .or.
     &        (pt1(dir) .gt. pt1_max)) then
!            write(6,*) 'failed onecl lept',idx1,pt1(dir),pt1_min,pt1_max
            return
         else if (twocl) then
            if ((idx2 == -1) .or. (pt2(dir) < pt2_min) .or.
     &           (pt2(dir) .gt. pt2_max)) then
!               write(6,*) 'failed twocl lept',idx2,pt2(dir),pt2_min,pt2_max
               return
            endif
         endif
      endif

c     Cut on photons if the process contains them.
      if (gamma) then
         if ((idxg == -1) .or. (ptg(dir) < ptg_min) .or.
     &        (ptg(dir) .gt. ptg_max)) then
!            write(6,*) 'failed gamma'
            return
         endif
      endif

c     Cut on jets if the process contains them.
      if (njets .gt. 0) then

c     Cut on hardest jet pT.
         if ((idxj == -1) .or. (ptj(dir) < ptj_min) .or.
     &        (ptj(dir) .gt. ptj_max)) then
!            write(6,*) 'failed hardest jet'
            return
         endif
         
c     Cut on the hardest heavy jet.
         if (heavy) then

c     Cut on heavy jet must be hardest.
            if ((hard_bjet) .and. (.not. hardest)) then
!               write(6,*) 'failed heavy jet hardest'
               return
            endif

c     Cut on hardest heavy jet pT.
            if ((idxb == -1) .or. (ptb(dir) < ptb_min) .or.
     &           (ptb(dir) .gt. ptb_max)) then
!               write(6,*) 'failed heavy pt'
               return
            endif
         endif
         
c     Cut on pT(lepton + jet).
         pt1j(dir) = sqrt((prts(idx1, 1) + prts(idxj, 1))**2 + 
     &        (prts(idx1, 2) + prts(idxj, 2))**2)
         if ((pt1j(dir) < pt1j_min) .or. (pt1j(dir) .gt. pt1j_max))
     &        then
!            write(6,*) 'failed pt(lep+jet)'
            return
         endif
         
c     Cut on dR(lepton + jet)
         if (R(prts, idx1, idxj) < Rjlmin) then
!            write(6,*) 'failed dR(lep+jet)'
            return
         else if (twocl .and. (R(prts, idx2, idxj) < Rjlmin)) then
!            write(6,*) 'failed dR(lep+jet)'
            return
         endif
      endif

c     Calculate m(lepton + jet).
      m1j(dir) = sqrt((prts(idx1, 4) + prts(idxj, 4))**2 - 
     &     (prts(idx1, 1) + prts(idxj, 1))**2 - 
     &     (prts(idx1, 2) + prts(idxj, 2))**2 - 
     &     (prts(idx1, 3) + prts(idxj, 3))**2)

c     Set cuts as passed.
      if (dir == 2) then
         eta1 = -eta1; eta2 = -eta2; etaj = -etaj
      endif
      lhcb_cuts      = .false.
      lhcb_pass(dir) = .true.
      
      return
      end

c     -------------------------------------------------------------------------
c     Fill the LHCb histograms.
c 
c     prts: four-vectors of the particles in the event.
c     wt:   event weight.
c     wt2:  event weight squared.
c     nd:   the dipole contribution.
c     -------------------------------------------------------------------------
      subroutine lhcb_plots(prts, wt, wt2, nd)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nplot.f'
      include 'lhcb.f'
c     Arguments.
      real(dp):: prts(mxpart, 4), wt, wt2
      integer nd
c     Local variables.
      integer dir, dir1, dir2
c     Saved variables.
      integer tag,nplotmax
      integer, parameter:: tagbook=1, tagplot=2
      logical, save::first=.true.
      common/nplotmax/nplotmax
!      character*4 tag
!      logical first
!      data tag/'book'/first/.true./
!      save tag, first

c     Return if not an LHCb cut mode.
      if (cut_mode == 0) then
         return
      endif

c     Set the directions and weights for forward-backward events.
      if (dir_mode == 0) then
         wt = wt / 2._dp; wt2 = wt2 / 4._dp; dir1 = 1; dir2 = 2
      else
         dir1 = dir_mode; dir2 = dir_mode 
      endif

c     Book the histograms.
      if (first) then
         tag=tagbook
         call bookplot(nextnplot, tag, '#sigma',
     &        0._dp, 0._dp, 0._dp, -0.5_dp, 0.5_dp, 1._dp, 'lin')
         call bookplot(nextnplot + 1, tag, 'p_{T}(lep)',
     &        0._dp, 0._dp, 0._dp, 0._dp, 200._dp, 5._dp, 'lin')
         call bookplot(nextnplot + 2, tag, '#eta(lep)',
     &        0._dp, 0._dp, 0._dp, -6._dp, 6._dp, 0.5_dp, 'lin')
         call bookplot(nextnplot + 3, tag, 'p_{T}(jet)',
     &        0._dp, 0._dp, 0._dp, 0._dp, 200._dp, 5._dp, 'lin')
         call bookplot(nextnplot + 4, tag, '#eta(jet)',
     &        0._dp, 0._dp, 0._dp, -6._dp, 6._dp, 0.5_dp, 'lin')
         call bookplot(nextnplot + 5, tag, 'p_{T}(lep + jet)',
     &        0._dp, 0._dp, 0._dp, 0._dp, 200._dp, 5._dp, 'lin')
         call bookplot(nextnplot + 6, tag, 'p_{T}(lep)',
     &        0._dp, 0._dp, 0._dp, 25._dp, 225._dp, 10._dp, 'lin')
         call bookplot(nextnplot + 7, tag, 'p_{T}(jet)',
     &        0._dp, 0._dp, 0._dp, 50._dp, 200._dp, 10._dp, 'lin')
         call bookplot(nextnplot + 8, tag, 'p_{T}(lep + jet)',
     &        0._dp, 0._dp, 0._dp, 20._dp, 220._dp, 25._dp, 'lin')
         call bookplot(nextnplot + 9, tag, 'm(lep + jet)',
     &        0._dp, 0._dp, 0._dp, 0._dp, 200._dp, 25._dp, 'lin')
      else
         tag=tagplot
      endif
      
c     Fill the histograms.
      do dir = dir1, dir2
         if (lhcb_pass(dir)) then
            call bookplot(nextnplot , tag, '', 0._dp, wt, wt2,
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 1, tag, '', pt1(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 2, tag, '', eta1(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 3, tag, '', ptj(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 4, tag, '', etaj(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 5, tag, '', pt1j(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 6, tag, '', pt1(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 7, tag, '', ptj(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 8, tag, '', pt1j(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
            call bookplot(nextnplot + 9, tag, '', m1j(dir), wt, wt2, 
     &           0._dp, 0._dp, 0._dp, '')
         endif
      enddo
      nextnplot = nextnplot + 10

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(nextnplot)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=nextnplot
      endif
      
      end

c     -------------------------------------------------------------------------
c     Configure MCFM for LHCb usage.
c     -------------------------------------------------------------------------
      subroutine lhcb_config()
      implicit none
      include 'types.f'
      include 'leptcuts.f'
      include 'jetcuts.f'
      include 'lhcb.f'
      include 'energy.f'
      
c     Return if not LHCb mode.
      if (cut_mode == 0) then
         return
      endif

c     Set the forward pseudo-rapidity limits.
      eta1_max(1) = leptrapmax
      eta2_max(1) = leptrap2max
      etaj_min(1) = etajetmin
      etaj_max(1) = etajetmax
      etab_max(1) = etabjetmax
      etag_max(1) = gammrapmax
      if (eta1_max(1) < eta1_min(1)) then
         eta1_max(1) = 99
      endif
      if (eta2_max(1) < eta2_min(1)) then
         eta2_max(1) = 99
      endif
      if (etaj_max(1) < etaj_min(1)) then
         etaj_max(1) = 99
      endif
      if (etab_max(1) < etab_min(1)) then
         etab_max(1) = 99
      endif
      if (etag_max(1) < etag_min(1)) then
         etag_max(1) = 99
      endif

c     Set the backward pseudo-rapidity limits.
      eta1_min(2) = -eta1_max(1)
      eta1_max(2) = -eta1_min(1)
      eta2_min(2) = -eta2_max(1)
      eta2_max(2) = -eta2_min(1)
      etaj_min(2) = -etaj_max(1)
      etaj_max(2) = -etaj_min(1)
      etab_min(2) = -etab_max(1)
      etab_max(2) = -etab_min(1)
      etag_min(2) = -etag_max(1)
      etag_max(2) = -etag_min(1)

c     Set the pT limits.
      pt1_min = leptptmin
      pt2_min = leptpt2min
      ptj_min = ptjetmin
      ptb_min = ptbjetmin
      ptg_min = gammptmin
      if (pt1_max .le. pt1_min) then
         pt1_max = sqrts
      endif
      if (pt2_max .le. pt2_min) then
         pt2_max = sqrts
      endif
      if (ptj_max .le. ptj_min) then
         ptj_max = sqrts
      endif
      if (ptb_max .le. ptb_min) then
         ptb_max = sqrts
      endif
      if (ptg_max .le. ptg_min) then
         ptg_max = sqrts
      endif
      if (pt1j_max .le. pt1j_min) then
         pt1j_max = sqrts
      endif

c     Set the MCFM limits.
      ptjetmin   = 2
      etajetmin  = 0
      etajetmax  = 5
      ptbjetmin  = 2
      etabjetmax = 5
      end
