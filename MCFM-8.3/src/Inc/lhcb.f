c     The pseudo-rapidity and pT requirements. The pseudo-rapidity
c     requirements are arrays of size 2, corresponding to the forward
c     backward requirement: 1 is the first lepton (hardest), 2 the
c     second lepton, j the jet, b the heavy jet, g the photon, and 1j
c     the combined first lepton and jet.
      real(dp):: 
     &     eta1_min(2), eta1_max(2), pt1_min, pt1_max,
     &     eta2_min(2), eta2_max(2), pt2_min, pt2_max,
     &     etaj_min(2), etaj_max(2), ptj_min, ptj_max,
     &     etab_min(2), etab_max(2), ptb_min, ptb_max,
     &     etag_min(2), etag_max(2), ptg_min, ptg_max,
     &     pt1j_min, pt1j_max
c     The variables stored for histogramming. All variables are arrays
c     of size 2 corresponding to forward and backward.
      real(dp)::
     &     eta1(2), eta2(2), etaj(2), etab(2), etag(2),
     &     pt1(2), pt2(2), ptj(2), ptb(2), ptg(2), pt1j(2), m1j(2)
c     The minimum number of letpons, jets, and heavy jets.
      integer
     &     nl_min, nj_min, nb_min
c     The cut mode and direction mode. The cut modes are as follows.
c     0: MCFM cuts.
c     1: LHCb-ANA-2014-076, cut on hardest lepton and jet in acceptance.
c        If the event contains heavy flavor and hard_bjet is true the hardest 
c        jet must be heavy.
c     2: Same as (1), but ignore the lepton from the W decay in nprocs 71, 
c        76, 181, and 186.
c     3: Same as (1), but ignore the leptons from the Z decay in nprocs 71
c        and 76.
c     The direction mode specifies how to handle asymmetric
c     pseudo-rapidity.
c     0: Use the event and flipped (in pseudo-rapidity) event.
c     1: Use only the event.
c     2: Use only the flipped event.
      integer
     &     cut_mode, dir_mode
c     Logical if the event passes the forward and backward requirements,
c     and logical if the hardest jet must be heavy.
      logical
     &     lhcb_pass(2), hard_bjet
      common/lhcb/
     &     eta1_min, eta1_max, pt1_min, pt1_max,
     &     eta2_min, eta2_max, pt2_min, pt2_max,
     &     etaj_min, etaj_max, ptj_min, ptj_max,
     &     etab_min, etab_max, ptb_min, ptb_max,
     &     etag_min, etag_max, ptg_min, ptg_max,
     &     pt1j_min, pt1j_max,
     &     eta1, eta2, etaj, etab, etag,
     &     pt1, pt2, ptj, ptb, ptg, pt1j, m1j,
     &     nl_min, nj_min, nb_min,
     &     cut_mode, dir_mode,
     &     lhcb_pass, hard_bjet
!$omp threadprivate(/lhcb/)
