      function gencuts(pjet,njets)
       use types
       implicit none
      logical:: gencuts
************************************************************************
*   Author: J.M. Campbell, 26th November 2013                          *
*                                                                      *
*   This routine calls specific cutting routines if requested,         *
*   otherwise just uses the cuts specified in the input file           *
*                                                                      *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'runstring.f'
      include 'first.f'
      include 'lhcb.f'
      include 'specialcuts.f'
      logical:: failed,gencuts_Zt,gencuts_WZjj,gencuts_input
      logical, save :: makeVBScuts,makeATLAS_sscuts,
     &     makeCMS_hzz,makeCMS_hzz_vbf,makeVHbb,makeVHWW,
     &     makeATLAS_gaga2
      integer:: njets
      logical gencuts_VHbb,gencuts_VHWW,gencuts_ATLAS_gaga2,lhcb_cuts
      real(dp):: pjet(mxpart,4)
!$omp threadprivate(makeVBScuts,makeATLAS_sscuts)
!$omp threadprivate(makeCMS_hzz,makeCMS_hzz_vbf)
      
      if (first) then
        first=.false.
        makeVBScuts=(index(runstring,'VBS') > 0)
        makeATLAS_sscuts=(index(runstring,'ATLAS_ss') > 0)
        makeCMS_hzz=(index(runstring,'CMS_hzz') > 0)
        makeCMS_hzz_vbf=(index(runstring,'CMS_hzz_vbf') > 0)
        makeVHbb=(index(runstring,'VHbb') > 0) 
        makeVHWW=(index(runstring,'VHWW') > 0) 
        makeATLAS_gaga2=(index(runstring,'ATLASgaga') > 0) 
      endif
      
      gencuts=.false.
      if (makeVHbb) then
        gencuts=gencuts_VHbb(pjet)
        return
      endif
      if (makeATLAS_gaga2) then
        gencuts=gencuts_ATLAS_gaga2(pjet)
        return
      endif
      if (makeVHWW) then
        gencuts=gencuts_VHWW(pjet)
        return
      endif
      
      if (makeVBScuts) then
        call VBS(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif

      if (makeATLAS_sscuts) then
        call ATLAS_ss(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif

      if (makeCMS_hzz_vbf) then
        call CMS_hzz_vbf(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif
      
      if (makeCMS_hzz) then
        call CMS_hzz(pjet,failed)
        if (failed) gencuts=.true.
        return
      endif
      
      if (cut_mode > 0) then
        lhcb_pass(1) = .false.
        lhcb_pass(2) = .false.
        if (dir_mode == 0) then
           gencuts = lhcb_cuts(pjet, njets, 1)
           gencuts = lhcb_cuts(pjet, njets, 2) .and. gencuts
        else
           gencuts = lhcb_cuts(pjet, njets, dir_mode)
        endif
        return
      endif

c     if (kcase==ktopanom) then
c         ! if there are more than 2 jets, return early
c         ! rest of checks (one b, one light) needs to be done in plotting routine
c         if (njets /= 2) then
c             gencuts = .true.
c             return
c         endif
c     endif

      
c--- Default: use the cuts from the input file
      gencuts=gencuts_input(pjet,njets)
        
      return


      end
 



 
 
c      if ( (kcase==kHZZ_4l)
c     & .or.(kcase==kHZZ_tb)
c     & .or.(kcase==kHZZint)
c     & .or.(kcase==kHZZHpi)
c     & .or.(kcase==kggZZ4l) 
c     & .or.(kcase==kHZZqgI)) then 
c        call CMS_hzz(pjet,failed)
c        if (failed) gencuts=.true.
c        return
c      endif

c      if ( (kcase==kHWW_4l)
c     & .or.(kcase==kHWW_tb)
c     & .or.(kcase==kHWWint)
c     & .or.(kcase==kHWWHpi)
c     & .or.(kcase==kggWW4l) 
c     & .or.(kcase==kWWqqbr)) then 
c        call ATLAS_hww2013(pjet,failed)
c        if (failed) gencuts=.true.
c        return
c      endif

c -- CMS FCNC cuts for Z_tdkj
c      if ( kcase==kZ_tdkj .and. runstring(1:3) == 'CMS' ) then
c         gencuts=gencuts_Zt(pjet,njets)
c         return
c      endif

c      if ( kcase==kWpmZjj .and. runstring(1:3) == 'CMS' ) then
c         gencuts=gencuts_WZjj(pjet,njets)
c         return
c      endif

 
