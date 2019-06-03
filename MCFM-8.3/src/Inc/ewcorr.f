      character*24 ewcorr
      integer kewcorr
      integer, parameter:: knone=1, ksudakov=2, kexact=3
      real(dp):: wt_noew
      common/kewcorr/kewcorr
      common/ewcorr/ewcorr
      common/wt_noew/wt_noew
!$omp threadprivate(/wt_noew/)
      
