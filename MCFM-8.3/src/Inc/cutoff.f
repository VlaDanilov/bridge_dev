      real(dp) :: cutoff,cutoff_s
      common/cutoff/cutoff,cutoff_s
!$omp threadprivate(/cutoff/)
