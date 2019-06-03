      real(dp):: msq_mix(0:3,fn:nf,fn:nf)
      common/msq_mix/msq_mix
!$omp threadprivate(/msq_mix/)
