      real(dp) :: b1scale,q2scale,q1scale,b2scale
      common/bqscale/b1scale,q2scale,q1scale,b2scale
!$omp threadprivate(/bqscale/)
