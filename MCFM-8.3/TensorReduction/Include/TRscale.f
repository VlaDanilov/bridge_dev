      real(dp)::scale,musq
      common/pvscale/scale,musq
!$omp threadprivate(/pvscale/)
