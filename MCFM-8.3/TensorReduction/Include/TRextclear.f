      logical clear(1:5)
      common/clearext/clear
!$omp threadprivate(/clearext/)
