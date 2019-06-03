      complex(dp):: Av(Naa*Namax,-2:0)
      common/Avext/Av
!$omp threadprivate(/Avext/)
