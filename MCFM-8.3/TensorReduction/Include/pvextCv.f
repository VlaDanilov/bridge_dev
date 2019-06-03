      complex(dp):: Cv(Ncc*Ncmax,-2:0)
      common/Cvext/Cv
!$omp threadprivate(/Cvext/)
