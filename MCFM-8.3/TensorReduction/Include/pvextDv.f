      complex(dp):: Dv(Ndd*Ndmax,-2:0)
      common/Dvext/Dv
!$omp threadprivate(/Dvext/)
