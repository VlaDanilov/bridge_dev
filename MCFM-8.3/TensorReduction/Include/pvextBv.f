      complex(dp):: Bv(Nbb*Nbmax,-2:0)
      common/Bvext/Bv
!$omp threadprivate(/Bvext/)
