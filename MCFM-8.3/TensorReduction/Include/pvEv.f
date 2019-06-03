      complex(dp):: Ev(Nee*Nemax,-2:0)
      common/Ev/Ev
!$omp threadprivate(/Ev/)
