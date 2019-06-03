      real(dp):: ptjetmin,ptjetmax,etajetmin,etajetmax,
     & ptbjetmin,ptbjetmax,etabjetmin,etabjetmax
      common/jetcuts/ptjetmin,ptjetmax,etajetmin,etajetmax,
     & ptbjetmin,ptbjetmax,etabjetmin,etabjetmax
!$omp threadprivate(/jetcuts/)
