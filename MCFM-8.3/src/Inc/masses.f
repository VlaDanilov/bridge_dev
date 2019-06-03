      real(dp):: 
     & md,mu,ms,mc,mb,mt,
     & mel,mmu,mtau,
     & hmass,hwidth,
     & wmass,wwidth,
     & zmass,zwidth,
     & twidth,
     & tauwidth,
     & mtausq,mcsq,mbsq
      common/masses/
     & md,mu,ms,mc,mb,mt,
     & mel,mmu,mtau,
     & hmass,hwidth,
     & wmass,wwidth,
     & zmass,zwidth,
     & twidth,
     & tauwidth,
     & mtausq,mcsq,mbsq
!$omp threadprivate(/masses/)
