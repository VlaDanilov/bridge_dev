      integer function pvextAcache(m1sq)
      implicit none
      include 'types.f'
      include 'pvAnames.f'
      include 'TRconstants.f'
      include 'TRextclear.f'
      include 'TRonshellcutoff.f'
      real(dp):: para(Paa),m1sq
      real(dp),save:: tableA(Paa,Namax)  
      integer,save:: Nstore=0
      integer::j,jtable,Ntrue    
!$omp threadprivate(tableA,Nstore)

      if (clear(1)) then
      clear(1)=.false.
      Nstore=0
      endif


      if (Nstore .gt. NAmax) then
      print * 
      print *, 'pvextAcache: Nstore .gt. Namax'
      print *, 'pvextAcache:Nstore,Namax',Nstore,Namax
      print *, 'Either adjust Namax in Anames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      para(1)=m1sq

C if parameter set is found set pvextAcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Paa
        if (abs(para(j)-tableA(j,jtable)) .lt. 1d-8) Ntrue=Ntrue+1 
        enddo
      if (Ntrue .eq. Paa) then
      pvextAcache=(jtable-1)*Naa
      return
      endif
      enddo

C    if parameter set is not found we have to calculate
 20   pvextAcache=Nstore*Naa
      Nstore=Nstore+1
      do j=1,Paa
        if(abs(para(j)) .lt. onshellcutoff) para(j)=zero
      enddo
      do j=1,Paa
      tableA(j,Nstore)=para(j)
      enddo
      call pvextAfill(para(1),pvextAcache)
      return
      end 
