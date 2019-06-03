      integer function pvextCcache(p1sq,p2sq,p3sq,m1s,m2s,m3s)
      implicit none
      include 'types.f'
      include 'pvCnames.f'
      include 'TRextclear.f'
      include 'TRconstants.f'
      include 'TRonshellcutoff.f'
      include 'pvRespectmaxcindex.f'
      include 'pvforcerecalc.f'
      include 'pvCitry.f'
      real(dp):: para(Pcc),p1sq,p2sq,p3sq,m1s,m2s,m3s
      integer:: j,jtable,Ntrue
      logical,save:: maxcindexrespected(Ncmax)
      real(dp),save:: tableC(Pcc,Ncmax)
      integer,save:: Nstore=0
!$omp threadprivate(tableC,Nstore,maxcindexrespected)
      
C--set the number of stored values to zero
      if (clear(3)) then
      clear(3)=.false.
      Nstore=0
      endif
      if (Nstore .ge. Ncmax) then
      print *
      print *, 'pvextCcache:Nstore .ge. Ncmax'
      print *, 'pvextCcache:Nstore,Ncmax',Nstore,Ncmax
      print *, 'Either adjust Ncmax in Cnames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      para(1)=p1sq
      para(2)=p2sq
      para(3)=p3sq
      para(4)=m1s
      para(5)=m2s
      para(6)=m3s
C if parameter set is found, set pvextCcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Pcc
        if (abs(para(j)-tableC(j,jtable)) .lt. 1d-8) Ntrue=Ntrue+1 
        enddo
        if (Ntrue .eq. Pcc) then
          pvextCcache=(jtable-1)*Ncc
          if    (((maxcindexrespected(jtable) .eqv. .true.)
     &     .and. (pvRespectmaxCindex .eqv. .false.))
     &           .or. (pvforcerecalc)) then
c--- although integral is cached, ranks higher than those
c--- set by pvmaxcindex are required so recalculate
            call pvextCfill(p1sq,p2sq,p3sq,m1s,m2s,m3s,pvextCcache)
            maxcindexrespected(jtable)=.false.  
          endif
          return
        endif
      enddo

C    if parameter set is not found, the values are 
C    not in the cache, so we have to calculate them
C    and put them in the cache starting at pvextCcache
 20   pvextCcache=Nstore*Ncc
      pvCitry(pvextCcache)=-1 ! label tensor as unchecked
      Nstore=Nstore+1
      do j=1,Pcc
        if(abs(para(j)) .lt. onshellcutoff) para(j)=zero
      enddo
      do j=1,Pcc
      tableC(j,Nstore)=para(j)
      enddo
      call pvextCfill(para(1),para(2),para(3),para(4),para(5),para(6),
     &             pvextCcache)
      maxcindexrespected(Nstore)=pvRespectmaxcindex  
      return
      end 
