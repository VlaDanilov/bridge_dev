      subroutine read_jetcuts(read_ptmin,read_etamin,read_etamax)
      implicit none
      include 'clustering.f'
      include 'jetcuts.f'
      integer nargs
      character*72 jetcutsfile
      double precision read_ptmin,read_etamin,read_etamax
      double precision ptmin,etamax,ptmin_tev,etamax_tev,
     . ptmin_lhc,etamax_lhc
      logical useTevcuts,useLHCcuts
      
c      if (newinput) then
        read_ptmin=ptjetmin
        read_etamin=etajetmin
        read_etamax=etajetmax
        return
c      endif
      
      nargs=iargc()
      if (nargs .eq. 2) then
      call getarg(2,jetcutsfile)
      else
      jetcutsfile='jetcuts.DAT'
      endif      
                                     
      open(unit=21,file=jetcutsfile,status='old',err=99)
      call checkversion(21,jetcutsfile)
      read(21,*) algorithm
      read(21,*) inclusive
      read(21,*) useTevcuts
      read(21,*) useLHCcuts
      read(21,*) ptmin
      read(21,*) etamax
      read(21,*) ptmin_tev
      read(21,*) etamax_tev
      read(21,*) ptmin_lhc
      read(21,*) etamax_lhc
      close(21)
      
      if     (useTevcuts) then
c--- preset cuts for the Tevatron
        read_ptmin=ptmin_tev         
        read_etamax=etamax_tev         
      elseif (useLHCcuts) then
c--- preset cuts for the LHC
        read_ptmin=ptmin_lhc         
        read_etamax=etamax_lhc         
      else
c--- generic cuts
        read_ptmin=ptmin         
        read_etamax=etamax
      endif         
      
      if (useTevcuts .and. useLHCcuts) then
       write(6,*) 'Cannot use both Tevatron and LHC cuts in jetcuts.DAT'
       stop
      endif
         
      if ((algorithm .ne. 'ktal') .and. (algorithm .ne. 'cone')) then
       write(6,*)
       write(6,*) 'Incorrect choice of algorithm in jetcuts.DAT, use:'
       write(6,*) '    '''//'ktal'//''' for kt algorithm'
       write(6,*) '    '''//'cone'//''' for cone algorithm'
       stop
      endif
         
      return
      
   99 write(6,*) 'Error reading ',jetcutsfile
      call flush(6)
      stop
      end
      
      
      
