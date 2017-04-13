      subroutine mcfmsub(r,er)
      implicit none
      include 'types.f'
c--- This is an entry point into MCFM (usually called by mcfm program)    
      
      real(dp):: r,er
      COMMON /check/tt2,tt3,tb,te,tt1,ti1,ti2,ttt1,ttt2,tm
      real(4):: tt2,tt3,tb,te,tt1,tsub1,tsub2,ti1,ti2,ttt1,ttt2,tm
      character*72 inputfile,workdir
      print *, "<<<<<   CALLING determinefilenames() and mcfmmain() HERE!!!!  >>>>>"           
      tt1 = secnds (0.0)
      call determinefilenames(inputfile,workdir)
      tt2 = secnds (tt1)
      tm = secnds (0.0)
      call mcfmmain(inputfile,workdir,r,er)
      tt3 = secnds (tm)
      return
      end
