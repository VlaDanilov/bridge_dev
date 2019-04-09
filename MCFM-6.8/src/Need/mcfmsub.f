      subroutine mcfmsub
c--- This is an entry point into MCFM (usually called by mcfm program)    
      implicit none
      character*72 inputfile,workdir
      call determinefilenames(inputfile,workdir)
      call mcfmmain(inputfile,workdir)
      end
