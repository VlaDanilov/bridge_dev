      subroutine mtxt(n,m)
          use types
          implicit none
          character (len=255) :: runname
          common/runname/runname
          integer, intent(in) :: n,m
          include 'histo.f'

          integer :: j

          if (book(n) .NE. 'YES') return

          open(unit=150, file=(trim(runname) // "_" // trim(title(n)) // ".txt"), action="write", status="replace")

            do j=1,nbin(n)
                write (150, *) xhis(n,j), hist(n,j), hist(m,j)
            enddo
          close(unit=150)

      end subroutine
