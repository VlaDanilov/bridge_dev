      module mathfun
          use types

      public :: crossprod
      public :: euclideanNorm

      contains

      function euclideanNorm(x)
          implicit none

          real(dp) :: euclideanNorm
          real(dp), intent(in) :: x(:)

          euclideanNorm = sqrt(sum(x**2))

      end function
      
      function crossprod(a,b)
          implicit none

          real(dp), intent(in) :: a(3), b(3)
          real(dp) :: crossprod(3)

          crossprod(1) = a(2)*b(3) - a(3)*b(2)
          crossprod(2) = a(3)*b(1) - a(1)*b(3)
          crossprod(3) = a(1)*b(2) - a(2)*b(1)
      end function

      end module
