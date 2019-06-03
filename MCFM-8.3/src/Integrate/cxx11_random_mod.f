      module cxx11random
        implicit none

        interface
          subroutine cxx11_init_random(seeds) bind(C, name="cxx11_init_random")
            use ISO_C_BINDING, only: c_int, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: seeds
          end subroutine

          function cxx11_random_number() bind(C, name="cxx11_random_number")
            use ISO_C_BINDING, only: c_double
            implicit none
            real(c_double) :: cxx11_random_number
          end function
        end interface

      end module
