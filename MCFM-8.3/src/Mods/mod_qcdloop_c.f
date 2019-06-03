      module mod_qcdloop_c
        implicit none

      interface
        subroutine qlcachesize(csize) bind(C,name="qlcachesize")
            use iso_c_binding
            integer(c_int), intent(in) :: csize
        end subroutine

        function cln(x,isig) bind(C,name="cln")
          use iso_c_binding
          complex(c_double_complex), intent(in) :: x
          real(c_double), intent(in) :: isig
          complex(c_double_complex) :: cln
        end function

        function qlzero(x) bind(C,name="qlzero")
          use iso_c_binding
          real(c_double), intent(in) :: x
          logical(c_bool) :: qlzero
        end function

        function qlnonzero(x) bind(C,name="qlnonzero")
          use iso_c_binding
          real(c_double), intent(in) :: x
          logical(c_bool) :: qlnonzero
        end function

        function qli1(m1,mu2,ep) bind(C,name="qli1")
          use iso_c_binding
          real(c_double), intent(in) :: m1,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli1
        end function

        function qli1c(m1,mu2,ep) bind(C,name="qli1c")
          use iso_c_binding
          complex(c_double_complex) :: m1
          real(c_double), intent(in) :: mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli1c
        end function

        function qli2(p1,m1,m2,mu2,ep) bind(C,name="qli2")
          use iso_c_binding
          real(c_double), intent(in) :: p1,m1,m2,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli2
        end function

        function qli2c(p1,m1,m2,mu2,ep) bind(C,name="qli2c")
          use iso_c_binding
          real(c_double), intent(in) :: p1,mu2
          complex(c_double_complex) :: m1,m2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli2c
        end function

        function qli3(p1,p2,p3,m1,m2,m3,mu2,ep) bind(C,name="qli3")
          use iso_c_binding
          real(c_double), intent(in) :: p1,p2,p3,m1,m2,m3,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli3
        end function

        function qli3c(p1,p2,p3,m1,m2,m3,mu2,ep) bind(C,name="qli3c")
          use iso_c_binding
          real(c_double), intent(in) :: p1,p2,p3,mu2
          complex(c_double_complex) :: m1,m2,m3
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli3c
        end function

        function qli4(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
     &      bind(C,name="qli4")
          use iso_c_binding
          real(c_double), intent(in) :: p1,p2,p3,p4,s12,s23
          real(c_double), intent(in) :: m1,m2,m3,m4,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli4
        end function

        function qli4c(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
     &      bind(C,name="qli4c")
          use iso_c_binding
          real(c_double), intent(in) :: p1,p2,p3,p4,s12,s23, mu2
          complex(c_double_complex), intent(in) :: m1,m2,m3,m4
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli4c
        end function

        ! I consider the quad precision interface to be broken:
        ! it probably only works with gfortran since there are no
        ! standard c bindings for quad precision types
        ! (on gfortran there is c_float128 and c_float128_complex)
        ! fortran2008 defines real128, but that is still no interface type

        function qli1q(m1,mu2,ep) bind(C,name="qli1q")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: m1,mu2
          integer(c_int), intent(in) :: ep
          complex*32 :: qli1q
        end function 

        function qli2q(p1,m1,m2,mu2,ep) bind(C,name="qli2q")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: p1,m1,m2,mu2
          integer(c_int), intent(in) :: ep
          complex*32 :: qli2q
        end function 

        function qli3q(p1,p2,p3,m1,m2,m3,mu2,ep)
     &      bind(C,name="qli3q")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: p1,p2,p3,m1,m2,m3,mu2
          integer(c_int), intent(in) :: ep
          complex*32 :: qli3q
        end function

        function qli4q(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
     &      bind(C,name="qli4q")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: p1,p2,p3,p4,s12,s23
          real*16, intent(in) :: m1,m2,m3,m4,mu2
          integer(c_int), intent(in) :: ep
          complex*32 :: qli4q
        end function

        function qli1qc(m1,mu2,ep) bind(C,name="qli1qc")
          use iso_c_binding
          use iso_fortran_env
          complex*32, intent(in) :: m1
          real*16, intent(in) :: mu2
          integer(c_int), intent(in) :: ep
          complex*32 :: qli1qc
        end function

        function qli2qc(p1,m1,m2,mu2,ep) bind(C,name="qli2qc")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: p1,mu2
          complex*32, intent(in) :: m1,m2
          integer(c_int), intent(in) :: ep
          complex*32 :: qli2qc
        end function

        function qli3qc(p1,p2,p3,m1,m2,m3,mu2,ep)
     &      bind(C,name="qli3qc")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: p1,p2,p3,mu2
          complex*32, intent(in) :: m1,m2,m3
          integer(c_int), intent(in) :: ep
          complex*32 :: qli3qc
        end function

        function qli4qc(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
     &      bind(C,name="qli4qc")
          use iso_c_binding
          use iso_fortran_env
          real*16, intent(in) :: p1,p2,p3,p4,s12,s23
          real*16, intent(in) :: mu2
          complex*32, intent(in) :: m1,m2,m3,m4
          integer(c_int), intent(in) :: ep
          complex*32 :: qli4qc
        end function

      end interface

      end module
