      module jettagging
        use types
        implicit none

        include 'mxpart.f'

        integer, save, public :: jetcontent(mxpart)
!$omp threadprivate(jetcontent)

        integer, parameter, public :: bitflag(14) = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192]

      end module
