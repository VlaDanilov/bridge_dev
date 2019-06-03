       function streal_lightResonant_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_P_SM =
     &     -((propT1267*propW167*propW34*za(jn,jc)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)
     &     *zb(jg,jb))*(za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg
     &     ,je)))/(za(jd,jg)*za(ju,jg)))
       end function streal_lightResonant_MMMM_P_SM

       function streal_lightResonant_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_M_SM =
     &     -((propT1267*propW167*propW34*za(jn,jc)*zb(jb,ju)*(za(jb,jd)*zb(jd,ju)*z
     &     b(je,jb) + za(ju,jd)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*zb(jd,ju)*zb(jg,je)
     &      - za(jd,jg)*zb(jd,je)*zb(jg,ju) + za(jb,jg)*zb(je,jb)*zb(jg,ju) + za(ju
     &     ,jg)*zb(je,ju)*zb(jg,ju)))/(zb(jg,jd)*zb(jg,ju)))
       end function streal_lightResonant_MMMM_M_SM

       function streal_heavyResonant_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MMMM_P_SM =
     &     (propW16*propW34*za(jn,jc)*(propT1267*za(jc,jg)*(-(za(jc,jd)*zb(jc,je)) 
     &     + za(jn,jd)*zb(je,jn))*(za(jb,jn)*zb(jb,ju) - za(jn,jg)*zb(jg,ju)) + pro
     &     pT126*za(jb,jg)*zb(jb,ju)*(za(jb,jd)*za(jn,jc)*(-zb(je,jb) + propT1267*z
     &     a(jc,jg)*zb(jc,je)*zb(jg,jb)) + mtsq*propT1267*za(jc,jg)*za(jn,jd)*zb(jg
     &     ,je) + za(jn,jc)*za(ju,jd)*(-zb(je,ju) + propT1267*za(jc,jg)*zb(jc,je)*z
     &     b(jg,ju)))))/(za(jb,jg)*za(jc,jg)*za(jn,jg))
       end function streal_heavyResonant_MMMM_P_SM

       function streal_heavyResonant_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MMMM_M_SM =
     &     (propT126*propW16*propW34*zb(jb,ju)*(-(zb(je,jb)*(za(jb,jd)*(za(jn,jc)*z
     &     b(jc,jb) + za(jn,jg)*zb(jg,jb)) - mtsq*propT1267*za(jd,jg)*za(jn,jc)*zb(
     &     jg,jc))) - za(ju,jd)*(za(jn,jg)*zb(je,ju)*zb(jg,jb) + za(jn,jc)*(zb(jc,j
     &     b)*zb(je,ju) + propT1267*zb(jb,ju)*(za(jc,jg)*zb(jc,je) - za(jn,jg)*zb(j
     &     e,jn))*zb(jg,jc)))))/(zb(jg,jb)*zb(jg,jc))
       end function streal_heavyResonant_MMMM_M_SM

       function streal_lightWWG_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_MMMM_P_SM =
     &     -(gb**2*propW167*propW34*(za(jn,jd)*za(ju,jc)*za(ju,jd)*zb(jb,ju)*zb(je,
     &     ju) - za(je,jd)*za(jn,jc)*za(ju,jd)*zb(je,jb)*zb(je,ju) + za(je,jc)*za(j
     &     n,jd)*za(ju,jd)*zb(je,jb)*zb(je,ju) + za(jc,jg)*za(jn,jd)*za(ju,jd)*zb(j
     &     e,ju)*zb(jg,jb) + za(jb,jd)*za(jn,jc)*zb(je,jb)*(za(ju,jd)*zb(jb,ju) + z
     &     a(jd,jg)*zb(jg,jb)) + za(jd,jg)*za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(jg,je) 
     &     - za(jd,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(jg,je) + za(jd,jg)*za(je,jc
     &     )*za(jn,jd)*zb(je,jb)*zb(jg,je) + za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jb
     &     )*zb(jg,je) + za(jc,jd)*(za(jn,jc)*(za(ju,jd)*(zb(jb,ju)*zb(jc,je) + zb(
     &     jc,ju)*zb(je,jb)) + za(jd,jg)*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc)
     &     )) + (za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*(za(jb,jn)*zb(je,jb) - 
     &     za(ju,jn)*zb(je,ju) - za(jn,jg)*zb(jg,je)) + za(jn,jd)*(za(ju,jd)*(-(zb(
     &     jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + za(jd,jg)*(-(zb(jd,je)*zb(jg,
     &     jb)) + zb(jd,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*
     &     zb(jg,jn) + za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jn,jb) + za(jd,jg
     &     )*za(jn,jc)*za(jn,jd)*zb(jg,je)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*za(ju,jd
     &     )*zb(je,jb)*zb(jn,ju)))/(3._dp*ecossin**2*za(jb,jc)*za(jd,jg)*za(ju,jg)*zb(
     &     jc,jb))
       end function streal_lightWWG_MMMM_P_SM

       function streal_lightWWG_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_MMMM_M_SM =
     &     -(gb**2*propW167*propW34*(za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(jd,ju)*zb(je,
     &     jb) + za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) - za(je,jd)*za(j
     &     n,jc)*zb(jd,ju)*zb(je,jb)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jd,ju)*zb(j
     &     e,jb)*zb(je,ju) + za(jc,jg)*za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) + za
     &     (jc,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,ju) - za(jc,jg)*za(jn,jd)*zb
     &     (jb,ju)*zb(jd,je)*zb(jg,ju) + za(jb,jn)*za(jc,jg)*zb(jb,ju)*zb(je,jb)*zb
     &     (jg,ju) + za(jb,jg)*za(jn,jc)*zb(jb,ju)*zb(je,jb)*zb(jg,ju) + za(jc,jg)*
     &     za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,ju) + za(jn,jg)*za(ju,jc)*zb(jb,ju)*
     &     zb(je,ju)*zb(jg,ju) - za(jc,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) 
     &     - za(je,jg)*za(jn,jc)*zb(je,jb)*zb(je,ju)*zb(jg,ju) + za(je,jc)*za(jn,jg
     &     )*zb(je,jb)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jb
     &     )*zb(jg,ju) - za(jc,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je)*zb(jg,ju) + za(jc,
     &     jd)*(za(jb,jn)*zb(jb,ju)*zb(jd,ju)*zb(je,jb) + za(jn,jc)*zb(jd,ju)*(zb(j
     &     b,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) - za(ju,jn)*zb(jb,ju)*zb(jd,ju)*z
     &     b(je,ju) + za(jn,jd)*zb(jd,ju)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je
     &     ,ju)) - za(jn,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,jb)*zb
     &     (je,ju)*zb(jg,ju)) + za(jn,jc)*za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jn,jb) +
     &      za(jn,jc)*za(jn,jg)*zb(je,ju)*zb(jg,ju)*zb(jn,jb) - za(jn,jc)*za(jn,jd)
     &     *zb(jd,ju)*zb(je,jb)*zb(jn,ju) - za(jn,jc)*za(jn,jg)*zb(je,jb)*zb(jg,ju)
     &     *zb(jn,ju)))/(3._dp*ecossin**2*za(jb,jc)*zb(jc,jb)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWG_MMMM_M_SM

       function streal_lightWWG_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_PPMM_P_SM =
     &     -(gb**2*propW167*propW34*(za(jn,jd)*(za(jb,je)*zb(jc,je) + za(ju,jb)*zb(
     &     jc,ju) + za(jb,jg)*zb(jg,jc))*(za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)
     &     ) + za(jb,jd)*(-(za(jn,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,je)) - za(ju,jd)*za
     &     (ju,jn)*zb(jc,ju)*zb(je,ju) + za(jn,jd)*za(ju,jd)*zb(jd,jc)*zb(je,ju) - 
     &     za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(je,ju)*
     &     zb(jg,jc) + za(jn,jc)*zb(jc,je)*(za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,j
     &     c)) + za(jb,jn)*(za(ju,jd)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) +
     &      za(jd,jg)*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc))) - za(jn,jg)*za(j
     &     u,jd)*zb(jc,ju)*zb(jg,je) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*zb(jg,je) - za
     &     (jd,jg)*za(jn,jg)*zb(jg,jc)*zb(jg,je)) + za(jb,jn)*(za(jc,jd)*zb(jc,je)*
     &     (za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc)) - za(je,jd)*zb(jc,je)*(za(ju
     &     ,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + za(jn,jd)*(za(jd,jg)*(zb(jc,jn)*
     &     zb(jg,je) - zb(jc,je)*zb(jg,jn)) + za(ju,jd)*(zb(jc,jn)*zb(je,ju) - zb(j
     &     c,je)*zb(jn,ju))))))/(3._dp*ecossin**2*za(jb,jc)*za(jd,jg)*za(ju,jg)*zb(jc,
     &     jb))
       end function streal_lightWWG_PPMM_P_SM

       function streal_lightWWG_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_PPMM_M_SM =
     &     -(gb**2*propW167*propW34*(za(jb,je)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(je,
     &     ju) + za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jb,jg)*za(j
     &     n,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jn,jc)*zb(jc,je)*zb(j
     &     c,ju)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,ju) + za
     &     (jb,je)*za(jn,jg)*zb(jc,je)*zb(je,ju)*zb(jg,ju) + za(jn,jg)*za(ju,jb)*zb
     &     (jc,ju)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb
     &     (jg,ju) + za(jb,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*
     &     za(jn,jg)*zb(jc,ju)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(za(jn,jc)*zb(jc,je)
     &     *zb(jc,ju)*zb(jd,ju) - za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jd,ju) + za(jb,j
     &     n)*zb(jd,ju)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) - za(ju,jn)*zb(
     &     jc,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,jc)*zb(jd,ju)*zb(je,ju) - z
     &     a(jn,jg)*zb(jc,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,jc)*zb(je,ju)*z
     &     b(jg,ju)) + za(jb,jn)*(za(jc,jd)*zb(jc,je)*zb(jc,ju)*zb(jd,ju) - za(je,j
     &     d)*zb(jc,je)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jc,jn)*zb(jd,ju)*zb(je,j
     &     u) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,ju) + za(jc,jg)*zb(jc,je)*zb(jc
     &     ,ju)*zb(jg,ju) + za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,ju) - za(je,jg)*zb(
     &     jc,je)*zb(je,ju)*zb(jg,ju) + za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,ju) - z
     &     a(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(jn,ju) - za(jn,jg)*zb(jc,je)*zb(jg,ju)*z
     &     b(jn,ju))))/(3._dp*ecossin**2*za(jb,jc)*zb(jc,jb)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWG_PPMM_M_SM

       function streal_heavyWWG_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_MMMM_P_SM =
     &     -(gb**2*propW16*propW34*(za(jb,jc)*(za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(je,
     &     jb) + za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(j
     &     e,jb)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju) + za(jc,jd)*(z
     &     a(jn,jc)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(-(zb(j
     &     b,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(za(jb,jn)*zb(je,jb)
     &      - za(ju,jn)*zb(je,ju) + za(jn,jg)*zb(jg,je))) - za(jd,jg)*za(jn,jc)*zb(
     &     je,jb)*zb(jg,ju) + za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jn,jb) - za(jn,jc)*z
     &     a(jn,jd)*zb(je,jb)*zb(jn,ju)) + za(jc,jg)*(za(jb,jd)*za(jn,jc)*zb(jb,ju)
     &     *zb(jg,je) - za(je,jd)*za(jn,jc)*zb(je,ju)*zb(jg,je) + za(je,jc)*za(jn,j
     &     d)*zb(je,ju)*zb(jg,je) + za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jn) - za(jn
     &     ,jd)*za(ju,jc)*zb(je,ju)*zb(jg,ju) - za(jd,jg)*za(jn,jc)*zb(jg,je)*zb(jg
     &     ,ju) + za(jc,jd)*(-((za(jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) + za(jn,j
     &     g)*zb(jg,je))*zb(jg,ju)) + za(jn,jc)*(zb(jc,ju)*zb(jg,je) - zb(jc,je)*zb
     &     (jg,ju)) + za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))) - za(j
     &     n,jc)*za(jn,jd)*zb(jg,je)*zb(jn,ju))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,j
     &     g) + s(jc,jg))*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWG_MMMM_P_SM

       function streal_heavyWWG_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_MMMM_M_SM =
     &     -(gb**2*propW16*propW34*(za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(jc,jb)*zb(je,j
     &     u) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(je,jc)*za(jn
     &     ,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) - za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc
     &     ,je)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jb) - za(
     &     jb,jn)*za(jd,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jb) + za(jn,jd)*za(ju,jg)*zb(
     &     jb,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(
     &     jg,jb) - za(jd,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jb) + za(je,jg)*z
     &     a(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*zb(je,jb)*z
     &     b(je,ju)*zb(jg,jb) + za(jb,jd)*zb(jb,ju)*zb(je,jb)*(za(jn,jc)*zb(jc,jb) 
     &     + za(jn,jg)*zb(jg,jb)) - za(jd,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,jb)*zb(jg,j
     &     e) + za(jc,jd)*(za(jb,jn)*zb(jb,ju)*zb(jc,jb)*zb(je,jb) + za(jn,jc)*zb(j
     &     c,jb)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) - za(ju,jn)*zb(jb,ju)*
     &     zb(jc,jb)*zb(je,ju) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,ju)*zb(jd,je)) + zb(j
     &     d,jb)*zb(je,ju)) + za(jn,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(jn,jg)*z
     &     b(jb,ju)*zb(jc,jb)*zb(jg,je)) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(je,jb)*
     &     zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju) + za(jn,jc
     &     )*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) + za(jn,jd)*za(jn,jg)*zb(je,ju
     &     )*zb(jg,jb)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju
     &     ) - za(jn,jd)*za(jn,jg)*zb(je,jb)*zb(jg,jb)*zb(jn,ju)))/(3._dp*ecossin**2*(
     &     s(jb,jc) + s(jb,jg) + s(jc,jg))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWG_MMMM_M_SM

       function streal_heavyWWG_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_PPMM_P_SM =
     &     -(gb**2*propW16*propW34*(za(jb,jc)*(za(jn,jd)*(za(jb,je)*zb(jc,je) + za(
     &     ju,jb)*zb(jc,ju))*zb(je,ju) + za(jb,jd)*(za(jn,jc)*zb(jc,je)*zb(jc,ju) -
     &      za(jn,jd)*zb(jc,ju)*zb(jd,je) + za(jb,jn)*(zb(jb,ju)*zb(jc,je) + zb(jc,
     &     ju)*zb(je,jb)) - za(ju,jn)*zb(jc,ju)*zb(je,ju) + za(jn,jd)*zb(jd,jc)*zb(
     &     je,ju) + za(jn,jg)*zb(jc,ju)*zb(jg,je)) + za(jb,jn)*(za(jc,jd)*zb(jc,je)
     &     *zb(jc,ju) - za(je,jd)*zb(jc,je)*zb(je,ju) + za(jn,jd)*zb(jc,jn)*zb(je,j
     &     u) - za(jd,jg)*zb(jc,je)*zb(jg,ju) - za(jn,jd)*zb(jc,je)*zb(jn,ju))) + z
     &     a(jb,jg)*(za(jn,jd)*zb(je,ju)*(za(jb,je)*zb(jg,je) + za(ju,jb)*zb(jg,ju)
     &     ) + za(jb,jd)*((za(jn,jc)*zb(jc,je) - za(ju,jn)*zb(je,ju) + za(jn,jg)*zb
     &     (jg,je))*zb(jg,ju) - za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju
     &     )) + za(jb,jn)*(zb(jb,ju)*zb(jg,je) + zb(je,jb)*zb(jg,ju))) + za(jb,jn)*
     &     (za(jc,jd)*zb(jc,ju)*zb(jg,je) - za(je,jd)*zb(je,ju)*zb(jg,je) + za(jn,j
     &     d)*zb(je,ju)*zb(jg,jn) - za(jd,jg)*zb(jg,je)*zb(jg,ju) - za(jn,jd)*zb(jg
     &     ,je)*zb(jn,ju)))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*za(jb
     &     ,jg)*za(jc,jg))
       end function streal_heavyWWG_PPMM_P_SM

       function streal_heavyWWG_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_PPMM_M_SM =
     &     -(gb**2*propW16*propW34*(za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,j
     &     u) + za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jd,jg)*za(jn
     &     ,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jc) + za(jc,jd)*za(jn,jg)*zb(jc,je)*zb(jc
     &     ,ju)*zb(jg,jc) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jc) + za(
     &     je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jc) - za(je,jd)*za(jn,jg)*zb(
     &     jc,je)*zb(je,ju)*zb(jg,jc) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(
     &     jg,jc) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) - za(jd,jg)*z
     &     a(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*z
     &     b(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc)*zb(jg,je) +
     &      za(jb,jd)*(za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jc,ju) - za(jn,jd)*zb(jc,jb
     &     )*zb(jc,ju)*zb(jd,je) + za(jb,jn)*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) + zb(jc
     &     ,ju)*zb(je,jb)) - za(ju,jn)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jn,jd)*zb
     &     (jc,jb)*zb(jd,jc)*zb(je,ju) + za(jn,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) + 
     &     za(jn,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je)) - za(jd,jg)*za(jn,jg)*zb(jc,je)
     &     *zb(jg,jc)*zb(jg,ju) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jc)*zb(jn,ju)
     &      + za(jb,jn)*(za(jc,jd)*zb(jc,jb)*zb(jc,je)*zb(jc,ju) - za(je,jd)*zb(jc,
     &     jb)*zb(jc,je)*zb(je,ju) + za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) + za(j
     &     d,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jd,jg)*zb(jc,jb)*zb(jc,je)*zb(j
     &     g,ju) - za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(jn,ju))))/(3._dp*ecossin**2*(s(jb,
     &     jc) + s(jb,jg) + s(jc,jg))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWG_PPMM_M_SM

       function streal_lightWWZ_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_MMMM_P_SM =
     &     -((gb**2 + 3*gw**2)*propW167*propW34*propZ25*(za(jn,jd)*za(ju,jc)*za(ju,
     &     jd)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*za(ju,jd)*zb(je,jb)*zb(je,
     &     ju) + za(je,jc)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(je,ju) + za(jc,jg)*za(j
     &     n,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jb) + za(jb,jd)*za(jn,jc)*zb(je,jb)*(za(
     &     ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)) + za(jd,jg)*za(jn,jd)*za(ju,jc)*
     &     zb(jb,ju)*zb(jg,je) - za(jd,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(jg,je) 
     &     + za(jd,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(jg,je) + za(jc,jg)*za(jd,jg
     &     )*za(jn,jd)*zb(jg,jb)*zb(jg,je) + za(jc,jd)*(za(jn,jc)*(za(ju,jd)*(zb(jb
     &     ,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) + za(jd,jg)*(zb(jc,je)*zb(jg,jb) +
     &      zb(je,jb)*zb(jg,jc))) + (za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*(za
     &     (jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) - za(jn,jg)*zb(jg,je)) + za(jn,j
     &     d)*(za(ju,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + za(jd,jg)
     &     *(-(zb(jd,je)*zb(jg,jb)) + zb(jd,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*
     &     za(jn,jd)*zb(je,jb)*zb(jg,jn) + za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(je,ju)*
     &     zb(jn,jb) + za(jd,jg)*za(jn,jc)*za(jn,jd)*zb(jg,je)*zb(jn,jb) - za(jn,jc
     &     )*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(jn,ju)))/(6._dp*ecossin**2*za(jd,jg)*za(
     &     ju,jg))
       end function streal_lightWWZ_MMMM_P_SM

       function streal_lightWWZ_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_MMMM_M_SM =
     &     -((gb**2 + 3*gw**2)*propW167*propW34*propZ25*(za(jb,jd)*za(jn,jc)*zb(jb,
     &     ju)*zb(jd,ju)*zb(je,jb) + za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(jd,ju)*zb(je,
     &     ju) - za(je,jd)*za(jn,jc)*zb(jd,ju)*zb(je,jb)*zb(je,ju) + za(je,jc)*za(j
     &     n,jd)*zb(jd,ju)*zb(je,jb)*zb(je,ju) + za(jc,jg)*za(jn,jd)*zb(jd,ju)*zb(j
     &     e,ju)*zb(jg,jb) + za(jc,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,ju) - za
     &     (jc,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,ju) + za(jb,jn)*za(jc,jg)*zb
     &     (jb,ju)*zb(je,jb)*zb(jg,ju) + za(jb,jg)*za(jn,jc)*zb(jb,ju)*zb(je,jb)*zb
     &     (jg,ju) + za(jc,jg)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,ju) + za(jn,jg)*
     &     za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za(jc,jg)*za(ju,jn)*zb(jb,ju)*
     &     zb(je,ju)*zb(jg,ju) - za(je,jg)*za(jn,jc)*zb(je,jb)*zb(je,ju)*zb(jg,ju) 
     &     + za(je,jc)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*za(jn,jg
     &     )*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jc,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je
     &     )*zb(jg,ju) + za(jc,jd)*(za(jb,jn)*zb(jb,ju)*zb(jd,ju)*zb(je,jb) + za(jn
     &     ,jc)*zb(jd,ju)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) - za(ju,jn)*z
     &     b(jb,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)*(-(zb(jb,ju)*zb(jd,je
     &     )) + zb(jd,jb)*zb(je,ju)) - za(jn,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za
     &     (jn,jg)*zb(jd,jb)*zb(je,ju)*zb(jg,ju)) + za(jn,jc)*za(jn,jd)*zb(jd,ju)*z
     &     b(je,ju)*zb(jn,jb) + za(jn,jc)*za(jn,jg)*zb(je,ju)*zb(jg,ju)*zb(jn,jb) -
     &      za(jn,jc)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(jn,ju) - za(jn,jc)*za(jn,jg)
     &     *zb(je,jb)*zb(jg,ju)*zb(jn,ju)))/(6._dp*ecossin**2*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWZ_MMMM_M_SM

       function streal_lightWWZ_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_PPMM_P_SM =
     &     (gb**2*propW167*propW34*propZ25*(za(jn,jd)*(za(jb,je)*zb(jc,je) + za(ju,
     &     jb)*zb(jc,ju) + za(jb,jg)*zb(jg,jc))*(za(ju,jd)*zb(je,ju) + za(jd,jg)*zb
     &     (jg,je)) + za(jb,jd)*(-(za(jn,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,je)) - za(ju
     &     ,jd)*za(ju,jn)*zb(jc,ju)*zb(je,ju) + za(jn,jd)*za(ju,jd)*zb(jd,jc)*zb(je
     &     ,ju) - za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(
     &     je,ju)*zb(jg,jc) + za(jn,jc)*zb(jc,je)*(za(ju,jd)*zb(jc,ju) + za(jd,jg)*
     &     zb(jg,jc)) + za(jb,jn)*(za(ju,jd)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je
     &     ,jb)) + za(jd,jg)*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc))) - za(jn,j
     &     g)*za(ju,jd)*zb(jc,ju)*zb(jg,je) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*zb(jg,j
     &     e) - za(jd,jg)*za(jn,jg)*zb(jg,jc)*zb(jg,je)) + za(jb,jn)*(za(jc,jd)*zb(
     &     jc,je)*(za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc)) - za(je,jd)*zb(jc,je)
     &     *(za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + za(jn,jd)*(za(jd,jg)*(zb(
     &     jc,jn)*zb(jg,je) - zb(jc,je)*zb(jg,jn)) + za(ju,jd)*(zb(jc,jn)*zb(je,ju)
     &      - zb(jc,je)*zb(jn,ju))))))/(3._dp*ecossin**2*za(jd,jg)*za(ju,jg))
       end function streal_lightWWZ_PPMM_P_SM

       function streal_lightWWZ_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_PPMM_M_SM =
     &     (gb**2*propW167*propW34*propZ25*(za(jb,je)*za(jn,jd)*zb(jc,je)*zb(jd,ju)
     &     *zb(je,ju) + za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jb,j
     &     g)*za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jn,jc)*zb(jc,j
     &     e)*zb(jc,ju)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,j
     &     u) + za(jb,je)*za(jn,jg)*zb(jc,je)*zb(je,ju)*zb(jg,ju) + za(jn,jg)*za(ju
     &     ,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(ju,jn)*zb(jc,ju)*zb(je
     &     ,ju)*zb(jg,ju) + za(jb,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) - za(
     &     jb,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(za(jn,jc)*zb
     &     (jc,je)*zb(jc,ju)*zb(jd,ju) - za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jd,ju) + 
     &     za(jb,jn)*zb(jd,ju)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) - za(ju,
     &     jn)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,jc)*zb(jd,ju)*zb(je,
     &     ju) - za(jn,jg)*zb(jc,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,jc)*zb(j
     &     e,ju)*zb(jg,ju)) + za(jb,jn)*(za(jc,jd)*zb(jc,je)*zb(jc,ju)*zb(jd,ju) - 
     &     za(je,jd)*zb(jc,je)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jc,jn)*zb(jd,ju)*
     &     zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,ju) + za(jc,jg)*zb(jc,je
     &     )*zb(jc,ju)*zb(jg,ju) + za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,ju) - za(je,
     &     jg)*zb(jc,je)*zb(je,ju)*zb(jg,ju) + za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,
     &     ju) - za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(jn,ju) - za(jn,jg)*zb(jc,je)*zb(j
     &     g,ju)*zb(jn,ju))))/(3._dp*ecossin**2*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWZ_PPMM_M_SM

       function streal_heavyWWZ_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_MMMM_P_SM =
     &     -((gb**2 + 3*gw**2)*propW16*propW34*propZ257*(za(jb,jc)*(za(jb,jd)*za(jn
     &     ,jc)*zb(jb,ju)*zb(je,jb) + za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju) - za(
     &     je,jd)*za(jn,jc)*zb(je,jb)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(je,jb)*zb(
     &     je,ju) + za(jc,jd)*(za(jn,jc)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)
     &     ) + za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)
     &     *(za(jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) + za(jn,jg)*zb(jg,je))) - za
     &     (jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,ju) + za(jn,jc)*za(jn,jd)*zb(je,ju)*zb
     &     (jn,jb) - za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jn,ju)) + za(jc,jg)*(za(jb,jd
     &     )*za(jn,jc)*zb(jb,ju)*zb(jg,je) - za(je,jd)*za(jn,jc)*zb(je,ju)*zb(jg,je
     &     ) + za(je,jc)*za(jn,jd)*zb(je,ju)*zb(jg,je) + za(jn,jc)*za(jn,jd)*zb(je,
     &     ju)*zb(jg,jn) - za(jn,jd)*za(ju,jc)*zb(je,ju)*zb(jg,ju) - za(jd,jg)*za(j
     &     n,jc)*zb(jg,je)*zb(jg,ju) + za(jc,jd)*(-((za(jb,jn)*zb(je,jb) - za(ju,jn
     &     )*zb(je,ju) + za(jn,jg)*zb(jg,je))*zb(jg,ju)) + za(jn,jc)*(zb(jc,ju)*zb(
     &     jg,je) - zb(jc,je)*zb(jg,ju)) + za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,j
     &     e)*zb(jg,ju))) - za(jn,jc)*za(jn,jd)*zb(jg,je)*zb(jn,ju))))/(6._dp*ecossin*
     &     *2*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWZ_MMMM_P_SM

       function streal_heavyWWZ_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_MMMM_M_SM =
     &     -((gb**2 + 3*gw**2)*propW16*propW34*propZ257*(za(jn,jd)*za(ju,jc)*zb(jb,
     &     ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,jb)*zb(je,
     &     ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) - za(jd,jg)*za(j
     &     n,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(j
     &     d,je)*zb(jg,jb) - za(jb,jn)*za(jd,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jb) + za
     &     (jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(ju,jn)*zb
     &     (jb,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb
     &     (jg,jb) + za(je,jg)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(je,jd)*
     &     za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,jb) + za(jb,jd)*zb(jb,ju)*zb(je,jb)*
     &     (za(jn,jc)*zb(jc,jb) + za(jn,jg)*zb(jg,jb)) - za(jd,jg)*za(jn,jg)*zb(jb,
     &     ju)*zb(jg,jb)*zb(jg,je) + za(jc,jd)*(za(jb,jn)*zb(jb,ju)*zb(jc,jb)*zb(je
     &     ,jb) + za(jn,jc)*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) + zb(jc,ju)*zb(je,jb)) -
     &      za(ju,jn)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,
     &     ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + za(jn,jg)*zb(jc,ju)*zb(je,jb)*zb
     &     (jg,jb) + za(jn,jg)*zb(jb,ju)*zb(jc,jb)*zb(jg,je)) - za(jd,jg)*za(jn,jc)
     &     *zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jb)
     &     *zb(jg,ju) + za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) + za(jn,j
     &     d)*za(jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,j
     &     b)*zb(je,jb)*zb(jn,ju) - za(jn,jd)*za(jn,jg)*zb(je,jb)*zb(jg,jb)*zb(jn,j
     &     u)))/(6._dp*ecossin**2*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWZ_MMMM_M_SM

       function streal_heavyWWZ_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_PPMM_P_SM =
     &     (gb**2*propW16*propW34*propZ257*(za(jb,jc)*(za(jn,jd)*(za(jb,je)*zb(jc,j
     &     e) + za(ju,jb)*zb(jc,ju))*zb(je,ju) + za(jb,jd)*(za(jn,jc)*zb(jc,je)*zb(
     &     jc,ju) - za(jn,jd)*zb(jc,ju)*zb(jd,je) + za(jb,jn)*(zb(jb,ju)*zb(jc,je) 
     &     + zb(jc,ju)*zb(je,jb)) - za(ju,jn)*zb(jc,ju)*zb(je,ju) + za(jn,jd)*zb(jd
     &     ,jc)*zb(je,ju) + za(jn,jg)*zb(jc,ju)*zb(jg,je)) + za(jb,jn)*(za(jc,jd)*z
     &     b(jc,je)*zb(jc,ju) - za(je,jd)*zb(jc,je)*zb(je,ju) + za(jn,jd)*zb(jc,jn)
     &     *zb(je,ju) - za(jd,jg)*zb(jc,je)*zb(jg,ju) - za(jn,jd)*zb(jc,je)*zb(jn,j
     &     u))) + za(jb,jg)*(za(jn,jd)*zb(je,ju)*(za(jb,je)*zb(jg,je) + za(ju,jb)*z
     &     b(jg,ju)) + za(jb,jd)*((za(jn,jc)*zb(jc,je) - za(ju,jn)*zb(je,ju) + za(j
     &     n,jg)*zb(jg,je))*zb(jg,ju) - za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*
     &     zb(jg,ju)) + za(jb,jn)*(zb(jb,ju)*zb(jg,je) + zb(je,jb)*zb(jg,ju))) + za
     &     (jb,jn)*(za(jc,jd)*zb(jc,ju)*zb(jg,je) - za(je,jd)*zb(je,ju)*zb(jg,je) +
     &      za(jn,jd)*zb(je,ju)*zb(jg,jn) - za(jd,jg)*zb(jg,je)*zb(jg,ju) - za(jn,j
     &     d)*zb(jg,je)*zb(jn,ju)))))/(3._dp*ecossin**2*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWZ_PPMM_P_SM

       function streal_heavyWWZ_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_PPMM_M_SM =
     &     (gb**2*propW16*propW34*propZ257*(za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(jc,je)
     &     *zb(je,ju) + za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jd,j
     &     g)*za(jn,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jc) + za(jc,jd)*za(jn,jg)*zb(jc,j
     &     e)*zb(jc,ju)*zb(jg,jc) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,j
     &     c) + za(je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jc) - za(je,jd)*za(jn
     &     ,jg)*zb(jc,je)*zb(je,ju)*zb(jg,jc) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je
     &     ,ju)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) - za(
     &     jd,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*zb(
     &     jd,jc)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc)*zb(
     &     jg,je) + za(jb,jd)*(za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jc,ju) - za(jn,jd)*
     &     zb(jc,jb)*zb(jc,ju)*zb(jd,je) + za(jb,jn)*zb(jc,jb)*(zb(jb,ju)*zb(jc,je)
     &      + zb(jc,ju)*zb(je,jb)) - za(ju,jn)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(j
     &     n,jd)*zb(jc,jb)*zb(jd,jc)*zb(je,ju) + za(jn,jg)*zb(jb,ju)*zb(jc,je)*zb(j
     &     g,jc) + za(jn,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je)) - za(jd,jg)*za(jn,jg)*z
     &     b(jc,je)*zb(jg,jc)*zb(jg,ju) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jc)*z
     &     b(jn,ju) + za(jb,jn)*(za(jc,jd)*zb(jc,jb)*zb(jc,je)*zb(jc,ju) - za(je,jd
     &     )*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju
     &     ) + za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jd,jg)*zb(jc,jb)*zb(jc,
     &     je)*zb(jg,ju) - za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(jn,ju))))/(3._dp*ecossin**
     &     2*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWZ_PPMM_M_SM

       function streal_lightGR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_MMMM_P_SM =
     &     (4*gb**2*propW34*za(jn,jd)*(za(jc,jd)*zb(jd,je)*(za(ju,jd)*zb(jb,ju) + z
     &     a(jd,jg)*zb(jg,jb)) + zb(je,jn)*(za(jn,jc)*(za(ju,jd)*zb(jb,ju) + za(jd,
     &     jg)*zb(jg,jb)) + (za(jn,jd)*za(ju,jg)*zb(jb,ju)*(za(jb,jc)*zb(jg,jb) + z
     &     a(ju,jc)*zb(jg,ju)))/(s(jb,jc) + s(jb,ju) + s(jc,ju)))))/(9._dp*ecossin**2*
     &     (s(jd,je) + s(jd,jn) + s(je,jn))*za(jb,jc)*za(jd,jg)*za(ju,jg)*zb(jc,jb)
     &     )
       end function streal_lightGR_MMMM_P_SM

       function streal_lightGR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_MMMM_M_SM =
     &     (-4*gb**2*propW34*zb(jb,ju)*(za(ju,jc)*zb(je,ju)*(za(jn,jd)*zb(jd,ju) + 
     &     za(jn,jg)*zb(jg,ju)) + za(jb,jc)*(za(jn,jd)*(zb(jd,ju)*zb(je,jb) + (zb(j
     &     b,ju)*(-(za(jd,jg)*zb(jd,je)) + za(jn,jg)*zb(je,jn))*zb(jg,jd))/(s(jd,je
     &     ) + s(jd,jn) + s(je,jn))) + za(jn,jg)*zb(je,jb)*zb(jg,ju))))/(9._dp*ecossin
     &     **2*(s(jb,jc) + s(jb,ju) + s(jc,ju))*za(jb,jc)*zb(jc,jb)*zb(jg,jd)*zb(jg
     &     ,ju))
       end function streal_lightGR_MMMM_M_SM

       function streal_lightGR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_PPMM_P_SM =
     &     (4*gb**2*propW34*za(jn,jd)*(za(jb,jd)*zb(jd,je)*(za(ju,jd)*zb(jc,ju) + z
     &     a(jd,jg)*zb(jg,jc)) + zb(je,jn)*(-(za(jb,jn)*(za(ju,jd)*zb(jc,ju) + za(j
     &     d,jg)*zb(jg,jc))) + (za(jn,jd)*za(ju,jg)*zb(jc,ju)*(-(za(jb,jc)*zb(jg,jc
     &     )) + za(ju,jb)*zb(jg,ju)))/(s(jb,jc) + s(jb,ju) + s(jc,ju)))))/(9._dp*ecoss
     &     in**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(jb,jc)*za(jd,jg)*za(ju,jg)*zb(
     &     jc,jb))
       end function streal_lightGR_PPMM_P_SM

       function streal_lightGR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_PPMM_M_SM =
     &     (-4*gb**2*propW34*zb(jc,ju)*(za(ju,jb)*zb(je,ju)*(za(jn,jd)*zb(jd,ju) + 
     &     za(jn,jg)*zb(jg,ju)) + za(jb,jc)*(za(jn,jd)*(zb(jc,je)*zb(jd,ju) + (zb(j
     &     c,ju)*(za(jd,jg)*zb(jd,je) - za(jn,jg)*zb(je,jn))*zb(jg,jd))/(s(jd,je) +
     &      s(jd,jn) + s(je,jn))) + za(jn,jg)*zb(jc,je)*zb(jg,ju))))/(9._dp*ecossin**2
     &     *(s(jb,jc) + s(jb,ju) + s(jc,ju))*za(jb,jc)*zb(jc,jb)*zb(jg,jd)*zb(jg,ju
     &     ))
       end function streal_lightGR_PPMM_M_SM

       function streal_lightGL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_MMMM_P_SM =
     &     (-2*gb**2*propW34*za(jc,jd)*(za(jn,jc)*zb(jc,jb)*(za(ju,jd)*zb(je,ju) + 
     &     za(jd,jg)*zb(jg,je)) + za(jn,jd)*zb(jd,jb)*(za(ju,jd)*zb(je,ju) + za(jd,
     &     jg)*zb(jg,je)) + (za(jc,jd)*za(ju,jg)*zb(jc,jb)*zb(je,ju)*(za(jn,je)*zb(
     &     jg,je) - za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(je,ju) + s(jn,ju))))/(9._dp*ec
     &     ossin**2*(s(jb,jc) + s(jb,jd) + s(jc,jd))*za(jb,jc)*za(jd,jg)*za(ju,jg)*
     &     zb(jc,jb))
       end function streal_lightGL_MMMM_P_SM

       function streal_lightGL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_MMMM_M_SM =
     &     (2*gb**2*propW34*zb(je,ju)*(za(ju,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,ju) + z
     &     a(jc,jg)*zb(jg,ju)) + za(jn,je)*(za(jc,jd)*(zb(jd,ju)*zb(je,jb) + ((za(j
     &     c,jg)*zb(jc,jb) + za(jd,jg)*zb(jd,jb))*zb(je,ju)*zb(jg,jd))/(s(jb,jc) + 
     &     s(jb,jd) + s(jc,jd))) + za(jc,jg)*zb(je,jb)*zb(jg,ju))))/(9._dp*ecossin**2*
     &     (s(je,jn) + s(je,ju) + s(jn,ju))*za(jb,jc)*zb(jc,jb)*zb(jg,jd)*zb(jg,ju)
     &     )
       end function streal_lightGL_MMMM_M_SM

       function streal_lightGL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_PPMM_P_SM =
     &     (-2*gb**2*propW34*za(jb,jd)*(za(jb,jn)*zb(jc,jb)*(za(ju,jd)*zb(je,ju) + 
     &     za(jd,jg)*zb(jg,je)) + za(jn,jd)*zb(jd,jc)*(za(ju,jd)*zb(je,ju) + za(jd,
     &     jg)*zb(jg,je)) + (za(jb,jd)*za(ju,jg)*zb(jc,jb)*zb(je,ju)*(-(za(jn,je)*z
     &     b(jg,je)) + za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(je,ju) + s(jn,ju))))/(9.
     &     *ecossin**2*(s(jb,jc) + s(jb,jd) + s(jc,jd))*za(jb,jc)*za(jd,jg)*za(ju,j
     &     g)*zb(jc,jb))
       end function streal_lightGL_PPMM_P_SM

       function streal_lightGL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_PPMM_M_SM =
     &     (2*gb**2*propW34*zb(je,ju)*(za(ju,jn)*zb(jc,ju)*(za(jb,jd)*zb(jd,ju) + z
     &     a(jb,jg)*zb(jg,ju)) - za(jn,je)*(za(jb,jd)*(zb(jc,je)*zb(jd,ju) + ((za(j
     &     b,jg)*zb(jc,jb) - za(jd,jg)*zb(jd,jc))*zb(je,ju)*zb(jg,jd))/(s(jb,jc) + 
     &     s(jb,jd) + s(jc,jd))) + za(jb,jg)*zb(jc,je)*zb(jg,ju))))/(9._dp*ecossin**2*
     &     (s(je,jn) + s(je,ju) + s(jn,ju))*za(jb,jc)*zb(jc,jb)*zb(jg,jd)*zb(jg,ju)
     &     )
       end function streal_lightGL_PPMM_M_SM

       function streal_heavyGR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_MMMM_P_SM =
     &     (4*gb**2*propW34*za(jn,jd)*(za(jb,jc)*zb(je,jb) + za(ju,jc)*zb(je,ju) + 
     &     za(jc,jg)*zb(jg,je))*(-(za(jb,jc)*zb(jb,ju)) + za(jc,jg)*zb(jg,ju)))/(9.
     &     *ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je) + s(jd,jn) + s(je
     &     ,jn))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGR_MMMM_P_SM

       function streal_heavyGR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_MMMM_M_SM =
     &     (-4*gb**2*propW34*za(jn,jd)*zb(jb,ju)*(za(jb,jc)*zb(jc,jb)*zb(je,jb) + z
     &     a(ju,jc)*zb(jc,jb)*zb(je,ju) - za(jc,jg)*zb(jc,je)*zb(jg,jb) + za(jb,jg)
     &     *zb(je,jb)*zb(jg,jb) + za(ju,jg)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jc,j
     &     b)*zb(jg,je)))/(9._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je)
     &      + s(jd,jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyGR_MMMM_M_SM

       function streal_heavyGR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_PPMM_P_SM =
     &     (-4*gb**2*propW34*za(jn,jd)*(za(jb,jc)*zb(jc,je) + za(ju,jb)*zb(je,ju) +
     &      za(jb,jg)*zb(jg,je))*(za(jb,jc)*zb(jc,ju) + za(jb,jg)*zb(jg,ju)))/(9._dp*e
     &     cossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je) + s(jd,jn) + s(je,j
     &     n))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGR_PPMM_P_SM

       function streal_heavyGR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_PPMM_M_SM =
     &     (-4*gb**2*propW34*za(jn,jd)*zb(jc,ju)*(za(jb,jc)*zb(jc,jb)*zb(jc,je) + z
     &     a(ju,jb)*zb(jc,jb)*zb(je,ju) + za(jc,jg)*zb(jc,je)*zb(jg,jc) - za(jb,jg)
     &     *zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*zb(jc,j
     &     b)*zb(jg,je)))/(9._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je)
     &      + s(jd,jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyGR_PPMM_M_SM

       function streal_heavyGL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_MMMM_P_SM =
     &     (2*gb**2*propW34*za(jc,jd)*zb(je,ju)*(za(jb,jc)*(za(ju,jn)*zb(jb,ju) + z
     &     a(jn,je)*zb(je,jb)) + za(jc,jg)*(za(jn,je)*zb(jg,je) - za(ju,jn)*zb(jg,j
     &     u))))/(9._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(je,jn) + s(je,j
     &     u) + s(jn,ju))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGL_MMMM_P_SM

       function streal_heavyGL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_MMMM_M_SM =
     &     (2*gb**2*propW34*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb))*zb(je,ju)*(
     &     za(ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb)))/(9._dp*ec
     &     ossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(je,jn) + s(je,ju) + s(jn,ju
     &     ))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyGL_MMMM_M_SM

       function streal_heavyGL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_PPMM_P_SM =
     &     (-2*gb**2*propW34*za(jb,jd)*zb(je,ju)*(-(za(jn,je)*(za(jb,jd)*zb(jd,je) 
     &     - za(jb,jn)*zb(je,jn) + za(ju,jb)*zb(je,ju))) + za(ju,jn)*(za(jb,jd)*zb(
     &     jd,ju) + za(jb,je)*zb(je,ju) + za(jb,jn)*zb(jn,ju))))/(9._dp*ecossin**2*(s(
     &     jb,jc) + s(jb,jg) + s(jc,jg))*(s(je,jn) + s(je,ju) + s(jn,ju))*za(jb,jg)
     &     *za(jc,jg))
       end function streal_heavyGL_PPMM_P_SM

       function streal_heavyGL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_PPMM_M_SM =
     &     (2*gb**2*propW34*(-(za(jn,je)*zb(jc,je)) + za(ju,jn)*zb(jc,ju))*zb(je,ju
     &     )*(za(jb,jd)*zb(jc,jb) + za(jd,jg)*zb(jg,jc)))/(9._dp*ecossin**2*(s(jb,jc) 
     &     + s(jb,jg) + s(jc,jg))*(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jb)*zb(jg,
     &     jc))
       end function streal_heavyGL_PPMM_M_SM

       function streal_lightZR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_MMMM_P_SM =
     &     -((gb**4 - 9*gw**4)*propW34*propZ25*za(jn,jd)*(za(jc,jd)*zb(jd,je)*(za(j
     &     u,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)) + zb(je,jn)*(za(jn,jc)*(za(ju,jd)
     &     *zb(jb,ju) + za(jd,jg)*zb(jg,jb)) + (za(jn,jd)*za(ju,jg)*zb(jb,ju)*(za(j
     &     b,jc)*zb(jg,jb) + za(ju,jc)*zb(jg,ju)))/(s(jb,jc) + s(jb,ju) + s(jc,ju))
     &     )))/(18._dp*ecossin**2*gw**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(jd,jg)*za(
     &     ju,jg))
       end function streal_lightZR_MMMM_P_SM

       function streal_lightZR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_MMMM_M_SM =
     &     ((gb**4 - 9*gw**4)*propW34*propZ25*zb(jb,ju)*(za(ju,jc)*zb(je,ju)*(za(jn
     &     ,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) + za(jb,jc)*(za(jn,jd)*(zb(jd,ju)*
     &     zb(je,jb) + (zb(jb,ju)*(-(za(jd,jg)*zb(jd,je)) + za(jn,jg)*zb(je,jn))*zb
     &     (jg,jd))/(s(jd,je) + s(jd,jn) + s(je,jn))) + za(jn,jg)*zb(je,jb)*zb(jg,j
     &     u))))/(18._dp*ecossin**2*gw**2*(s(jb,jc) + s(jb,ju) + s(jc,ju))*zb(jg,jd)*z
     &     b(jg,ju))
       end function streal_lightZR_MMMM_M_SM

       function streal_lightZR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_PPMM_P_SM =
     &     (gb**2*(gb**2 - 3*gw**2)*propW34*propZ25*za(jn,jd)*(za(jb,jd)*zb(jd,je)*
     &     (za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc)) + zb(je,jn)*(-(za(jb,jn)*(za
     &     (ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc))) + (za(jn,jd)*za(ju,jg)*zb(jc,j
     &     u)*(-(za(jb,jc)*zb(jg,jc)) + za(ju,jb)*zb(jg,ju)))/(s(jb,jc) + s(jb,ju) 
     &     + s(jc,ju)))))/(9._dp*ecossin**2*gw**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(
     &     jd,jg)*za(ju,jg))
       end function streal_lightZR_PPMM_P_SM

       function streal_lightZR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_PPMM_M_SM =
     &     -(gb**2*(gb**2 - 3*gw**2)*propW34*propZ25*zb(jc,ju)*(za(ju,jb)*zb(je,ju)
     &     *(za(jn,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) + za(jb,jc)*(za(jn,jd)*(zb(
     &     jc,je)*zb(jd,ju) + (zb(jc,ju)*(za(jd,jg)*zb(jd,je) - za(jn,jg)*zb(je,jn)
     &     )*zb(jg,jd))/(s(jd,je) + s(jd,jn) + s(je,jn))) + za(jn,jg)*zb(jc,je)*zb(
     &     jg,ju))))/(9._dp*ecossin**2*gw**2*(s(jb,jc) + s(jb,ju) + s(jc,ju))*zb(jg,jd
     &     )*zb(jg,ju))
       end function streal_lightZR_PPMM_M_SM

       function streal_lightZL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_MMMM_P_SM =
     &     -((gb**2 + 3*gw**2)**2*propW34*propZ25*za(jc,jd)*(za(jn,jc)*zb(jc,jb)*(z
     &     a(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + za(jn,jd)*zb(jd,jb)*(za(ju,j
     &     d)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + (za(jc,jd)*za(ju,jg)*zb(jc,jb)*zb(
     &     je,ju)*(za(jn,je)*zb(jg,je) - za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(je,ju)
     &      + s(jn,ju))))/(18._dp*ecossin**2*gw**2*(s(jb,jc) + s(jb,jd) + s(jc,jd))*za
     &     (jd,jg)*za(ju,jg))
       end function streal_lightZL_MMMM_P_SM

       function streal_lightZL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_MMMM_M_SM =
     &     ((gb**2 + 3*gw**2)**2*propW34*propZ25*zb(je,ju)*(za(ju,jn)*zb(jb,ju)*(za
     &     (jc,jd)*zb(jd,ju) + za(jc,jg)*zb(jg,ju)) + za(jn,je)*(za(jc,jd)*(zb(jd,j
     &     u)*zb(je,jb) + ((za(jc,jg)*zb(jc,jb) + za(jd,jg)*zb(jd,jb))*zb(je,ju)*zb
     &     (jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) + za(jc,jg)*zb(je,jb)*zb(jg,j
     &     u))))/(18._dp*ecossin**2*gw**2*(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jd)*z
     &     b(jg,ju))
       end function streal_lightZL_MMMM_M_SM

       function streal_lightZL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_PPMM_P_SM =
     &     (gb**2*(gb**2 + 3*gw**2)*propW34*propZ25*za(jb,jd)*(za(jb,jn)*zb(jc,jb)*
     &     (za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + za(jn,jd)*zb(jd,jc)*(za(ju
     &     ,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + (za(jb,jd)*za(ju,jg)*zb(jc,jb)*z
     &     b(je,ju)*(-(za(jn,je)*zb(jg,je)) + za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(j
     &     e,ju) + s(jn,ju))))/(9._dp*ecossin**2*gw**2*(s(jb,jc) + s(jb,jd) + s(jc,jd)
     &     )*za(jd,jg)*za(ju,jg))
       end function streal_lightZL_PPMM_P_SM

       function streal_lightZL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_PPMM_M_SM =
     &     (gb**2*(gb**2 + 3*gw**2)*propW34*propZ25*zb(je,ju)*(-(za(ju,jn)*zb(jc,ju
     &     )*(za(jb,jd)*zb(jd,ju) + za(jb,jg)*zb(jg,ju))) + za(jn,je)*(za(jb,jd)*(z
     &     b(jc,je)*zb(jd,ju) + ((za(jb,jg)*zb(jc,jb) - za(jd,jg)*zb(jd,jc))*zb(je,
     &     ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) + za(jb,jg)*zb(jc,je)*z
     &     b(jg,ju))))/(9._dp*ecossin**2*gw**2*(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,
     &     jd)*zb(jg,ju))
       end function streal_lightZL_PPMM_M_SM

       function streal_heavyZR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_MMMM_P_SM =
     &     ((gb**4 - 9*gw**4)*propW34*propZ257*za(jn,jd)*(za(jb,jc)*zb(je,jb) + za(
     &     ju,jc)*zb(je,ju) + za(jc,jg)*zb(jg,je))*(za(jb,jc)*zb(jb,ju) - za(jc,jg)
     &     *zb(jg,ju)))/(18._dp*ecossin**2*gw**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(j
     &     b,jg)*za(jc,jg))
       end function streal_heavyZR_MMMM_P_SM

       function streal_heavyZR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_MMMM_M_SM =
     &     ((gb**4 - 9*gw**4)*propW34*propZ257*za(jn,jd)*zb(jb,ju)*(za(jb,jc)*zb(jc
     &     ,jb)*zb(je,jb) + za(ju,jc)*zb(jc,jb)*zb(je,ju) - za(jc,jg)*zb(jc,je)*zb(
     &     jg,jb) + za(jb,jg)*zb(je,jb)*zb(jg,jb) + za(ju,jg)*zb(je,ju)*zb(jg,jb) +
     &      za(jc,jg)*zb(jc,jb)*zb(jg,je)))/(18._dp*ecossin**2*gw**2*(s(jd,je) + s(jd,
     &     jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZR_MMMM_M_SM

       function streal_heavyZR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_PPMM_P_SM =
     &     -(gb**2*(gb**2 - 3*gw**2)*propW34*propZ257*za(jn,jd)*(za(jb,jc)*zb(jc,je
     &     ) + za(ju,jb)*zb(je,ju) + za(jb,jg)*zb(jg,je))*(za(jb,jc)*zb(jc,ju) + za
     &     (jb,jg)*zb(jg,ju)))/(9._dp*ecossin**2*gw**2*(s(jd,je) + s(jd,jn) + s(je,jn)
     &     )*za(jb,jg)*za(jc,jg))
       end function streal_heavyZR_PPMM_P_SM

       function streal_heavyZR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_PPMM_M_SM =
     &     -(gb**2*(gb**2 - 3*gw**2)*propW34*propZ257*za(jn,jd)*zb(jc,ju)*(za(jb,jc
     &     )*zb(jc,jb)*zb(jc,je) + za(ju,jb)*zb(jc,jb)*zb(je,ju) + za(jc,jg)*zb(jc,
     &     je)*zb(jg,jc) - za(jb,jg)*zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(je,ju)*zb(j
     &     g,jc) + za(jb,jg)*zb(jc,jb)*zb(jg,je)))/(9._dp*ecossin**2*gw**2*(s(jd,je) +
     &      s(jd,jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZR_PPMM_M_SM

       function streal_heavyZL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_MMMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_MMMM_P_SM =
     &     ((gb**2 + 3*gw**2)**2*propW34*propZ257*za(jc,jd)*zb(je,ju)*(za(jb,jc)*(z
     &     a(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb)) + za(jc,jg)*(za(jn,je)*zb(jg,j
     &     e) - za(ju,jn)*zb(jg,ju))))/(18._dp*ecossin**2*gw**2*(s(je,jn) + s(je,ju) +
     &      s(jn,ju))*za(jb,jg)*za(jc,jg))
       end function streal_heavyZL_MMMM_P_SM

       function streal_heavyZL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_MMMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_MMMM_M_SM =
     &     ((gb**2 + 3*gw**2)**2*propW34*propZ257*(za(ju,jn)*zb(jb,ju) + za(jn,je)*
     &     zb(je,jb))*zb(je,ju)*(za(ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,
     &     jd)*zb(jn,jb)))/(18._dp*ecossin**2*gw**2*(s(je,jn) + s(je,ju) + s(jn,ju))*z
     &     b(jg,jb)*zb(jg,jc))
       end function streal_heavyZL_MMMM_M_SM

       function streal_heavyZL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_PPMM_P_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_PPMM_P_SM =
     &     (gb**2*(gb**2 + 3*gw**2)*propW34*propZ257*za(jb,jd)*zb(je,ju)*(-(za(jn,j
     &     e)*(za(jb,jd)*zb(jd,je) - za(jb,jn)*zb(je,jn) + za(ju,jb)*zb(je,ju))) + 
     &     za(ju,jn)*(za(jb,jd)*zb(jd,ju) + za(jb,je)*zb(je,ju) + za(jb,jn)*zb(jn,j
     &     u))))/(9._dp*ecossin**2*gw**2*(s(je,jn) + s(je,ju) + s(jn,ju))*za(jb,jg)*za
     &     (jc,jg))
       end function streal_heavyZL_PPMM_P_SM

       function streal_heavyZL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_PPMM_M_SM
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_PPMM_M_SM =
     &     (gb**2*(gb**2 + 3*gw**2)*propW34*propZ257*(za(jn,je)*zb(jc,je) - za(ju,j
     &     n)*zb(jc,ju))*zb(je,ju)*(za(jb,jd)*zb(jc,jb) + za(jd,jg)*zb(jg,jc)))/(9.
     &     *ecossin**2*gw**2*(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZL_PPMM_M_SM

       function streal_lightResonant_MMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_P_L2 =
     &     (-2*propT1267*propW34*za(jn,jc)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,j
     &     b))*(propW167*real(anomc1)*(za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) + 
     &     za(jd,jg)*zb(jg,je)) + 2*real(anomc8)*(za(jb,jd)*zb(je,jb) + za(ju,jd)*z
     &     b(je,ju) + za(jd,jg)*zb(jg,je)) - 2*Sqrt(mtsq)*propW167*(real(anomc3)*(z
     &     a(jn,jd)*zb(je,jn) - za(ju,jd)*zb(je,ju) - za(jd,jg)*zb(jg,je)) + im*ima
     &     g(anomc3)*(za(jn,jd)*zb(je,jn) + za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,j
     &     e)))))/(za(jd,jg)*za(ju,jg))
       end function streal_lightResonant_MMMM_P_L2

       function streal_lightResonant_MMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_M_L2 =
     &     (-2*propT1267*propW34*za(jn,jc)*(propW167*real(anomc1)*zb(jb,ju)*(za(jb,
     &     jd)*zb(jd,ju)*zb(je,jb) + za(ju,jd)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*zb(j
     &     d,ju)*zb(jg,je) - za(jd,jg)*zb(jd,je)*zb(jg,ju) + za(jb,jg)*zb(je,jb)*zb
     &     (jg,ju) + za(ju,jg)*zb(je,ju)*zb(jg,ju)) + 2*real(anomc8)*zb(jb,ju)*(za(
     &     jb,jd)*zb(jd,ju)*zb(je,jb) + za(ju,jd)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*z
     &     b(jd,ju)*zb(jg,je) - za(jd,jg)*zb(jd,je)*zb(jg,ju) + za(jb,jg)*zb(je,jb)
     &     *zb(jg,ju) + za(ju,jg)*zb(je,ju)*zb(jg,ju)) - Sqrt(mtsq)*propW167*(im*im
     &     ag(anomc3)*(2*za(jn,jd)*zb(jb,ju)*zb(jd,ju)*zb(je,jn) + 2*za(ju,jd)*zb(j
     &     b,ju)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) + za
     &     (jd,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) - za(jd,jg)*zb(jb,ju)*zb(jd,je)*zb
     &     (jg,ju) + 2*za(jn,jg)*zb(jb,ju)*zb(je,jn)*zb(jg,ju) + 2*za(ju,jg)*zb(jb,
     &     ju)*zb(je,ju)*zb(jg,ju) - za(jd,jg)*zb(jd,jb)*zb(je,ju)*zb(jg,ju)) + rea
     &     l(anomc3)*(2*za(jn,jd)*zb(jb,ju)*zb(jd,ju)*zb(je,jn) - 2*za(ju,jd)*zb(jb
     &     ,ju)*zb(jd,ju)*zb(je,ju) - za(jd,jg)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) - za(
     &     jd,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za(jd,jg)*zb(jb,ju)*zb(jd,je)*zb(
     &     jg,ju) + 2*za(jn,jg)*zb(jb,ju)*zb(je,jn)*zb(jg,ju) - 2*za(ju,jg)*zb(jb,j
     &     u)*zb(je,ju)*zb(jg,ju) + za(jd,jg)*zb(jd,jb)*zb(je,ju)*zb(jg,ju)))))/(zb
     &     (jg,jd)*zb(jg,ju))
       end function streal_lightResonant_MMMM_M_L2

       function streal_lightResonant_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MPMM_P_L2 =
     &     (propW167*propW34*(Sqrt(mtsq)*propT1267*(im*imag(anomc2) + real(anomc2))
     &     *za(jb,jd)*za(jn,jc)*(za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)) + 2*im*
     &     imag(anomc4)*(2*propT1267*za(jb,jd)**2*za(jn,jc)*zb(je,jb)*(za(ju,jd)*zb
     &     (jd,ju) + za(jd,jg)*zb(jg,jd)) + (za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,
     &     je))*(za(jb,jn)*za(jc,jd) + propT1267*za(jn,jc)*(-(za(jd,jg)*za(ju,jb)) 
     &     + za(jb,jg)*za(ju,jd))*zb(jg,ju)) + za(jb,jd)*za(jn,jc)*(2*propT1267*za(
     &     ju,jd)**2*zb(jd,ju)*zb(je,ju) + za(jd,jg)*(-2*propT1267*za(ju,jb)*zb(je,
     &     jb)*zb(jg,ju) + zb(jg,je)*(-1 + 2*propT1267*za(jd,jg)*zb(jg,jd) + propT1
     &     267*za(ju,jg)*zb(jg,ju))) + za(ju,jd)*(zb(je,ju)*(-1 + 2*propT1267*za(jd
     &     ,jg)*zb(jg,jd) + propT1267*za(ju,jg)*zb(jg,ju)) + 2*propT1267*(za(jd,jg)
     &     *zb(jd,ju)*zb(jg,je) + za(jb,jg)*zb(je,jb)*zb(jg,ju))))) + 2*real(anomc4
     &     )*(2*propT1267*za(jb,jd)**2*za(jn,jc)*zb(je,jb)*(za(ju,jd)*zb(jd,ju) + z
     &     a(jd,jg)*zb(jg,jd)) + (za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je))*(za(jb
     &     ,jn)*za(jc,jd) + propT1267*za(jn,jc)*(-(za(jd,jg)*za(ju,jb)) + za(jb,jg)
     &     *za(ju,jd))*zb(jg,ju)) + za(jb,jd)*za(jn,jc)*(2*propT1267*za(ju,jd)**2*z
     &     b(jd,ju)*zb(je,ju) + za(jd,jg)*(-2*propT1267*za(ju,jb)*zb(je,jb)*zb(jg,j
     &     u) + zb(jg,je)*(-1 + 2*propT1267*za(jd,jg)*zb(jg,jd) + propT1267*za(ju,j
     &     g)*zb(jg,ju))) + za(ju,jd)*(zb(je,ju)*(-1 + 2*propT1267*za(jd,jg)*zb(jg,
     &     jd) + propT1267*za(ju,jg)*zb(jg,ju)) + 2*propT1267*(za(jd,jg)*zb(jd,ju)*
     &     zb(jg,je) + za(jb,jg)*zb(je,jb)*zb(jg,ju)))))))/(za(jd,jg)*za(ju,jg))
       end function streal_lightResonant_MPMM_P_L2

       function streal_lightResonant_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MPMM_M_L2 =
     &     (propW167*(im*Sqrt(mtsq)*propT1267*propW34*imag(anomc2)*za(jn,jc)*zb(je,
     &     ju)*(za(jb,jd)*zb(jd,ju) + za(jb,jg)*zb(jg,ju)) + Sqrt(mtsq)*propT1267*p
     &     ropW34*real(anomc2)*za(jn,jc)*zb(je,ju)*(za(jb,jd)*zb(jd,ju) + za(jb,jg)
     &     *zb(jg,ju)) + 2*im*imag(anomc4)*(2*propT1267*propW34*za(jn,jc)*(za(jb,jd
     &     )*zb(jd,ju) + za(jb,jg)*zb(jg,ju))*(za(jb,jd)*zb(jd,ju)*zb(je,jb) + za(j
     &     u,jd)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*zb(jd,ju)*zb(jg,je) - za(jd,jg)*zb
     &     (jd,je)*zb(jg,ju) + za(jb,jg)*zb(je,jb)*zb(jg,ju) + za(ju,jg)*zb(je,ju)*
     &     zb(jg,ju)) + propW34*zb(je,ju)*(-(za(jn,jc)*(za(jb,jd)*zb(jd,ju) + za(jb
     &     ,jg)*zb(jg,ju))) + za(jb,jn)*(za(jc,jd)*zb(jd,ju) + za(jc,jg)*zb(jg,ju))
     &     )) + 2*real(anomc4)*(2*propT1267*propW34*za(jn,jc)*(za(jb,jd)*zb(jd,ju) 
     &     + za(jb,jg)*zb(jg,ju))*(za(jb,jd)*zb(jd,ju)*zb(je,jb) + za(ju,jd)*zb(jd,
     &     ju)*zb(je,ju) + za(jd,jg)*zb(jd,ju)*zb(jg,je) - za(jd,jg)*zb(jd,je)*zb(j
     &     g,ju) + za(jb,jg)*zb(je,jb)*zb(jg,ju) + za(ju,jg)*zb(je,ju)*zb(jg,ju)) +
     &      propW34*zb(je,ju)*(-(za(jn,jc)*(za(jb,jd)*zb(jd,ju) + za(jb,jg)*zb(jg,j
     &     u))) + za(jb,jn)*(za(jc,jd)*zb(jd,ju) + za(jc,jg)*zb(jg,ju))))))/(zb(jg,
     &     jd)*zb(jg,ju))
       end function streal_lightResonant_MPMM_M_L2

       function streal_lightResonant_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_PMMM_P_L2 =
     &     (propW167*propW34*(im*Sqrt(mtsq)*propT1267*imag(anomc2)*za(jn,jd)*zb(jc,
     &     je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)) - Sqrt(mtsq)*propT1267*r
     &     eal(anomc2)*za(jn,jd)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,j
     &     b)) - 2*im*imag(anomc4)*(za(jn,jd)*(za(ju,jd)*(zb(jb,ju)*zb(jc,je) - zb(
     &     jc,ju)*zb(je,jb)) + za(jd,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc)
     &     )) + 2*propT1267*za(jn,je)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb
     &     (jg,jb))*(za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je
     &     ))) + 2*real(anomc4)*(za(jn,jd)*(za(ju,jd)*(zb(jb,ju)*zb(jc,je) - zb(jc,
     &     ju)*zb(je,jb)) + za(jd,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc))) 
     &     + 2*propT1267*za(jn,je)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg
     &     ,jb))*(za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) + za(jd,jg)*zb(jg,je)))
     &     ))/(za(jd,jg)*za(ju,jg))
       end function streal_lightResonant_PMMM_P_L2

       function streal_lightResonant_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_PMMM_M_L2 =
     &     (propW167*(im*Sqrt(mtsq)*propT1267*propW34*imag(anomc2)*zb(jb,ju)*zb(jc,
     &     je)*(za(jn,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) - Sqrt(mtsq)*propT1267*p
     &     ropW34*real(anomc2)*zb(jb,ju)*zb(jc,je)*(za(jn,jd)*zb(jd,ju) + za(jn,jg)
     &     *zb(jg,ju)) + 2*im*imag(anomc4)*(propW34*(-(zb(jb,ju)*zb(jc,je)) + zb(jc
     &     ,ju)*zb(je,jb))*(za(jn,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) - 2*propT126
     &     7*propW34*za(jn,je)*zb(jb,ju)*zb(jc,je)*(za(jb,jd)*zb(jd,ju)*zb(je,jb) +
     &      za(ju,jd)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*zb(jd,ju)*zb(jg,je) - za(jd,j
     &     g)*zb(jd,je)*zb(jg,ju) + za(jb,jg)*zb(je,jb)*zb(jg,ju) + za(ju,jg)*zb(je
     &     ,ju)*zb(jg,ju))) + 2*real(anomc4)*(propW34*(zb(jb,ju)*zb(jc,je) - zb(jc,
     &     ju)*zb(je,jb))*(za(jn,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) + 2*propT1267
     &     *propW34*za(jn,je)*zb(jb,ju)*zb(jc,je)*(za(jb,jd)*zb(jd,ju)*zb(je,jb) + 
     &     za(ju,jd)*zb(jd,ju)*zb(je,ju) + za(jd,jg)*zb(jd,ju)*zb(jg,je) - za(jd,jg
     &     )*zb(jd,je)*zb(jg,ju) + za(jb,jg)*zb(je,jb)*zb(jg,ju) + za(ju,jg)*zb(je,
     &     ju)*zb(jg,ju)))))/(zb(jg,jd)*zb(jg,ju))
       end function streal_lightResonant_PMMM_M_L2

       function streal_lightResonant_MPPP_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MPPP_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MPPP_P_L2 =
     &     (2*Sqrt(mtsq)*propT1267*propW34*real(anomc9)*za(jn,jc)*za(ju,jb)*(za(ju,
     &     jd)*zb(jd,je) + za(ju,jg)*zb(jg,je)))/(za(jd,jg)*za(ju,jg))
       end function streal_lightResonant_MPPP_P_L2

       function streal_lightResonant_MPPP_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MPPP_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MPPP_M_L2 =
     &     (2*Sqrt(mtsq)*propT1267*propW34*real(anomc9)*za(jn,jc)*zb(jd,je)*(za(ju,
     &     jb)*zb(jd,ju) + za(jb,jg)*zb(jg,jd)))/(zb(jg,jd)*zb(jg,ju))
       end function streal_lightResonant_MPPP_M_L2

       function streal_heavyResonant_MMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MMMM_P_L2 =
     &     (2*propW34*za(jn,jc)*(im*Sqrt(mtsq)*propT126*propT1267*propW16*imag(anom
     &     c6)*za(jb,jg)*za(jc,jg)*za(jd,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je)*(za(jb,j
     &     d)*zb(jg,jb) + za(ju,jd)*zb(jg,ju)) - Sqrt(mtsq)*propT126*propT1267*prop
     &     W16*real(anomc6)*za(jb,jg)*za(jc,jg)*za(jd,jg)*za(jn,jg)*zb(jb,ju)*zb(jg
     &     ,je)*(za(jb,jd)*zb(jg,jb) + za(ju,jd)*zb(jg,ju)) + propW16*real(anomc1)*
     &     za(jd,jg)*(propT1267*za(jc,jg)*(-(za(jc,jd)*zb(jc,je)) + za(jn,jd)*zb(je
     &     ,jn))*(za(jb,jn)*zb(jb,ju) - za(jn,jg)*zb(jg,ju)) + propT126*za(jb,jg)*z
     &     b(jb,ju)*(za(jb,jd)*za(jn,jc)*(-zb(je,jb) + propT1267*za(jc,jg)*zb(jc,je
     &     )*zb(jg,jb)) + mtsq*propT1267*za(jc,jg)*za(jn,jd)*zb(jg,je) + za(jn,jc)*
     &     za(ju,jd)*(-zb(je,ju) + propT1267*za(jc,jg)*zb(jc,je)*zb(jg,ju)))) + 2*r
     &     eal(anomc8)*za(jd,jg)*(propT1267*za(jc,jg)*(-(za(jc,jd)*zb(jc,je)) + za(
     &     jn,jd)*zb(je,jn))*(za(jb,jn)*zb(jb,ju) - za(jn,jg)*zb(jg,ju)) + propT126
     &     *za(jb,jg)*zb(jb,ju)*(za(jb,jd)*za(jn,jc)*(-zb(je,jb) + propT1267*za(jc,
     &     jg)*zb(jc,je)*zb(jg,jb)) + mtsq*propT1267*za(jc,jg)*za(jn,jd)*zb(jg,je) 
     &     + za(jn,jc)*za(ju,jd)*(-zb(je,ju) + propT1267*za(jc,jg)*zb(jc,je)*zb(jg,
     &     ju)))) + 2*im*Sqrt(mtsq)*propW16*imag(anomc3)*za(jn,jg)*(propT1267*za(jc
     &     ,jg)*(za(jn,jd)*zb(je,jn) + za(ju,jd)*zb(je,ju))*(za(jb,jd)*zb(jb,ju) - 
     &     za(jd,jg)*zb(jg,ju)) + propT126*za(jb,jg)*zb(jb,ju)*(propT1267*za(jb,jd)
     &     *za(jc,jg)*(za(jn,jd)*zb(je,jn)*zb(jg,jb) + za(ju,jd)*zb(jb,ju)*zb(jg,je
     &     )) - za(jc,jd)*(za(jn,jd)*zb(je,jn) + za(ju,jd)*(zb(je,ju) - propT1267*z
     &     a(jc,jg)*zb(jc,je)*zb(jg,ju))))) + 2*Sqrt(mtsq)*propW16*real(anomc3)*za(
     &     jn,jg)*(propT1267*za(jc,jg)*(za(jn,jd)*zb(je,jn) - za(ju,jd)*zb(je,ju))*
     &     (za(jb,jd)*zb(jb,ju) - za(jd,jg)*zb(jg,ju)) + propT126*za(jb,jg)*zb(jb,j
     &     u)*(za(jn,jd)*zb(je,jn)*(-za(jc,jd) + propT1267*za(jb,jd)*za(jc,jg)*zb(j
     &     g,jb)) + za(ju,jd)*(za(jc,jd)*(zb(je,ju) - propT1267*za(jc,jg)*zb(jc,je)
     &     *zb(jg,ju)) + propT1267*za(jc,jg)*(-(za(jb,jd)*zb(jb,ju)*zb(jg,je)) + 2*
     &     za(jn,jd)*zb(je,jn)*zb(jg,ju)))))))/(za(jb,jg)*za(jc,jg)*za(jd,jg)*za(jn
     &     ,jg))
       end function streal_heavyResonant_MMMM_P_L2

       function streal_heavyResonant_MMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MMMM_M_L2 =
     &     (-2*propW34*zb(jb,ju)*(-(Sqrt(mtsq)*propW16*(2*im*imag(anomc3)*(propT126
     &     *za(jn,jg)*(za(jn,jd)*zb(je,jn) + za(ju,jd)*zb(je,ju))*zb(jg,jb) + propT
     &     1267*za(jn,jc)*(za(jn,jd)*zb(jc,jb)*zb(je,jn) + propT126*(za(jb,jd)*za(j
     &     n,jg)*zb(jc,jb) - za(jd,jg)*za(jn,je)*zb(jc,je))*zb(je,jn)*zb(jg,jb) + z
     &     a(ju,jd)*(zb(jc,jb)*zb(je,ju) + propT126*zb(jc,je)*(za(jb,jg)*zb(jb,ju) 
     &     + za(jc,jg)*zb(jc,ju) + za(jd,jg)*zb(jd,ju))*zb(jg,jb)))) + 2*real(anomc
     &     3)*(propT126*za(jn,jg)*(za(jn,jd)*zb(je,jn) - za(ju,jd)*zb(je,ju))*zb(jg
     &     ,jb) + propT1267*za(jn,jc)*(za(jn,jd)*zb(jc,jb)*zb(je,jn) + propT126*(za
     &     (jb,jd)*za(jn,jg)*zb(jc,jb) - za(jd,jg)*za(jn,je)*zb(jc,je))*zb(je,jn)*z
     &     b(jg,jb) - za(ju,jd)*(zb(jc,jb)*zb(je,ju) + propT126*(za(jb,jg)*zb(jb,ju
     &     )*zb(jc,je) + za(jc,jg)*zb(jc,je)*zb(jc,ju) + za(jd,jg)*zb(jc,je)*zb(jd,
     &     ju) - 2*za(jn,jg)*zb(jc,ju)*zb(je,jn))*zb(jg,jb)))) + propT126*propT1267
     &     *(im*imag(anomc6) + real(anomc6))*za(jd,jg)*za(jn,jc)*(za(jc,jg)*zb(jc,j
     &     e) - za(jn,jg)*zb(je,jn))*zb(jg,jb)*zb(jg,jc))) + propT126*propW16*real(
     &     anomc1)*(za(jn,jg)*za(ju,jd)*zb(je,ju)*zb(jg,jb) + za(jb,jd)*zb(je,jb)*(
     &     za(jn,jc)*zb(jc,jb) + za(jn,jg)*zb(jg,jb)) + za(jn,jc)*(-(mtsq*propT1267
     &     *za(jd,jg)*zb(je,jb)*zb(jg,jc)) + za(ju,jd)*(zb(jc,jb)*zb(je,ju) + propT
     &     1267*zb(jb,ju)*(za(jc,jg)*zb(jc,je) - za(jn,jg)*zb(je,jn))*zb(jg,jc)))) 
     &     + 2*propT126*real(anomc8)*(za(jn,jg)*za(ju,jd)*zb(je,ju)*zb(jg,jb) + za(
     &     jb,jd)*zb(je,jb)*(za(jn,jc)*zb(jc,jb) + za(jn,jg)*zb(jg,jb)) + za(jn,jc)
     &     *(-(mtsq*propT1267*za(jd,jg)*zb(je,jb)*zb(jg,jc)) + za(ju,jd)*(zb(jc,jb)
     &     *zb(je,ju) + propT1267*zb(jb,ju)*(za(jc,jg)*zb(jc,je) - za(jn,jg)*zb(je,
     &     jn))*zb(jg,jc))))))/(zb(jg,jb)*zb(jg,jc))
       end function streal_heavyResonant_MMMM_M_L2

       function streal_heavyResonant_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MPMM_P_L2 =
     &     -((propW16*propW34*(Sqrt(mtsq)*propT126*(im*imag(anomc2) + real(anomc2))
     &     *za(jb,jd)*za(jn,jc)*(za(jb,jc)*(-zb(je,ju) + propT1267*za(jc,jg)*zb(jc,
     &     je)*zb(jg,ju)) + propT1267*za(jc,jg)*(za(jb,jd)*zb(jd,ju)*zb(jg,je) - za
     &     (jb,jn)*zb(je,jn)*zb(jg,ju))) + 2*im*imag(anomc4)*(za(jb,jc)*(-(za(jb,jn
     &     )*za(jc,jd)*zb(je,ju)) + za(jb,jd)*za(jn,jc)*(-2*propT1267*za(jc,jd)*zb(
     &     jc,je)*zb(jd,ju) + 2*propT1267*za(jn,jd)*zb(jd,ju)*zb(je,jn) + zb(je,ju)
     &     )) + 2*propT126*propT1267*za(jb,jd)*za(jb,jg)*za(jn,jc)*zb(jd,ju)*(za(jb
     &     ,jd)*za(jn,jc)*zb(je,jn)*zb(jg,jb) + mtsq*za(jc,jd)*zb(jg,je) + za(jn,jc
     &     )*za(ju,jd)*zb(je,jn)*zb(jg,ju))) + real(anomc4)*(za(jb,jc)*(-2*za(jb,jn
     &     )*za(jc,jd)*zb(je,ju) + 2*za(jb,jd)*za(jn,jc)*(-2*propT1267*za(jc,jd)*zb
     &     (jc,je)*zb(jd,ju) + 2*propT1267*za(jn,jd)*zb(jd,ju)*zb(je,jn) + zb(je,ju
     &     ))) + 4*propT126*propT1267*za(jb,jd)*za(jb,jg)*za(jn,jc)*zb(jd,ju)*(za(j
     &     b,jd)*za(jn,jc)*zb(je,jn)*zb(jg,jb) + mtsq*za(jc,jd)*zb(jg,je) + za(jn,j
     &     c)*za(ju,jd)*zb(je,jn)*zb(jg,ju)))))/(za(jb,jg)*za(jc,jg)))
       end function streal_heavyResonant_MPMM_P_L2

       function streal_heavyResonant_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MPMM_M_L2 =
     &     (propW16*propW34*(im*Sqrt(mtsq)*imag(anomc2)*zb(je,ju)*(propT1267*za(jd,
     &     jg)*za(jn,jc)*zb(jg,jc)*zb(jg,ju) + za(jb,jd)*(za(jn,jc)*(propT126*zb(jc
     &     ,ju)*zb(jg,jb) - propT1267*(zb(jb,ju) + propT126*za(jb,jg)*zb(jb,ju)*zb(
     &     jg,jb) + propT126*za(jd,jg)*zb(jd,ju)*zb(jg,jb))*zb(jg,jc)) + propT126*z
     &     a(jn,jg)*zb(jg,jb)*zb(jg,ju))) + Sqrt(mtsq)*real(anomc2)*zb(je,ju)*(prop
     &     T1267*za(jd,jg)*za(jn,jc)*zb(jg,jc)*zb(jg,ju) + za(jb,jd)*(za(jn,jc)*(pr
     &     opT126*zb(jc,ju)*zb(jg,jb) - propT1267*(zb(jb,ju) + propT126*za(jb,jg)*z
     &     b(jb,ju)*zb(jg,jb) + propT126*za(jd,jg)*zb(jd,ju)*zb(jg,jb))*zb(jg,jc)) 
     &     + propT126*za(jn,jg)*zb(jg,jb)*zb(jg,ju))) + 2*zb(jg,ju)*(-(propT1267*(i
     &     m*imag(anomc7) + real(anomc7))*za(jb,jg)*za(jn,jc)*zb(jb,ju)*(za(jc,jd)*
     &     zb(jc,je) - za(jn,jd)*zb(je,jn))*zb(jg,jc)) + im*imag(anomc4)*(2*propT12
     &     6*za(jb,jd)**2*zb(jd,ju)*zb(je,jb)*(za(jn,jc)*zb(jc,jb) + za(jn,jg)*zb(j
     &     g,jb)) + (-(za(jd,jg)*za(jn,jc)*(2*propT1267*za(jn,jd)*zb(jd,ju)*zb(je,j
     &     n) + zb(je,ju))) + za(jc,jd)*(2*propT1267*za(jd,jg)*za(jn,jc)*zb(jc,je)*
     &     zb(jd,ju) + za(jn,jg)*zb(je,ju)))*zb(jg,jc) + za(jb,jd)*(za(ju,jn)*zb(jb
     &     ,ju)*zb(je,ju) + za(jn,jd)*zb(jd,jb)*zb(je,ju) + 2*propT126*za(jn,jc)*za
     &     (ju,jd)*zb(jc,jb)*zb(jd,ju)*zb(je,ju) + za(jn,je)*zb(je,jb)*zb(je,ju) + 
     &     2*propT126*za(jn,jg)*za(ju,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) + 2*propT12
     &     6*propT1267*za(jc,jg)*za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(jc,je)*zb(jd,ju)*
     &     zb(jg,jc) - 2*mtsq*propT126*propT1267*za(jd,jg)*za(jn,jc)*zb(jd,ju)*zb(j
     &     e,jb)*zb(jg,jc) - 2*propT126*propT1267*za(jn,jc)*za(jn,jg)*za(ju,jd)*zb(
     &     jb,ju)*zb(jd,ju)*zb(je,jn)*zb(jg,jc)) + za(jb,jn)*zb(je,ju)*(za(ju,jd)*z
     &     b(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb))) + real(anomc4)*(2
     &     *propT126*za(jb,jd)**2*zb(jd,ju)*zb(je,jb)*(za(jn,jc)*zb(jc,jb) + za(jn,
     &     jg)*zb(jg,jb)) + (-(za(jd,jg)*za(jn,jc)*(2*propT1267*za(jn,jd)*zb(jd,ju)
     &     *zb(je,jn) + zb(je,ju))) + za(jc,jd)*(2*propT1267*za(jd,jg)*za(jn,jc)*zb
     &     (jc,je)*zb(jd,ju) + za(jn,jg)*zb(je,ju)))*zb(jg,jc) + za(jb,jd)*(za(ju,j
     &     n)*zb(jb,ju)*zb(je,ju) + za(jn,jd)*zb(jd,jb)*zb(je,ju) + 2*propT126*za(j
     &     n,jc)*za(ju,jd)*zb(jc,jb)*zb(jd,ju)*zb(je,ju) + za(jn,je)*zb(je,jb)*zb(j
     &     e,ju) + 2*propT126*za(jn,jg)*za(ju,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) + 2
     &     *propT126*propT1267*za(jc,jg)*za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(jc,je)*zb
     &     (jd,ju)*zb(jg,jc) - 2*mtsq*propT126*propT1267*za(jd,jg)*za(jn,jc)*zb(jd,
     &     ju)*zb(je,jb)*zb(jg,jc) - 2*propT126*propT1267*za(jn,jc)*za(jn,jg)*za(ju
     &     ,jd)*zb(jb,ju)*zb(jd,ju)*zb(je,jn)*zb(jg,jc)) + za(jb,jn)*zb(je,ju)*(za(
     &     ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb))))))/(zb(jg
     &     ,jb)*zb(jg,jc)*zb(jg,ju))
       end function streal_heavyResonant_MPMM_M_L2

       function streal_heavyResonant_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_PMMM_P_L2 =
     &     (propW16*propW34*(im*Sqrt(mtsq)*imag(anomc2)*za(jd,jg)*za(jn,jd)*(propT1
     &     267*za(jb,jn)*za(jc,jg)*zb(jb,ju)*zb(jc,je) + propT126*za(jb,jg)*zb(jb,j
     &     u)*(za(jn,jc)*zb(jc,je)*(1 + propT1267*za(jc,jg)*zb(jg,jc)) + (za(jn,jg)
     &      + propT1267*za(jc,jg)*za(jn,je)*zb(jc,je))*zb(jg,je)) - propT1267*za(jc
     &     ,jg)*za(jn,jg)*zb(jc,je)*zb(jg,ju)) - Sqrt(mtsq)*real(anomc2)*za(jd,jg)*
     &     za(jn,jd)*(propT1267*za(jb,jn)*za(jc,jg)*zb(jb,ju)*zb(jc,je) + propT126*
     &     za(jb,jg)*zb(jb,ju)*(za(jn,jc)*zb(jc,je)*(1 + propT1267*za(jc,jg)*zb(jg,
     &     jc)) + (za(jn,jg) + propT1267*za(jc,jg)*za(jn,je)*zb(jc,je))*zb(jg,je)) 
     &     - propT1267*za(jc,jg)*za(jn,jg)*zb(jc,je)*zb(jg,ju)) + 2*za(jn,jg)*(prop
     &     T126*(im*imag(anomc7) - real(anomc7))*za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(j
     &     b,ju)*(za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju))*zb(jg,jc) + im*imag(an
     &     omc4)*(za(jb,jd)*(za(jc,jg)*(-2*propT1267*za(jc,jd)*za(jn,je)*zb(jb,ju)*
     &     zb(jc,je)**2*(1 + propT126*za(jb,jg)*zb(jg,jb)) + za(jn,jd)*(zb(jc,ju)*z
     &     b(je,jb) + zb(jb,ju)*zb(jc,je)*(-1 + 2*propT1267*za(jn,je)*zb(je,jn)*(1 
     &     + propT126*za(jb,jg)*zb(jg,jb))))) + 2*za(jb,jg)*zb(jb,ju)*zb(je,jb)*(za
     &     (jn,jd) + propT126*za(jn,je)*(za(jc,jd)*zb(jc,je) - za(jd,jg)*zb(jg,je))
     &     )) + za(jc,jg)*za(jd,jg)*(2*propT1267*za(jc,jd)*za(jn,je)*zb(jc,je)**2*z
     &     b(jg,ju) + za(jn,jd)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*(1 - 2*propT1267*z
     &     a(jn,je)*zb(je,jn))*zb(jg,ju))) + za(jb,jg)*(-2*propT126*za(jn,je)*za(ju
     &     ,jd)*zb(jb,ju)*(za(jd,jg)*zb(je,ju)*zb(jg,je) + za(jc,jd)*zb(jc,je)*(-zb
     &     (je,ju) + propT1267*za(jc,jg)*zb(jc,je)*zb(jg,ju))) + za(jn,jd)*(za(je,j
     &     d)*zb(je,jb)*zb(je,ju) + za(ju,jd)*zb(jb,ju)*(zb(je,ju) + 2*propT126*pro
     &     pT1267*za(jc,jg)*za(jn,je)*zb(jc,je)*zb(je,jn)*zb(jg,ju))) + za(jn,jd)**
     &     2*(zb(jb,ju)*zb(je,jn) + zb(je,jb)*zb(jn,ju)))) - real(anomc4)*(za(jb,jd
     &     )*(za(jc,jg)*(-2*propT1267*za(jc,jd)*za(jn,je)*zb(jb,ju)*zb(jc,je)**2*(1
     &      + propT126*za(jb,jg)*zb(jg,jb)) + za(jn,jd)*(zb(jc,ju)*zb(je,jb) + zb(j
     &     b,ju)*zb(jc,je)*(-1 + 2*propT1267*za(jn,je)*zb(je,jn)*(1 + propT126*za(j
     &     b,jg)*zb(jg,jb))))) + 2*za(jb,jg)*zb(jb,ju)*zb(je,jb)*(za(jn,jd) + propT
     &     126*za(jn,je)*(za(jc,jd)*zb(jc,je) - za(jd,jg)*zb(jg,je)))) + za(jc,jg)*
     &     za(jd,jg)*(2*propT1267*za(jc,jd)*za(jn,je)*zb(jc,je)**2*zb(jg,ju) + za(j
     &     n,jd)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*(1 - 2*propT1267*za(jn,je)*zb(je,
     &     jn))*zb(jg,ju))) + za(jb,jg)*(-2*propT126*za(jn,je)*za(ju,jd)*zb(jb,ju)*
     &     (za(jd,jg)*zb(je,ju)*zb(jg,je) + za(jc,jd)*zb(jc,je)*(-zb(je,ju) + propT
     &     1267*za(jc,jg)*zb(jc,je)*zb(jg,ju))) + za(jn,jd)*(za(je,jd)*zb(je,jb)*zb
     &     (je,ju) + za(ju,jd)*zb(jb,ju)*(zb(je,ju) + 2*propT126*propT1267*za(jc,jg
     &     )*za(jn,je)*zb(jc,je)*zb(je,jn)*zb(jg,ju))) + za(jn,jd)**2*(zb(jb,ju)*zb
     &     (je,jn) + zb(je,jb)*zb(jn,ju)))))))/(za(jb,jg)*za(jc,jg)*za(jd,jg)*za(jn
     &     ,jg))
       end function streal_heavyResonant_PMMM_P_L2

       function streal_heavyResonant_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34
           real(dp) :: propW16

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_PMMM_M_L2 =
     &     (propW16*propW34*(2*im*imag(anomc4)*(za(jn,jd)*zb(jc,jb)*(zb(jc,ju)*zb(j
     &     e,jb) + zb(jb,ju)*zb(jc,je)*(-1 + 2*propT1267*za(jn,je)*zb(je,jn))) - 2*
     &     propT1267*za(jn,je)*zb(jb,ju)*zb(jc,je)*(za(jc,jd)*zb(jc,jb)*zb(jc,je) +
     &      propT126*(mtsq*za(jd,jg)*zb(jc,je) + (za(jb,jd)*zb(jc,jb) + za(ju,jd)*z
     &     b(jc,ju))*(za(jc,jg)*zb(jc,je) - za(jn,jg)*zb(je,jn)))*zb(jg,jb))) + 2*r
     &     eal(anomc4)*(za(jn,jd)*zb(jc,jb)*(-(zb(jc,ju)*zb(je,jb)) + zb(jb,ju)*zb(
     &     jc,je)*(1 - 2*propT1267*za(jn,je)*zb(je,jn))) + 2*propT1267*za(jn,je)*zb
     &     (jb,ju)*zb(jc,je)*(za(jc,jd)*zb(jc,jb)*zb(jc,je) + propT126*(mtsq*za(jd,
     &     jg)*zb(jc,je) + (za(jb,jd)*zb(jc,jb) + za(ju,jd)*zb(jc,ju))*(za(jc,jg)*z
     &     b(jc,je) - za(jn,jg)*zb(je,jn)))*zb(jg,jb))) + im*Sqrt(mtsq)*propT126*im
     &     ag(anomc2)*zb(jb,ju)*zb(jc,je)*(za(jn,jd)*zb(jc,jb) + propT1267*(za(jn,j
     &     g)*za(ju,jd)*zb(jb,ju) + za(jd,jg)*(za(jn,jc)*zb(jc,jb) + za(jn,je)*zb(j
     &     e,jb)))*zb(jg,jc)) - Sqrt(mtsq)*propT126*real(anomc2)*zb(jb,ju)*zb(jc,je
     &     )*(za(jn,jd)*zb(jc,jb) + propT1267*(za(jn,jg)*za(ju,jd)*zb(jb,ju) + za(j
     &     d,jg)*(za(jn,jc)*zb(jc,jb) + za(jn,je)*zb(je,jb)))*zb(jg,jc))))/(zb(jg,j
     &     b)*zb(jg,jc))
       end function streal_heavyResonant_PMMM_M_L2

       function streal_heavyResonant_MPPP_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MPPP_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MPPP_P_L2 =
     &     (2*Sqrt(mtsq)*propT1267*propW34*real(anomc9)*za(jn,jc)*za(ju,jb)*(za(jb,
     &     jc)*(zb(jd,je) + propT126*za(jb,jg)*zb(jd,jb)*zb(jg,je)) + propT126*za(j
     &     b,jg)*(za(jn,jc)*zb(je,jn)*zb(jg,jd) + za(ju,jc)*zb(jd,ju)*zb(jg,je))))/
     &     (za(jb,jg)*za(jc,jg))
       end function streal_heavyResonant_MPPP_P_L2

       function streal_heavyResonant_MPPP_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyResonant_MPPP_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT126, propT1267, propW34

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propT126 = 1._dp / (s(jb,jd)+s(jb,ju)+s(jd,ju) - mtsq)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_heavyResonant_MPPP_M_L2 =
     &     (2*Sqrt(mtsq)*propW34*real(anomc9)*zb(jd,je)*(propT126*za(ju,jb)*zb(jg,j
     &     b)*(za(jn,jc)*(-zb(jd,jc) + propT1267*(za(jb,jg)*zb(jd,jb) + za(ju,jg)*z
     &     b(jd,ju))*zb(jg,jc)) + za(jn,jg)*zb(jg,jd)) + propT1267*za(jn,jc)*zb(jg,
     &     jc)*(za(ju,jb)*zb(jd,jb) - za(ju,jg)*zb(jg,jd))))/(zb(jg,jb)*zb(jg,jc)*z
     &     b(jg,jd))
       end function streal_heavyResonant_MPPP_M_L2

       function streal_lightWWG_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_MPMM_P_L2 =
     &     (gb**2*propW167*propW34*(im*imag(anomc4) + real(anomc4))*(za(jn,jd)*za(j
     &     u,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*za(ju,jd)*zb(j
     &     c,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*za(ju,jd)*zb(jc,je)*zb(je,ju) + za
     &     (jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*za
     &     (ju,jd)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*za(ju,jd)*zb(je,jb)*zb
     &     (je,ju) - za(jb,je)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(je,ju) + za(jb,jg)*
     &     za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za(jn,jd)*za(ju,jd)*
     &     zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(jg,je) 
     &     - za(jd,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(jg,je) + za(jd,jg)*za(je,jc
     &     )*za(jn,jd)*zb(jc,je)*zb(jg,je) + za(jd,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn
     &     )*zb(jg,je) - za(jd,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(jg,je) + za(jb,
     &     jn)*za(jd,jg)*za(je,jd)*zb(je,jb)*zb(jg,je) - za(jb,je)*za(jd,jg)*za(jn,
     &     jd)*zb(je,jb)*zb(jg,je) + za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jb)*zb(jg,
     &     je) - za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jc)*zb(jg,je) + za(jb,jd)*(2*z
     &     a(jn,jc)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)) - (za(ju,
     &     jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*(za(ju,jn)*zb(je,ju) + za(jn,jg)*zb
     &     (jg,je)) + za(jn,jd)*(za(ju,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(j
     &     e,ju)) + za(jd,jg)*(-(zb(jd,je)*zb(jg,jb)) + zb(jd,jb)*zb(jg,je)))) + za
     &     (jc,jd)*(-((za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc))*(2*za(jb,jn)*zb(j
     &     e,jb) - za(ju,jn)*zb(je,ju) - za(jn,jg)*zb(jg,je))) + za(jn,jd)*(za(ju,j
     &     d)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + za(jd,jg)*(zb(jd,je)*zb
     &     (jg,jc) - zb(jd,jc)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*za(jn,jd)*zb(jc,j
     &     e)*zb(jg,jn) + za(jb,jn)*za(jd,jg)*za(jn,jd)*zb(je,jb)*zb(jg,jn) - za(jb
     &     ,jn)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jn,jb) - za(jb,jn)*za(jd,jg)*za(jn
     &     ,jd)*zb(jg,je)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,je)*zb(jn
     &     ,ju) + za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(jn,ju)))/(ecossin**2*z
     &     a(jd,jg)*za(ju,jg)*zb(jc,jb))
       end function streal_lightWWG_MPMM_P_L2

       function streal_lightWWG_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_MPMM_M_L2 =
     &     (gb**2*propW167*propW34*(im*imag(anomc4) + real(anomc4))*(za(jn,jd)*za(j
     &     u,jb)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,je)*zb(j
     &     d,ju)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(je,ju) + za
     &     (jn,jc)*za(jn,jd)*zb(jc,jn)*zb(jd,ju)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb
     &     (jc,ju)*zb(jd,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jd,ju)*zb(je,jb)*zb
     &     (je,ju) - za(jb,je)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(je,ju) + za(jb,jg)*
     &     za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za(jn,jd)*zb(jd,ju)*
     &     zb(je,ju)*zb(jg,jc) + 2*za(jb,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,ju
     &     ) - za(jb,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,ju) + za(jc,jg)*za(jn,
     &     jd)*zb(jc,ju)*zb(jd,je)*zb(jg,ju) - 2*za(jb,jn)*za(jc,jg)*zb(jc,ju)*zb(j
     &     e,jb)*zb(jg,ju) + za(jn,jg)*za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za
     &     (jb,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za(je,jg)*za(jn,jc)*zb
     &     (jc,je)*zb(je,ju)*zb(jg,ju) + za(je,jc)*za(jn,jg)*zb(jc,je)*zb(je,ju)*zb
     &     (jg,ju) + za(jn,jc)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,ju) - za(jn,jg)*
     &     za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*za(ju,jn)*zb(jc,ju)*
     &     zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(je,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) 
     &     - za(jb,je)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jn,jg
     &     )*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jc,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jc
     &     )*zb(jg,ju) - za(jb,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je)*zb(jg,ju) + za(jc,
     &     jg)*za(jn,jg)*zb(jc,ju)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(2*za(jn,jc)*zb(
     &     jb,ju)*zb(jc,je)*zb(jd,ju) - za(ju,jn)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) + z
     &     a(jn,jd)*zb(jd,ju)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) - za(j
     &     n,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,jb)*zb(je,ju)*zb(j
     &     g,ju)) + za(jc,jd)*(-2*za(jb,jn)*zb(jc,ju)*zb(jd,ju)*zb(je,jb) + za(ju,j
     &     n)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)*(zb(jc,ju)*zb(jd,
     &     je) - zb(jd,jc)*zb(je,ju)) + za(jn,jg)*zb(jc,ju)*zb(jd,ju)*zb(jg,je) - z
     &     a(jn,jg)*zb(jd,jc)*zb(je,ju)*zb(jg,ju)) - za(jb,jn)*za(jn,jd)*zb(jd,ju)*
     &     zb(je,ju)*zb(jn,jb) - za(jb,jn)*za(jn,jg)*zb(je,ju)*zb(jg,ju)*zb(jn,jb) 
     &     - za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(jn,ju) + za(jb,jn)*za(jn,jd
     &     )*zb(jd,ju)*zb(je,jb)*zb(jn,ju) - za(jn,jc)*za(jn,jg)*zb(jc,je)*zb(jg,ju
     &     )*zb(jn,ju) + za(jb,jn)*za(jn,jg)*zb(je,jb)*zb(jg,ju)*zb(jn,ju)))/(ecoss
     &     in**2*zb(jc,jb)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWG_MPMM_M_L2

       function streal_lightWWG_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_PMMM_P_L2 =
     &     (gb**2*propW167*propW34*(im*imag(anomc4) - real(anomc4))*(za(jn,jd)*za(j
     &     u,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*za(ju,jd)*zb(j
     &     c,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*za(ju,jd)*zb(jc,je)*zb(je,ju) + za
     &     (jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*za
     &     (ju,jd)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*za(ju,jd)*zb(je,jb)*zb
     &     (je,ju) - za(jb,je)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(je,ju) + za(jb,jg)*
     &     za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za(jn,jd)*za(ju,jd)*
     &     zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(jg,je) 
     &     - za(jd,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(jg,je) + za(jd,jg)*za(je,jc
     &     )*za(jn,jd)*zb(jc,je)*zb(jg,je) + za(jd,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn
     &     )*zb(jg,je) - za(jd,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(jg,je) + za(jb,
     &     jn)*za(jd,jg)*za(je,jd)*zb(je,jb)*zb(jg,je) - za(jb,je)*za(jd,jg)*za(jn,
     &     jd)*zb(je,jb)*zb(jg,je) + za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jb)*zb(jg,
     &     je) - za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jc)*zb(jg,je) + za(jb,jd)*(2*z
     &     a(jn,jc)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)) - (za(ju,
     &     jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*(za(ju,jn)*zb(je,ju) + za(jn,jg)*zb
     &     (jg,je)) + za(jn,jd)*(za(ju,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(j
     &     e,ju)) + za(jd,jg)*(-(zb(jd,je)*zb(jg,jb)) + zb(jd,jb)*zb(jg,je)))) + za
     &     (jc,jd)*(-((za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc))*(2*za(jb,jn)*zb(j
     &     e,jb) - za(ju,jn)*zb(je,ju) - za(jn,jg)*zb(jg,je))) + za(jn,jd)*(za(ju,j
     &     d)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + za(jd,jg)*(zb(jd,je)*zb
     &     (jg,jc) - zb(jd,jc)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*za(jn,jd)*zb(jc,j
     &     e)*zb(jg,jn) + za(jb,jn)*za(jd,jg)*za(jn,jd)*zb(je,jb)*zb(jg,jn) - za(jb
     &     ,jn)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jn,jb) - za(jb,jn)*za(jd,jg)*za(jn
     &     ,jd)*zb(jg,je)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,je)*zb(jn
     &     ,ju) + za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(jn,ju)))/(ecossin**2*z
     &     a(jb,jc)*za(jd,jg)*za(ju,jg))
       end function streal_lightWWG_PMMM_P_L2

       function streal_lightWWG_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWG_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)

           streal_lightWWG_PMMM_M_L2 =
     &     (gb**2*propW167*propW34*(im*imag(anomc4) - real(anomc4))*(za(jn,jd)*za(j
     &     u,jb)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,je)*zb(j
     &     d,ju)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(je,ju) + za
     &     (jn,jc)*za(jn,jd)*zb(jc,jn)*zb(jd,ju)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb
     &     (jc,ju)*zb(jd,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jd,ju)*zb(je,jb)*zb
     &     (je,ju) - za(jb,je)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(je,ju) + za(jb,jg)*
     &     za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za(jn,jd)*zb(jd,ju)*
     &     zb(je,ju)*zb(jg,jc) + 2*za(jb,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,ju
     &     ) - za(jb,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,ju) + za(jc,jg)*za(jn,
     &     jd)*zb(jc,ju)*zb(jd,je)*zb(jg,ju) - 2*za(jb,jn)*za(jc,jg)*zb(jc,ju)*zb(j
     &     e,jb)*zb(jg,ju) + za(jn,jg)*za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za
     &     (jb,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za(je,jg)*za(jn,jc)*zb
     &     (jc,je)*zb(je,ju)*zb(jg,ju) + za(je,jc)*za(jn,jg)*zb(jc,je)*zb(je,ju)*zb
     &     (jg,ju) + za(jn,jc)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,ju) - za(jn,jg)*
     &     za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*za(ju,jn)*zb(jc,ju)*
     &     zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(je,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) 
     &     - za(jb,je)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jn,jg
     &     )*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jc,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jc
     &     )*zb(jg,ju) - za(jb,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je)*zb(jg,ju) + za(jc,
     &     jg)*za(jn,jg)*zb(jc,ju)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(2*za(jn,jc)*zb(
     &     jb,ju)*zb(jc,je)*zb(jd,ju) - za(ju,jn)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) + z
     &     a(jn,jd)*zb(jd,ju)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) - za(j
     &     n,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,jb)*zb(je,ju)*zb(j
     &     g,ju)) + za(jc,jd)*(-2*za(jb,jn)*zb(jc,ju)*zb(jd,ju)*zb(je,jb) + za(ju,j
     &     n)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)*(zb(jc,ju)*zb(jd,
     &     je) - zb(jd,jc)*zb(je,ju)) + za(jn,jg)*zb(jc,ju)*zb(jd,ju)*zb(jg,je) - z
     &     a(jn,jg)*zb(jd,jc)*zb(je,ju)*zb(jg,ju)) - za(jb,jn)*za(jn,jd)*zb(jd,ju)*
     &     zb(je,ju)*zb(jn,jb) - za(jb,jn)*za(jn,jg)*zb(je,ju)*zb(jg,ju)*zb(jn,jb) 
     &     - za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(jn,ju) + za(jb,jn)*za(jn,jd
     &     )*zb(jd,ju)*zb(je,jb)*zb(jn,ju) - za(jn,jc)*za(jn,jg)*zb(jc,je)*zb(jg,ju
     &     )*zb(jn,ju) + za(jb,jn)*za(jn,jg)*zb(je,jb)*zb(jg,ju)*zb(jn,ju)))/(ecoss
     &     in**2*za(jb,jc)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWG_PMMM_M_L2

       function streal_heavyWWG_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_MPMM_P_L2 =
     &     (gb**2*propW16*propW34*(im*imag(anomc4) + real(anomc4))*za(jb,jc)*(za(jb
     &     ,jc)*(za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(j
     &     c,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju) + za(jn,jc)*za
     &     (jn,jd)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju) + 
     &     za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*zb(je,jb)*
     &     zb(je,ju) + za(jc,jd)*(za(jn,jd)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,
     &     ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju) - za(jn,j
     &     g)*zb(jg,je))) + za(jb,jd)*(2*za(jn,jc)*zb(jb,ju)*zb(jc,je) + za(jn,jd)*
     &     (-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*
     &     zb(je,ju)) + za(jn,jg)*zb(jg,je))) - za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg
     &     ,ju) + za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,ju) - za(jb,jn)*za(jn,jd)*zb(
     &     je,ju)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*z
     &     a(jn,jd)*zb(je,jb)*zb(jn,ju)) + za(jb,jg)*(za(jb,jd)*za(jn,jc)*zb(jb,ju)
     &     *zb(jg,je) - za(je,jd)*za(jn,jc)*zb(je,ju)*zb(jg,je) + za(je,jc)*za(jn,j
     &     d)*zb(je,ju)*zb(jg,je) + za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jn) - za(jn
     &     ,jd)*za(ju,jc)*zb(je,ju)*zb(jg,ju) - za(jd,jg)*za(jn,jc)*zb(jg,je)*zb(jg
     &     ,ju) + za(jc,jd)*(-((za(jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) + za(jn,j
     &     g)*zb(jg,je))*zb(jg,ju)) + za(jn,jc)*(zb(jc,ju)*zb(jg,je) - zb(jc,je)*zb
     &     (jg,ju)) + za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))) - za(j
     &     n,jc)*za(jn,jd)*zb(jg,je)*zb(jn,ju)) + za(jc,jg)*(-(za(jn,jd)*zb(je,ju)*
     &     (za(jb,je)*zb(jg,je) + za(ju,jb)*zb(jg,ju))) + za(jb,jd)*(-((za(jn,jc)*z
     &     b(jc,je) - za(ju,jn)*zb(je,ju) + za(jn,jg)*zb(jg,je))*zb(jg,ju)) + za(jn
     &     ,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)) - za(jb,jn)*(zb(jb,ju)*
     &     zb(jg,je) + zb(je,jb)*zb(jg,ju))) + za(jb,jn)*(-(za(jc,jd)*zb(jc,ju)*zb(
     &     jg,je)) + za(je,jd)*zb(je,ju)*zb(jg,je) - za(jn,jd)*zb(je,ju)*zb(jg,jn) 
     &     + za(jd,jg)*zb(jg,je)*zb(jg,ju) + za(jn,jd)*zb(jg,je)*zb(jn,ju)))))/(eco
     &     ssin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWG_MPMM_P_L2

       function streal_heavyWWG_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_MPMM_M_L2 =
     &     (gb**2*propW16*propW34*(3*im*imag(anomc4)*(za(jb,jn)*za(jc,jd)*za(jc,jg)
     &     *zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,ju)
     &     *zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(je,ju)
     &     *zb(jg,jb) - za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb)
     &      + za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jb,j
     &     n)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*za(jn,j
     &     d)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jb,jn)*za(je,j
     &     d)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za(jn,jd)*zb(je,j
     &     b)*zb(je,ju)*zb(jg,jb) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,j
     &     e)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,j
     &     c) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jc) - za(jb
     &     ,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jn)*za(jc
     &     ,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn
     &     ,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jn,jd)*za(ju,jc)*zb(jb
     &     ,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,ju)*zb(je
     &     ,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,ju)*zb(jg
     &     ,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jc) + za(
     &     jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(
     &     jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jc,jg)*za(
     &     ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(
     &     jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jd,jc)*zb(
     &     je,ju)*zb(jg,jc) + za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(je,ju)*zb(
     &     jg,jc) - za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + z
     &     a(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd) - za(jb,jg)*z
     &     a(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,jg)*za(jd,jg)*z
     &     a(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,jd)*za(jc,jg)*z
     &     b(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*z
     &     b(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,jb)*zb(je,ju)*z
     &     b(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,je) -
     &      za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)
     &     *za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)
     &     *za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(je,jg)*za(jn,jd)
     &     *zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)
     &     *zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*zb(jg,jc)
     &     *zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc)*zb(jg,j
     &     e) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,je) + za(jc
     &     ,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - za(jb,jn)*za(jc
     &     ,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jc)*za(jn
     &     ,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*za(jn,jg)*zb(je
     &     ,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*zb(je,ju)*zb(jg
     &     ,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg
     &     ,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*zb(jg,ju) - za(
     &     jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jc,jg)*za(
     &     jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*za(
     &     ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(
     &     jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,je)*zb(
     &     jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg,jb)*zb(
     &     jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb)*zb(jg,ju) + 2
     &     *za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)
     &     *za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*za(jd,jg)
     &     *za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za(jd,jg)*za(jn,j
     &     c)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jd,j
     &     e)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(je,jb)*zb(jg,j
     &     c)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,j
     &     u) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jc
     &     ,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jb,jn)*za(jc
     &     ,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn
     &     ,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jc
     &     ,jb)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(za(jb,jg)*(za(jn,jc)*zb(jb,ju)*(zb
     &     (jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,je)) + zb(jg,j
     &     b)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)
     &     *(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je)))) + za(jb,jn)*za(jc,jg
     &     )*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) - zb(jc,jb)*zb(jg,je)) + zb(je,jb)*(zb
     &     (jc,ju)*zb(jg,jb) - zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(za(ju,jn)*zb(jc
     &     ,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) - za
     &     (jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(ju,jn)*zb(jc,jb)*zb(je,ju)*zb
     &     (jg,ju) - za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn,jc)*zb(jc,je)*
     &     (zb(jc,ju)*zb(jg,jb) + zb(jb,ju)*zb(jg,jc) - zb(jc,jb)*zb(jg,ju)) + za(j
     &     n,jd)*(-(zb(jc,ju)*zb(jd,je)*zb(jg,jb)) + zb(jd,jc)*zb(je,ju)*zb(jg,jb) 
     &     + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))) - za(jb,jg)*z
     &     a(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) - za(jb,jg)*za(jn,jc)*z
     &     a(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jb,jn)*za(jc,jg)*za(jn,jd)*z
     &     b(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,jb)*z
     &     b(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jg,jc)*z
     &     b(jn,ju) + za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju) +
     &      za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)
     &     *za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)*za(jn,jd)
     &     *za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(jn,jd)*za(jn,jg)
     &     *zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(ju,jb)*zb(jb,ju
     &     )*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(je,ju
     &     ) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(jn,jc)*za(jn,
     &     jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,jb)*zb(jc,
     &     ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) - za(j
     &     b,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)*za(jn,jc)*zb(j
     &     c,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(j
     &     g,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(je,jg)*za
     &     (jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*zb(jc,je)*zb
     &     (je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) - 
     &     za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg)*za(ju,jn)*
     &     zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*zb(je,ju)*
     &     zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jd,jg
     &     )*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(jd,jg)*zb(jb,ju
     &     )*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc
     &     ) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,
     &     jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*zb(je,jb)*zb(je,
     &     ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + za(j
     &     d,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)*za(jn,jg)*zb(j
     &     b,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb)*(zb(jc,ju)*zb
     &     (jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(jc,jb)*zb(je
     &     ,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - 
     &     zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*(2*za(jn,jc)*zb
     &     (jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,ju)*zb(jd,je)
     &     ) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,jb)*zb(je,ju)) +
     &      za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,
     &     je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) + za(jb,jn)*z
     &     a(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb(jc,je)*z
     &     b(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) -
     &      za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(jn,jd)*za(jn,jg)
     &     *zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)
     &     *zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju) - za(jn,j
     &     d)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za(jn,jg)*zb(je,j
     &     b)*zb(jg,jc)*zb(jn,ju))) + 3*real(anomc4)*(za(jb,jn)*za(jc,jd)*za(jc,jg)
     &     *zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,ju)
     &     *zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(je,ju)
     &     *zb(jg,jb) - za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb)
     &      + za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jb,j
     &     n)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*za(jn,j
     &     d)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jb,jn)*za(je,j
     &     d)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za(jn,jd)*zb(je,j
     &     b)*zb(je,ju)*zb(jg,jb) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,j
     &     e)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,j
     &     c) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jc) - za(jb
     &     ,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jn)*za(jc
     &     ,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn
     &     ,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jn,jd)*za(ju,jc)*zb(jb
     &     ,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,ju)*zb(je
     &     ,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,ju)*zb(jg
     &     ,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jc) + za(
     &     jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(
     &     jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jc,jg)*za(
     &     ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(
     &     jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jd,jc)*zb(
     &     je,ju)*zb(jg,jc) + za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(je,ju)*zb(
     &     jg,jc) - za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + z
     &     a(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd) - za(jb,jg)*z
     &     a(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,jg)*za(jd,jg)*z
     &     a(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,jd)*za(jc,jg)*z
     &     b(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*z
     &     b(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,jb)*zb(je,ju)*z
     &     b(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,je) -
     &      za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)
     &     *za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)
     &     *za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(je,jg)*za(jn,jd)
     &     *zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)
     &     *zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*zb(jg,jc)
     &     *zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc)*zb(jg,j
     &     e) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,je) + za(jc
     &     ,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - za(jb,jn)*za(jc
     &     ,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jc)*za(jn
     &     ,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*za(jn,jg)*zb(je
     &     ,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*zb(je,ju)*zb(jg
     &     ,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg
     &     ,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*zb(jg,ju) - za(
     &     jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jc,jg)*za(
     &     jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*za(
     &     ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(
     &     jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,je)*zb(
     &     jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg,jb)*zb(
     &     jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb)*zb(jg,ju) + 2
     &     *za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)
     &     *za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*za(jd,jg)
     &     *za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za(jd,jg)*za(jn,j
     &     c)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jd,j
     &     e)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(je,jb)*zb(jg,j
     &     c)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,j
     &     u) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jc
     &     ,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jb,jn)*za(jc
     &     ,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn
     &     ,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jc
     &     ,jb)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(za(jb,jg)*(za(jn,jc)*zb(jb,ju)*(zb
     &     (jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,je)) + zb(jg,j
     &     b)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)
     &     *(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je)))) + za(jb,jn)*za(jc,jg
     &     )*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) - zb(jc,jb)*zb(jg,je)) + zb(je,jb)*(zb
     &     (jc,ju)*zb(jg,jb) - zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(za(ju,jn)*zb(jc
     &     ,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) - za
     &     (jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(ju,jn)*zb(jc,jb)*zb(je,ju)*zb
     &     (jg,ju) - za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn,jc)*zb(jc,je)*
     &     (zb(jc,ju)*zb(jg,jb) + zb(jb,ju)*zb(jg,jc) - zb(jc,jb)*zb(jg,ju)) + za(j
     &     n,jd)*(-(zb(jc,ju)*zb(jd,je)*zb(jg,jb)) + zb(jd,jc)*zb(je,ju)*zb(jg,jb) 
     &     + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))) - za(jb,jg)*z
     &     a(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) - za(jb,jg)*za(jn,jc)*z
     &     a(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jb,jn)*za(jc,jg)*za(jn,jd)*z
     &     b(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,jb)*z
     &     b(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jg,jc)*z
     &     b(jn,ju) + za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju) +
     &      za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)
     &     *za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)*za(jn,jd)
     &     *za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(jn,jd)*za(jn,jg)
     &     *zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(ju,jb)*zb(jb,ju
     &     )*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(je,ju
     &     ) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(jn,jc)*za(jn,
     &     jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,jb)*zb(jc,
     &     ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) - za(j
     &     b,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)*za(jn,jc)*zb(j
     &     c,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(j
     &     g,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(je,jg)*za
     &     (jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*zb(jc,je)*zb
     &     (je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) - 
     &     za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg)*za(ju,jn)*
     &     zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*zb(je,ju)*
     &     zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jd,jg
     &     )*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(jd,jg)*zb(jb,ju
     &     )*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc
     &     ) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,
     &     jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*zb(je,jb)*zb(je,
     &     ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + za(j
     &     d,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)*za(jn,jg)*zb(j
     &     b,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb)*(zb(jc,ju)*zb
     &     (jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(jc,jb)*zb(je
     &     ,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - 
     &     zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*(2*za(jn,jc)*zb
     &     (jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,ju)*zb(jd,je)
     &     ) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,jb)*zb(je,ju)) +
     &      za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,
     &     je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) + za(jb,jn)*z
     &     a(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb(jc,je)*z
     &     b(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) -
     &      za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(jn,jd)*za(jn,jg)
     &     *zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)
     &     *zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju) - za(jn,j
     &     d)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za(jn,jg)*zb(je,j
     &     b)*zb(jg,jc)*zb(jn,ju))) + 2*(im*imag(anomc7) + real(anomc7))*(za(jb,je)
     &     *za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*za(jn,jd)
     &     *za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jb,jg)*za(jc,jd)*za(jn,jc)
     &     *zb(jb,ju)*zb(jc,je)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)
     &     *zb(jd,je)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)
     &     *zb(jg,jc) - za(jb,jg)*za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc)
     &      + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jb,j
     &     g)*za(jc,jd)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(je,j
     &     d)*za(jn,jc)*zb(je,jb)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jc)*za(jn,j
     &     d)*zb(je,jb)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,j
     &     u)*zb(jg,jc)*zb(jg,je) + za(jb,jd)*(za(jb,jn)*za(jc,jg)*(zb(jb,ju)*zb(jc
     &     ,je) + zb(jc,ju)*zb(je,jb))*zb(jg,jb) - za(jb,jg)*za(jn,jc)*zb(jb,ju)*zb
     &     (je,jb)*zb(jg,jc) + za(jc,jg)*zb(jg,jb)*(za(jn,jc)*zb(jc,je)*zb(jc,ju) +
     &      za(jn,jd)*(-(zb(jc,ju)*zb(jd,je)) + zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-
     &     (za(ju,jn)*zb(je,ju)) + za(jn,jg)*zb(jg,je)))) + za(jb,jg)*za(jd,jg)*za(
     &     jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(
     &     je,ju)*zb(jg,jc)*zb(jn,jb) + za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(
     &     jg,jc)*zb(jn,ju) + za(jb,jn)*(za(jc,jd)*(za(jc,jg)*zb(jc,je)*zb(jc,ju)*z
     &     b(jg,jb) - za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) - za(jc,jg)*zb(jg,jb
     &     )*(za(je,jd)*zb(jc,je)*zb(je,ju) + za(jd,jg)*zb(jc,je)*zb(jg,ju) + za(jn
     &     ,jd)*(-(zb(jc,jn)*zb(je,ju)) + zb(jc,je)*zb(jn,ju)))))))/(3._dp*ecossin**2*
     &     (s(jb,jc) + s(jb,jg) + s(jc,jg))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWG_MPMM_M_L2

       function streal_heavyWWG_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_PMMM_P_L2 =
     &     (gb**2*propW16*propW34*(3*im*imag(anomc4)*(-(za(jb,jn)*za(jc,jd)*za(jc,j
     &     g)*zb(jc,je)*zb(jc,ju)*zb(jg,jb)) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,
     &     ju)*zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(je,
     &     ju)*zb(jg,jb) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(je,ju)*zb(jg,
     &     jb) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(j
     &     b,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za(j
     &     n,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jb,jn)*za(j
     &     e,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za(jn,jd)*zb(j
     &     e,jb)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(j
     &     c,je)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(j
     &     g,jc) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jc) + za
     &     (jb,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jn)*za
     &     (jc,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za
     &     (jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(jn,jd)*za(ju,jc)*zb
     &     (jb,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,ju)*zb
     &     (je,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,ju)*zb
     &     (jg,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jc) + 
     &     za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*
     &     za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jc,jg)*
     &     za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*
     &     zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jd,jc)*
     &     zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(je,ju)*
     &     zb(jg,jc) + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) 
     &     + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd) - za(jb,jg
     &     )*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,jg)*za(jd,jg
     &     )*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,jd)*za(jc,jg
     &     )*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb
     &     )*zb(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,jb)*zb(je,ju
     &     )*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,je
     &     ) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,
     &     jg)*za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)*za(jc,
     &     jd)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(je,jg)*za(jn,
     &     jd)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jg)*zb(je,
     &     ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*zb(jg,
     &     jc)*zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc)*zb(j
     &     g,je) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,je) + za
     &     (jc,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - za(jb,jn)*za
     &     (jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jc)*za
     &     (jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*za(jn,jg)*zb
     &     (je,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*zb(je,ju)*zb
     &     (jg,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb
     &     (jg,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*zb(jg,ju) - 
     &     za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jc,jg)*
     &     za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*
     &     za(ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(ju,jn)*
     &     zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,je)*
     &     zb(jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg,jb)*
     &     zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb)*zb(jg,ju) 
     &     + 2*za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju) - za(jb,
     &     jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*za(jd,
     &     jg)*za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za(jd,jg)*za(j
     &     n,jc)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(j
     &     d,je)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(je,jb)*zb(j
     &     g,jc)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(j
     &     g,ju) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za
     &     (jc,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jb,jn)*za
     &     (jc,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za
     &     (jn,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb
     &     (jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jd)*(-(za(jb,jg)*(za(jn,jc)*zb(jb,ju
     &     )*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,je)) + zb
     &     (jg,jb)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(j
     &     b,ju)*(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je))))) + za(jb,jn)*za
     &     (jc,jg)*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) + zb(jc,jb)*zb(jg,je)) + zb(je,j
     &     b)*(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(za(ju,jn)
     &     *zb(jc,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je
     &     ) + za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) - za(ju,jn)*zb(jc,jb)*zb(je,
     &     ju)*zb(jg,ju) + za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn,jc)*zb(j
     &     c,je)*(zb(jc,ju)*zb(jg,jb) - zb(jb,ju)*zb(jg,jc) + zb(jc,jb)*zb(jg,ju)) 
     &     - za(jn,jd)*(zb(jc,ju)*zb(jd,je)*zb(jg,jb) - zb(jd,jc)*zb(je,ju)*zb(jg,j
     &     b) + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))) - za(jb,jg
     &     )*za(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) + za(jb,jg)*za(jn,jc
     &     )*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + za(jb,jn)*za(jc,jg)*za(jn,jd
     &     )*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,jb
     &     )*zb(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jg,jc
     &     )*zb(jn,ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju
     &     ) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,
     &     jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)*za(jn,
     &     jd)*za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(jn,jd)*za(jn,
     &     jg)*zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(ju,jb)*zb(jb
     &     ,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(je
     &     ,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(jn,jc)*za(
     &     jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,jb)*zb(
     &     jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) - z
     &     a(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)*za(jn,jc)*z
     &     b(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*z
     &     b(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(je,jg)
     &     *za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*zb(jc,je)
     &     *zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,jb)
     &      - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg)*za(ju,j
     &     n)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*zb(je,j
     &     u)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jd
     &     ,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(jd,jg)*zb(jb
     &     ,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg
     &     ,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(
     &     jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*zb(je,jb)*zb(
     &     je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + z
     &     a(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)*za(jn,jg)*z
     &     b(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb)*(zb(jc,ju)
     &     *zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(jc,jb)*zb
     &     (je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,je)*zb(jg,jb)
     &      - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*(2*za(jn,jc)
     &     *zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,ju)*zb(jd,
     &     je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,jb)*zb(je,ju)
     &     ) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(
     &     jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) + za(jb,jn
     &     )*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb(jc,je
     &     )*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju
     &     ) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(jn,jd)*za(jn,
     &     jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,
     &     je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju) - za(j
     &     n,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za(jn,jg)*zb(j
     &     e,jb)*zb(jg,jc)*zb(jn,ju))) - 3*real(anomc4)*(-(za(jb,jn)*za(jc,jd)*za(j
     &     c,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb)) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(
     &     jc,ju)*zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(
     &     je,ju)*zb(jg,jb) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(je,ju)*zb(
     &     jg,jb) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - z
     &     a(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*z
     &     a(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jb,jn)*z
     &     a(je,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za(jn,jd)*z
     &     b(je,jb)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*z
     &     b(jc,je)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*z
     &     b(jg,jc) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jc) +
     &      za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jn)
     &     *za(jc,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(jc,jd)
     &     *za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(jn,jd)*za(ju,jc)
     &     *zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,ju)
     &     *zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,ju)
     &     *zb(jg,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jc)
     &      + za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) - za(jc,j
     &     g)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jc,j
     &     g)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,j
     &     d)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jd,j
     &     c)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(je,j
     &     u)*zb(jg,jc) + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,j
     &     c) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd) - za(jb
     &     ,jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,jg)*za(jd
     &     ,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,jd)*za(jc
     &     ,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc
     &     ,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,jb)*zb(je
     &     ,ju)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,ju)*zb(jg
     &     ,je) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(
     &     jb,jg)*za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)*za(
     &     jc,jd)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(je,jg)*za(
     &     jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jg)*zb(
     &     je,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*zb(
     &     jg,jc)*zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc)*z
     &     b(jg,je) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,je) +
     &      za(jc,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - za(jb,jn)
     &     *za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jc)
     &     *za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*za(jn,jg)
     &     *zb(je,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*zb(je,ju)
     &     *zb(jg,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)
     &     *zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*zb(jg,ju)
     &      - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jc,j
     &     g)*za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(jn,j
     &     d)*za(ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(ju,j
     &     n)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,j
     &     e)*zb(jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg,j
     &     b)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb)*zb(jg,j
     &     u) + 2*za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju) - za(
     &     jb,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*za(
     &     jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za(jd,jg)*z
     &     a(jn,jc)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(jn,jd)*z
     &     b(jd,je)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(je,jb)*z
     &     b(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*z
     &     b(jg,ju) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) +
     &      za(jc,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jb,jn)
     &     *za(jc,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jd,jg)
     &     *za(jn,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*za(jn,jg)
     &     *zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jd)*(-(za(jb,jg)*(za(jn,jc)*zb(jb
     &     ,ju)*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,je)) +
     &      zb(jg,jb)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + z
     &     b(jb,ju)*(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je))))) + za(jb,jn)
     &     *za(jc,jg)*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) + zb(jc,jb)*zb(jg,je)) + zb(j
     &     e,jb)*(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(za(ju,
     &     jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg
     &     ,je) + za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) - za(ju,jn)*zb(jc,jb)*zb(
     &     je,ju)*zb(jg,ju) + za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn,jc)*z
     &     b(jc,je)*(zb(jc,ju)*zb(jg,jb) - zb(jb,ju)*zb(jg,jc) + zb(jc,jb)*zb(jg,ju
     &     )) - za(jn,jd)*(zb(jc,ju)*zb(jd,je)*zb(jg,jb) - zb(jd,jc)*zb(je,ju)*zb(j
     &     g,jb) + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))) - za(jb
     &     ,jg)*za(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) + za(jb,jg)*za(jn
     &     ,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + za(jb,jn)*za(jc,jg)*za(jn
     &     ,jd)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je
     &     ,jb)*zb(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jg
     &     ,jc)*zb(jn,ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn
     &     ,ju) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(
     &     jb,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)*za(
     &     jn,jd)*za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(jn,jd)*za(
     &     jn,jg)*zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(ju,jb)*zb
     &     (jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb
     &     (je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(jn,jc)*
     &     za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,jb)*
     &     zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) 
     &     - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)*za(jn,jc
     &     )*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je
     &     )*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(je,
     &     jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*zb(jc,
     &     je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(jg,
     &     jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg)*za(j
     &     u,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*zb(j
     &     e,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za
     &     (jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(jd,jg)*zb
     &     (jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb
     &     (jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*
     &     za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*zb(je,jb)*
     &     zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,jc) 
     &     + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)*za(jn,jg
     &     )*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb)*(zb(jc,
     &     ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(jc,jb)
     &     *zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,je)*zb(jg,
     &     jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*(2*za(jn,
     &     jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,ju)*zb(
     &     jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,jb)*zb(je,
     &     ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*
     &     zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) + za(jb
     &     ,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb(jc
     &     ,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jg
     &     ,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(jn,jd)*za(
     &     jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(
     &     jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju) - z
     &     a(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za(jn,jg)*z
     &     b(je,jb)*zb(jg,jc)*zb(jn,ju))) - 2*(im*imag(anomc7) - real(anomc7))*(za(
     &     jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*za(
     &     jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jb,jg)*za(jc,jd)*za(
     &     jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(
     &     jb,ju)*zb(jd,je)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(
     &     je,jb)*zb(jg,jc) - za(jb,jg)*za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(
     &     jg,jc) + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - z
     &     a(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*z
     &     a(je,jd)*za(jn,jc)*zb(je,jb)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jc)*z
     &     a(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jg)*z
     &     b(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jb,jd)*(za(jb,jn)*za(jc,jg)*(zb(jb,ju)
     &     *zb(jc,je) + zb(jc,ju)*zb(je,jb))*zb(jg,jb) - za(jb,jg)*za(jn,jc)*zb(jb,
     &     ju)*zb(je,jb)*zb(jg,jc) + za(jc,jg)*zb(jg,jb)*(za(jn,jc)*zb(jc,je)*zb(jc
     &     ,ju) + za(jn,jd)*(-(zb(jc,ju)*zb(jd,je)) + zb(jd,jc)*zb(je,ju)) + zb(jc,
     &     ju)*(-(za(ju,jn)*zb(je,ju)) + za(jn,jg)*zb(jg,je)))) + za(jb,jg)*za(jd,j
     &     g)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jn,jc)*za(jn,j
     &     d)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,j
     &     b)*zb(jg,jc)*zb(jn,ju) + za(jb,jn)*(za(jc,jd)*(za(jc,jg)*zb(jc,je)*zb(jc
     &     ,ju)*zb(jg,jb) - za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) - za(jc,jg)*zb
     &     (jg,jb)*(za(je,jd)*zb(jc,je)*zb(je,ju) + za(jd,jg)*zb(jc,je)*zb(jg,ju) +
     &      za(jn,jd)*(-(zb(jc,jn)*zb(je,ju)) + zb(jc,je)*zb(jn,ju)))))))/(3._dp*ecoss
     &     in**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWG_PMMM_P_L2

       function streal_heavyWWG_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWG_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)

           streal_heavyWWG_PMMM_M_L2 =
     &     (gb**2*propW16*propW34*(im*imag(anomc4) - real(anomc4))*zb(jc,jb)*(za(jn
     &     ,jd)*za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc
     &     ,jb)*zb(jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je
     &     ,ju) + za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(
     &     ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(
     &     je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + z
     &     a(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*z
     &     b(jc,ju)*zb(jd,je)*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*z
     &     b(jg,jb) + za(je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)
     &     *za(jn,jg)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)
     &     *zb(je,ju)*zb(jg,jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb)
     &      - za(jd,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,j
     &     d)*zb(jd,jc)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,j
     &     e)*zb(jg,jc) - za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb
     &     ,jn)*za(jd,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb
     &     ,ju)*zb(je,ju)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg
     &     ,jc) + za(jd,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(
     &     jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(
     &     je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + z
     &     a(jd,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*
     &     zb(jc,jb)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za
     &     (jb,jn)*zb(jc,jb)*zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*
     &     (zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za
     &     (jb,jd)*(2*za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)
     &     *(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)
     &     *zb(jc,jb)*zb(je,ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg
     &     ,jc) + zb(jc,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*
     &     zb(jg,ju) + za(jb,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg
     &     )*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb
     &     )*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb
     &     ) - za(jn,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,
     &     jd)*zb(jc,jb)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,
     &     jb)*zb(jn,ju) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(j
     &     n,jd)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jn,ju)))/(ecossin**2*(s(jb,jc) + 
     &     s(jb,jg) + s(jc,jg))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWG_PMMM_M_L2

       function streal_lightWWZ_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_MPMM_P_L2 =
     &     (gw**2*propW167*propW34*propZ25*(im*imag(anomc4) + real(anomc4))*za(jb,j
     &     c)*(za(jn,jd)*za(ju,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,
     &     jc)*za(ju,jd)*zb(jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*za(ju,jd)*zb(jc,
     &     je)*zb(je,ju) + za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,jn)*zb(je,ju) - za(j
     &     n,jd)*za(ju,jc)*za(ju,jd)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*za(j
     &     u,jd)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(j
     &     e,ju) + za(jb,jg)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za
     &     (jn,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*za(ju,jb)*zb
     &     (jb,ju)*zb(jg,je) - za(jd,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(jg,je) + 
     &     za(jd,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(jg,je) + za(jd,jg)*za(jn,jc)*
     &     za(jn,jd)*zb(jc,jn)*zb(jg,je) - za(jd,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*
     &     zb(jg,je) + za(jb,jn)*za(jd,jg)*za(je,jd)*zb(je,jb)*zb(jg,je) - za(jb,je
     &     )*za(jd,jg)*za(jn,jd)*zb(je,jb)*zb(jg,je) + za(jb,jg)*za(jd,jg)*za(jn,jd
     &     )*zb(jg,jb)*zb(jg,je) - za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jc)*zb(jg,je
     &     ) + za(jb,jd)*(2*za(jn,jc)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb
     &     (jg,jb)) - (za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*(za(ju,jn)*zb(je,
     &     ju) + za(jn,jg)*zb(jg,je)) + za(jn,jd)*(za(ju,jd)*(-(zb(jb,ju)*zb(jd,je)
     &     ) + zb(jd,jb)*zb(je,ju)) + za(jd,jg)*(-(zb(jd,je)*zb(jg,jb)) + zb(jd,jb)
     &     *zb(jg,je)))) + za(jc,jd)*(-((za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc))
     &     *(2*za(jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) - za(jn,jg)*zb(jg,je))) + 
     &     za(jn,jd)*(za(ju,jd)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + za(jd
     &     ,jg)*(zb(jd,je)*zb(jg,jc) - zb(jd,jc)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)
     &     *za(jn,jd)*zb(jc,je)*zb(jg,jn) + za(jb,jn)*za(jd,jg)*za(jn,jd)*zb(je,jb)
     &     *zb(jg,jn) - za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jn,jb) - za(jb,j
     &     n)*za(jd,jg)*za(jn,jd)*zb(jg,je)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*za(ju,j
     &     d)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(jn,j
     &     u)))/(ecossin**2*za(jd,jg)*za(ju,jg))
       end function streal_lightWWZ_MPMM_P_L2

       function streal_lightWWZ_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_MPMM_M_L2 =
     &     (gw**2*propW167*propW34*propZ25*(im*imag(anomc4) + real(anomc4))*za(jb,j
     &     c)*(za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) - za(je,jd)*za(jn,
     &     jc)*zb(jc,je)*zb(jd,ju)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,je)*zb(jd,
     &     ju)*zb(je,ju) + za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(jd,ju)*zb(je,ju) - za(j
     &     n,jd)*za(ju,jc)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(j
     &     d,ju)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(j
     &     e,ju) + za(jb,jg)*za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za
     &     (jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jc) + 2*za(jb,jg)*za(jn,jc)*zb(jb,ju)*
     &     zb(jc,je)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,ju) 
     &     + za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,ju) - 2*za(jb,jn)*za(jc,
     &     jg)*zb(jc,ju)*zb(je,jb)*zb(jg,ju) + za(jn,jg)*za(ju,jb)*zb(jb,ju)*zb(je,
     &     ju)*zb(jg,ju) - za(jb,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za(j
     &     e,jg)*za(jn,jc)*zb(jc,je)*zb(je,ju)*zb(jg,ju) + za(je,jc)*za(jn,jg)*zb(j
     &     c,je)*zb(je,ju)*zb(jg,ju) + za(jn,jc)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(j
     &     g,ju) - za(jn,jg)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*za
     &     (ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(je,jg)*zb(je,jb)*zb
     &     (je,ju)*zb(jg,ju) - za(jb,je)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) + 
     &     za(jb,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jc,jg)*za(jn,jg)*
     &     zb(je,ju)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je)*
     &     zb(jg,ju) + za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,je)*zb(jg,ju) + za(jb,jd
     &     )*(2*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jd,ju) - za(ju,jn)*zb(jb,ju)*zb(jd
     &     ,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)
     &     *zb(je,ju)) - za(jn,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,
     &     jb)*zb(je,ju)*zb(jg,ju)) + za(jc,jd)*(-2*za(jb,jn)*zb(jc,ju)*zb(jd,ju)*z
     &     b(je,jb) + za(ju,jn)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)
     &     *(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + za(jn,jg)*zb(jc,ju)*zb(jd
     &     ,ju)*zb(jg,je) - za(jn,jg)*zb(jd,jc)*zb(je,ju)*zb(jg,ju)) - za(jb,jn)*za
     &     (jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jn,jb) - za(jb,jn)*za(jn,jg)*zb(je,ju)*zb
     &     (jg,ju)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(jn,ju) + 
     &     za(jb,jn)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(jn,ju) - za(jn,jc)*za(jn,jg)*
     &     zb(jc,je)*zb(jg,ju)*zb(jn,ju) + za(jb,jn)*za(jn,jg)*zb(je,jb)*zb(jg,ju)*
     &     zb(jn,ju)))/(ecossin**2*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWZ_MPMM_M_L2

       function streal_lightWWZ_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_PMMM_P_L2 =
     &     (gw**2*propW167*propW34*propZ25*(im*imag(anomc4) - real(anomc4))*zb(jc,j
     &     b)*(za(jn,jd)*za(ju,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn,
     &     jc)*za(ju,jd)*zb(jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*za(ju,jd)*zb(jc,
     &     je)*zb(je,ju) + za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,jn)*zb(je,ju) - za(j
     &     n,jd)*za(ju,jc)*za(ju,jd)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*za(j
     &     u,jd)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(j
     &     e,ju) + za(jb,jg)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za
     &     (jn,jd)*za(ju,jd)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*za(ju,jb)*zb
     &     (jb,ju)*zb(jg,je) - za(jd,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(jg,je) + 
     &     za(jd,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(jg,je) + za(jd,jg)*za(jn,jc)*
     &     za(jn,jd)*zb(jc,jn)*zb(jg,je) - za(jd,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*
     &     zb(jg,je) + za(jb,jn)*za(jd,jg)*za(je,jd)*zb(je,jb)*zb(jg,je) - za(jb,je
     &     )*za(jd,jg)*za(jn,jd)*zb(je,jb)*zb(jg,je) + za(jb,jg)*za(jd,jg)*za(jn,jd
     &     )*zb(jg,jb)*zb(jg,je) - za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(jg,jc)*zb(jg,je
     &     ) + za(jb,jd)*(2*za(jn,jc)*zb(jc,je)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb
     &     (jg,jb)) - (za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*(za(ju,jn)*zb(je,
     &     ju) + za(jn,jg)*zb(jg,je)) + za(jn,jd)*(za(ju,jd)*(-(zb(jb,ju)*zb(jd,je)
     &     ) + zb(jd,jb)*zb(je,ju)) + za(jd,jg)*(-(zb(jd,je)*zb(jg,jb)) + zb(jd,jb)
     &     *zb(jg,je)))) + za(jc,jd)*(-((za(ju,jd)*zb(jc,ju) + za(jd,jg)*zb(jg,jc))
     &     *(2*za(jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) - za(jn,jg)*zb(jg,je))) + 
     &     za(jn,jd)*(za(ju,jd)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + za(jd
     &     ,jg)*(zb(jd,je)*zb(jg,jc) - zb(jd,jc)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)
     &     *za(jn,jd)*zb(jc,je)*zb(jg,jn) + za(jb,jn)*za(jd,jg)*za(jn,jd)*zb(je,jb)
     &     *zb(jg,jn) - za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(je,ju)*zb(jn,jb) - za(jb,j
     &     n)*za(jd,jg)*za(jn,jd)*zb(jg,je)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*za(ju,j
     &     d)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(je,jb)*zb(jn,j
     &     u)))/(ecossin**2*za(jd,jg)*za(ju,jg))
       end function streal_lightWWZ_PMMM_P_L2

       function streal_lightWWZ_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightWWZ_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW167, propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightWWZ_PMMM_M_L2 =
     &     (gw**2*propW167*propW34*propZ25*(im*imag(anomc4) - real(anomc4))*zb(jc,j
     &     b)*(za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(jd,ju)*zb(je,ju) - za(je,jd)*za(jn,
     &     jc)*zb(jc,je)*zb(jd,ju)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,je)*zb(jd,
     &     ju)*zb(je,ju) + za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(jd,ju)*zb(je,ju) - za(j
     &     n,jd)*za(ju,jc)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(j
     &     d,ju)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(j
     &     e,ju) + za(jb,jg)*za(jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*za
     &     (jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jg,jc) + 2*za(jb,jg)*za(jn,jc)*zb(jb,ju)*
     &     zb(jc,je)*zb(jg,ju) - za(jb,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,ju) 
     &     + za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,ju) - 2*za(jb,jn)*za(jc,
     &     jg)*zb(jc,ju)*zb(je,jb)*zb(jg,ju) + za(jn,jg)*za(ju,jb)*zb(jb,ju)*zb(je,
     &     ju)*zb(jg,ju) - za(jb,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,ju) - za(j
     &     e,jg)*za(jn,jc)*zb(jc,je)*zb(je,ju)*zb(jg,ju) + za(je,jc)*za(jn,jg)*zb(j
     &     c,je)*zb(je,ju)*zb(jg,ju) + za(jn,jc)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(j
     &     g,ju) - za(jn,jg)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*za
     &     (ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(je,jg)*zb(je,jb)*zb
     &     (je,ju)*zb(jg,ju) - za(jb,je)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,ju) + 
     &     za(jb,jg)*za(jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jc,jg)*za(jn,jg)*
     &     zb(je,ju)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,je)*
     &     zb(jg,ju) + za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,je)*zb(jg,ju) + za(jb,jd
     &     )*(2*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jd,ju) - za(ju,jn)*zb(jb,ju)*zb(jd
     &     ,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)
     &     *zb(je,ju)) - za(jn,jg)*zb(jb,ju)*zb(jd,ju)*zb(jg,je) + za(jn,jg)*zb(jd,
     &     jb)*zb(je,ju)*zb(jg,ju)) + za(jc,jd)*(-2*za(jb,jn)*zb(jc,ju)*zb(jd,ju)*z
     &     b(je,jb) + za(ju,jn)*zb(jc,ju)*zb(jd,ju)*zb(je,ju) + za(jn,jd)*zb(jd,ju)
     &     *(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + za(jn,jg)*zb(jc,ju)*zb(jd
     &     ,ju)*zb(jg,je) - za(jn,jg)*zb(jd,jc)*zb(je,ju)*zb(jg,ju)) - za(jb,jn)*za
     &     (jn,jd)*zb(jd,ju)*zb(je,ju)*zb(jn,jb) - za(jb,jn)*za(jn,jg)*zb(je,ju)*zb
     &     (jg,ju)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,ju)*zb(jn,ju) + 
     &     za(jb,jn)*za(jn,jd)*zb(jd,ju)*zb(je,jb)*zb(jn,ju) - za(jn,jc)*za(jn,jg)*
     &     zb(jc,je)*zb(jg,ju)*zb(jn,ju) + za(jb,jn)*za(jn,jg)*zb(je,jb)*zb(jg,ju)*
     &     zb(jn,ju)))/(ecossin**2*zb(jg,jd)*zb(jg,ju))
       end function streal_lightWWZ_PMMM_M_L2

       function streal_heavyWWZ_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_MPMM_P_L2 =
     &     (gw**2*propW16*propW34*propZ257*(im*imag(anomc4) + real(anomc4))*za(jb,j
     &     c)*(za(jb,jc)*(za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(je,ju) - za(je,jd)*za(jn
     &     ,jc)*zb(jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju) + za(
     &     jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(
     &     je,ju) + za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*z
     &     b(je,jb)*zb(je,ju) + za(jc,jd)*(za(jn,jd)*(zb(jc,ju)*zb(jd,je) - zb(jd,j
     &     c)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju) 
     &     - za(jn,jg)*zb(jg,je))) + za(jb,jd)*(2*za(jn,jc)*zb(jb,ju)*zb(jc,je) + z
     &     a(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(z
     &     a(ju,jn)*zb(je,ju)) + za(jn,jg)*zb(jg,je))) - za(jd,jg)*za(jn,jc)*zb(jc,
     &     je)*zb(jg,ju) + za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,ju) - za(jb,jn)*za(j
     &     n,jd)*zb(je,ju)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jn,ju) + za
     &     (jb,jn)*za(jn,jd)*zb(je,jb)*zb(jn,ju)) + za(jb,jg)*(za(jb,jd)*za(jn,jc)*
     &     zb(jb,ju)*zb(jg,je) - za(je,jd)*za(jn,jc)*zb(je,ju)*zb(jg,je) + za(je,jc
     &     )*za(jn,jd)*zb(je,ju)*zb(jg,je) + za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jn
     &     ) - za(jn,jd)*za(ju,jc)*zb(je,ju)*zb(jg,ju) - za(jd,jg)*za(jn,jc)*zb(jg,
     &     je)*zb(jg,ju) + za(jc,jd)*(-((za(jb,jn)*zb(je,jb) - za(ju,jn)*zb(je,ju) 
     &     + za(jn,jg)*zb(jg,je))*zb(jg,ju)) + za(jn,jc)*(zb(jc,ju)*zb(jg,je) - zb(
     &     jc,je)*zb(jg,ju)) + za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)
     &     )) - za(jn,jc)*za(jn,jd)*zb(jg,je)*zb(jn,ju)) + za(jc,jg)*(-(za(jn,jd)*z
     &     b(je,ju)*(za(jb,je)*zb(jg,je) + za(ju,jb)*zb(jg,ju))) + za(jb,jd)*(-((za
     &     (jn,jc)*zb(jc,je) - za(ju,jn)*zb(je,ju) + za(jn,jg)*zb(jg,je))*zb(jg,ju)
     &     ) + za(jn,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)) - za(jb,jn)*(z
     &     b(jb,ju)*zb(jg,je) + zb(je,jb)*zb(jg,ju))) + za(jb,jn)*(-(za(jc,jd)*zb(j
     &     c,ju)*zb(jg,je)) + za(je,jd)*zb(je,ju)*zb(jg,je) - za(jn,jd)*zb(je,ju)*z
     &     b(jg,jn) + za(jd,jg)*zb(jg,je)*zb(jg,ju) + za(jn,jd)*zb(jg,je)*zb(jn,ju)
     &     ))))/(ecossin**2*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWZ_MPMM_P_L2

       function streal_heavyWWZ_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_MPMM_M_L2 =
     &     (propW16*propW34*propZ257*(3*gw**2*im*imag(anomc4)*(za(jb,jn)*za(jc,jd)*
     &     za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*za(jb,jn)*za(jc,jd)*
     &     zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(jb,ju)*
     &     zb(je,ju)*zb(jg,jb) - za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(je,ju)*
     &     zb(jg,jb) + za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) 
     &     + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) + za(jc,jg
     &     )*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jb,jn
     &     )*za(je,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za(jn,jd
     &     )*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju
     &     )*zb(jc,je)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je
     &     )*zb(jg,jc) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jc
     &     ) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jb,
     &     jn)*za(jc,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jc,
     &     jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jn,jd)*za(ju,
     &     jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(jb,
     &     ju)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,
     &     ju)*zb(jg,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,
     &     jc) + za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) - za(j
     &     c,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(j
     &     c,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(j
     &     n,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(j
     &     d,jc)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*zb(j
     &     e,ju)*zb(jg,jc) - za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(j
     &     g,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd) - za
     &     (jb,jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,jg)*za
     &     (jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,jd)*za
     &     (jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb
     &     (jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,jb)*zb
     &     (je,ju)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,ju)*zb
     &     (jg,je) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + 
     &     za(jb,jg)*za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,jg)*
     &     za(jc,jd)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(je,jg)*
     &     za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jg)*
     &     zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*
     &     zb(jg,jc)*zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jc
     &     )*zb(jg,je) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,je
     &     ) + za(jc,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - za(jb,
     &     jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,
     &     jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*za(jn,
     &     jg)*zb(je,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*zb(je,
     &     ju)*zb(jg,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,
     &     je)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*zb(jg,
     &     ju) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(j
     &     c,jg)*za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*za(j
     &     n,jd)*za(ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(j
     &     u,jn)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(j
     &     c,je)*zb(jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(j
     &     g,jb)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb)*zb(j
     &     g,ju) + 2*za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju) - 
     &     za(jb,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*
     &     za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za(jd,jg
     &     )*za(jn,jc)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(jn,jd
     &     )*zb(jd,je)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(je,jb
     &     )*zb(jg,jc)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc
     &     )*zb(jg,ju) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,ju
     &     ) + za(jc,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(jb,
     &     jn)*za(jc,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jd,
     &     jg)*za(jn,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*za(jn,
     &     jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(za(jb,jg)*(za(jn,jc)*zb(j
     &     b,ju)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,je)) 
     &     + zb(jg,jb)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + 
     &     zb(jb,ju)*(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je)))) + za(jb,jn)
     &     *za(jc,jg)*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) - zb(jc,jb)*zb(jg,je)) + zb(j
     &     e,jb)*(zb(jc,ju)*zb(jg,jb) - zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(za(ju,
     &     jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg
     &     ,je) - za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(ju,jn)*zb(jc,jb)*zb(
     &     je,ju)*zb(jg,ju) - za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn,jc)*z
     &     b(jc,je)*(zb(jc,ju)*zb(jg,jb) + zb(jb,ju)*zb(jg,jc) - zb(jc,jb)*zb(jg,ju
     &     )) + za(jn,jd)*(-(zb(jc,ju)*zb(jd,je)*zb(jg,jb)) + zb(jd,jc)*zb(je,ju)*z
     &     b(jg,jb) + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))) - za
     &     (jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) - za(jb,jg)*za
     &     (jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jb,jn)*za(jc,jg)*za
     &     (jn,jd)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,jd)*zb
     &     (je,jb)*zb(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb
     &     (jg,jc)*zb(jn,ju) + za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb
     &     (jn,ju) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - 
     &     za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,jg)*
     &     za(jn,jd)*za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(jn,jd)*
     &     za(jn,jg)*zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(ju,jb)
     &     *zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)
     &     *zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(jn,j
     &     c)*za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(jc,j
     &     b)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,j
     &     u) - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)*za(jn
     &     ,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb(jd
     &     ,je)*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(
     &     je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*zb(
     &     jc,je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*zb(
     &     jg,jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg)*z
     &     a(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd,jc)*z
     &     b(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) -
     &      za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(jd,jg)
     &     *zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,ju)
     &     *zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jd,j
     &     g)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*zb(je,j
     &     b)*zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(jg,j
     &     c) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)*za(jn
     &     ,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb)*(zb(
     &     jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(jc,
     &     jb)*zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,je)*zb(
     &     jg,jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*(2*za(
     &     jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,ju)*
     &     zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,jb)*zb(
     &     je,ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,j
     &     b)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) + za
     &     (jb,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)*zb
     &     (jc,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb
     &     (jg,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(jn,jd)*
     &     za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,jb)*
     &     zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju) 
     &     - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za(jn,jg
     &     )*zb(je,jb)*zb(jg,jc)*zb(jn,ju))) + 3*gw**2*real(anomc4)*(za(jb,jn)*za(j
     &     c,jd)*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*za(jb,jn)*za(j
     &     c,jd)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(j
     &     b,ju)*zb(je,ju)*zb(jg,jb) - za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(j
     &     e,ju)*zb(jg,jb) + za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(j
     &     g,jb) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) + za
     &     (jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za
     &     (jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za
     &     (jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb
     &     (jb,ju)*zb(jc,je)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb
     &     (jd,je)*zb(jg,jc) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb
     &     (jg,jc) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - 
     &     za(jb,jn)*za(jc,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*
     &     za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(jb,jg)*za(jn,jd)*
     &     za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(ju,jn)*
     &     zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*
     &     zb(je,ju)*zb(jg,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*
     &     zb(jg,jc) + za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) 
     &     - za(jc,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd
     &     )*za(jc,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd
     &     )*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd
     &     )*zb(jd,jc)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb
     &     )*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju
     &     )*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd
     &     ) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,
     &     jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,
     &     jd)*za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,
     &     jc)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,
     &     jb)*zb(je,ju)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,
     &     ju)*zb(jg,je) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,
     &     je) + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(j
     &     b,jg)*za(jc,jd)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(j
     &     e,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(j
     &     n,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(j
     &     b,ju)*zb(jg,jc)*zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb
     &     (jg,jc)*zb(jg,je) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb
     &     (jg,je) + za(jc,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - 
     &     za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*
     &     za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*
     &     za(jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*
     &     zb(je,ju)*zb(jg,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*
     &     zb(jc,je)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*
     &     zb(jg,ju) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) 
     &     - za(jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg
     &     )*za(jn,jd)*za(ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd
     &     )*za(ju,jn)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg
     &     )*zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je
     &     )*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb
     &     )*zb(jg,ju) + 2*za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,
     &     ju) - za(jb,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(j
     &     b,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za
     &     (jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za
     &     (jn,jd)*zb(jd,je)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb
     &     (je,jb)*zb(jg,jc)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb
     &     (jg,jc)*zb(jg,ju) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb
     &     (jg,ju) + za(jc,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + 
     &     za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*
     &     za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*
     &     za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jb,jd)*(za(jb,jg)*(za(jn,jc
     &     )*zb(jb,ju)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg
     &     ,je)) + zb(jg,jb)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,j
     &     u)) + zb(jb,ju)*(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je)))) + za(
     &     jb,jn)*za(jc,jg)*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) - zb(jc,jb)*zb(jg,je)) 
     &     + zb(je,jb)*(zb(jc,ju)*zb(jg,jb) - zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(
     &     za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)
     &     *zb(jg,je) - za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(ju,jn)*zb(jc,j
     &     b)*zb(je,ju)*zb(jg,ju) - za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn
     &     ,jc)*zb(jc,je)*(zb(jc,ju)*zb(jg,jb) + zb(jb,ju)*zb(jg,jc) - zb(jc,jb)*zb
     &     (jg,ju)) + za(jn,jd)*(-(zb(jc,ju)*zb(jd,je)*zb(jg,jb)) + zb(jd,jc)*zb(je
     &     ,ju)*zb(jg,jb) + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))
     &     ) - za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) - za(jb,
     &     jg)*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jb,jn)*za(jc,
     &     jg)*za(jn,jd)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,
     &     jd)*zb(je,jb)*zb(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,
     &     je)*zb(jg,jc)*zb(jn,ju) + za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,
     &     jc)*zb(jn,ju) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,
     &     ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(j
     &     b,jg)*za(jn,jd)*za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(j
     &     n,jd)*za(jn,jg)*zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(
     &     ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(
     &     jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + z
     &     a(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*z
     &     b(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*z
     &     b(je,ju) - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)
     &     *za(jn,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)
     &     *zb(jd,je)*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb)
     &      + za(je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,j
     &     g)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,j
     &     u)*zb(jg,jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd
     &     ,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd
     &     ,jc)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg
     &     ,jc) - za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(
     &     jd,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(
     &     je,ju)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + z
     &     a(jd,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*z
     &     b(je,jb)*zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*z
     &     b(jg,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)
     &     *za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb
     &     )*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*
     &     zb(jc,jb)*zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,j
     &     e)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*
     &     (2*za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(j
     &     b,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,j
     &     b)*zb(je,ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + z
     &     b(jc,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju
     &     ) + za(jb,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,
     &     jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,
     &     jc)*zb(jg,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(j
     &     n,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(j
     &     c,jb)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(j
     &     n,ju) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za
     &     (jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jn,ju))) - (im*imag(anomc7) + real(anomc7
     &     ))*(2*gb**2*za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) 
     &     + 2*gb**2*za(jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + 
     &     gb**2*za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) + 3*gw
     &     **2*za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - gb**2*
     &     za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) - 3*gw**2*za
     &     (jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + gb**2*za(jb,
     &     jg)*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + 3*gw**2*za(jb,jg
     &     )*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + gb**2*za(jb,jg)*za
     &     (jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + 3*gw**2*za(jb,jg)*za(j
     &     n,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - gb**2*za(jb,jg)*za(jc,jd
     &     )*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - 3*gw**2*za(jb,jg)*za(jc,jd)*
     &     za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + gb**2*za(jb,jg)*za(jc,jd)*za(j
     &     n,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) + 3*gw**2*za(jb,jg)*za(jc,jd)*za(jn,
     &     jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - gb**2*za(jb,jg)*za(je,jd)*za(jn,jc)*
     &     zb(je,jb)*zb(je,ju)*zb(jg,jc) - 3*gw**2*za(jb,jg)*za(je,jd)*za(jn,jc)*zb
     &     (je,jb)*zb(je,ju)*zb(jg,jc) + gb**2*za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,
     &     jb)*zb(je,ju)*zb(jg,jc) + 3*gw**2*za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb
     &     )*zb(je,ju)*zb(jg,jc) + gb**2*za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*zb
     &     (jg,jc)*zb(jg,je) + 3*gw**2*za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju)*zb(j
     &     g,jc)*zb(jg,je) + za(jb,jd)*(2*gb**2*za(jb,jn)*za(jc,jg)*(zb(jb,ju)*zb(j
     &     c,je) + zb(jc,ju)*zb(je,jb))*zb(jg,jb) + (gb**2 + 3*gw**2)*za(jb,jg)*za(
     &     jn,jc)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) + 2*gb**2*za(jc,jg)*zb(jg,jb)*(za(j
     &     n,jc)*zb(jc,je)*zb(jc,ju) + za(jn,jd)*(-(zb(jc,ju)*zb(jd,je)) + zb(jd,jc
     &     )*zb(je,ju)) + zb(jc,ju)*(-(za(ju,jn)*zb(je,ju)) + za(jn,jg)*zb(jg,je)))
     &     ) - gb**2*za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - 
     &     3*gw**2*za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) + gb
     &     **2*za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + 3*gw**
     &     2*za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - gb**2*za
     &     (jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju) - 3*gw**2*za(j
     &     b,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju) + za(jb,jn)*(za(
     &     jc,jd)*(2*gb**2*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) + (gb**2 + 3*gw*
     &     *2)*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) - 2*gb**2*za(jc,jg)*zb(jg,j
     &     b)*(za(je,jd)*zb(jc,je)*zb(je,ju) + za(jd,jg)*zb(jc,je)*zb(jg,ju) + za(j
     &     n,jd)*(-(zb(jc,jn)*zb(je,ju)) + zb(jc,je)*zb(jn,ju)))))))/(3._dp*ecossin**2
     &     *zb(jg,jb)*zb(jg,jc))
       end function streal_heavyWWZ_MPMM_M_L2

       function streal_heavyWWZ_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_PMMM_P_L2 =
     &     (propW16*propW34*propZ257*(3*gw**2*im*imag(anomc4)*(-(za(jb,jn)*za(jc,jd
     &     )*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb)) - za(jb,jg)*za(jb,jn)*za(jc,j
     &     d)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb)*zb(jb,j
     &     u)*zb(je,ju)*zb(jg,jb) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je)*zb(je,j
     &     u)*zb(jg,jb) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,j
     &     b) - za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb) - za(jc
     &     ,jg)*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jb
     &     ,jn)*za(je,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,jg)*za(jn
     &     ,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb
     &     ,ju)*zb(jc,je)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd
     &     ,je)*zb(jg,jc) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg
     &     ,jc) + za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(
     &     jb,jn)*za(jc,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(
     &     jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(jn,jd)*za(
     &     ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(ju,jn)*zb(
     &     jb,ju)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(
     &     je,ju)*zb(jg,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(
     &     jg,jc) + za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jc) - z
     &     a(jc,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*z
     &     a(jc,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*z
     &     a(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za(jn,jd)*z
     &     b(jd,jc)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(je,jb)*z
     &     b(je,ju)*zb(jg,jc) + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb(je,ju)*z
     &     b(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jd) -
     &      za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + za(jc,jg)
     &     *za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*za(jc,jd)
     &     *za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jc)
     &     *zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,jb)
     &     *zb(je,ju)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(je,ju)
     &     *zb(jg,je) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je)
     &      + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) + za(jb,j
     &     g)*za(jc,jd)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(je,j
     &     g)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,j
     &     g)*zb(je,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,j
     &     u)*zb(jg,jc)*zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,ju)*zb(jg
     &     ,jc)*zb(jg,je) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg
     &     ,je) + za(jc,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,je) - za(
     &     jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(
     &     jn,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(jn,jd)*za(
     &     jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(jn,jg)*zb(
     &     je,ju)*zb(jg,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,jb)*zb(
     &     jc,je)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(jd,je)*zb(
     &     jg,ju) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - z
     &     a(jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za(jb,jg)*z
     &     a(jn,jd)*za(ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*z
     &     a(ju,jn)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(jc,jg)*za(jd,jg)*z
     &     b(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*z
     &     b(jg,jb)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb(jg,jb)*z
     &     b(jg,ju) + 2*za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*zb(jg,ju)
     &      - za(jb,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - za(jb,j
     &     g)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,jg)*za(jd
     &     ,jg)*za(jn,jc)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(jn
     &     ,jd)*zb(jd,je)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(je
     &     ,jb)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg
     &     ,jc)*zb(jg,ju) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jc)*zb(jg
     &     ,ju) + za(jc,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,ju) + za(
     &     jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(
     &     jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(jc,jd)*za(
     &     jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jd)*(-(za(jb,jg)*(za(jn,jc)
     &     *zb(jb,ju)*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,
     &     je)) + zb(jg,jb)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju
     &     )) + zb(jb,ju)*(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je))))) + za(
     &     jb,jn)*za(jc,jg)*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) + zb(jc,jb)*zb(jg,je)) 
     &     + zb(je,jb)*(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*zb(jg,ju))) + za(jc,jg)*(-(
     &     za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*zb(jg,jb)
     &     *zb(jg,je) + za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) - za(ju,jn)*zb(jc,j
     &     b)*zb(je,ju)*zb(jg,ju) + za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(jn
     &     ,jc)*zb(jc,je)*(zb(jc,ju)*zb(jg,jb) - zb(jb,ju)*zb(jg,jc) + zb(jc,jb)*zb
     &     (jg,ju)) - za(jn,jd)*(zb(jc,ju)*zb(jd,je)*zb(jg,jb) - zb(jd,jc)*zb(je,ju
     &     )*zb(jg,jb) + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))) -
     &      za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) + za(jb,jg)
     &     *za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + za(jb,jn)*za(jc,jg)
     &     *za(jn,jd)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*za(jn,jd)
     &     *zb(je,jb)*zb(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,je)
     &     *zb(jg,jc)*zb(jn,ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)
     &     *zb(jn,ju) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju)
     &      - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) - za(jb,j
     &     g)*za(jn,jd)*za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg)*za(jn,j
     &     d)*za(jn,jg)*zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,jd)*za(ju,
     &     jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,
     &     je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(j
     &     n,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju,jc)*zb(j
     &     c,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(j
     &     e,ju) - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(jd,jg)*za
     &     (jn,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(jc,ju)*zb
     &     (jd,je)*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + 
     &     za(je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*za(jn,jg)*
     &     zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*zb(je,ju)*
     &     zb(jg,jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) - za(jd,jg
     &     )*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)*zb(jd,jc
     &     )*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc
     &     ) - za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,jn)*za(jd,
     &     jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,ju)*zb(je,
     &     ju)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(j
     &     d,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn,jd)*zb(j
     &     e,jb)*zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je,ju)*zb(j
     &     g,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jd,jg)*za
     &     (jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb(jc,jb)*(
     &     zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(jb,jn)*zb(
     &     jc,jb)*zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(zb(jc,je)*
     &     zb(jg,jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(jb,jd)*(2*
     &     za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(-(zb(jb,j
     &     u)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(jc,jb)*
     &     zb(je,ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) + zb(j
     &     c,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) +
     &      za(jb,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*za(jn,jg)
     &     *zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*zb(jg,jc)
     &     *zb(jg,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) - za(jn,j
     &     d)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd)*zb(jc,j
     &     b)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,j
     &     u) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,jd)*za(jn
     &     ,jg)*zb(je,jb)*zb(jg,jc)*zb(jn,ju))) - 3*gw**2*real(anomc4)*(-(za(jb,jn)
     &     *za(jc,jd)*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb)) - za(jb,jg)*za(jb,jn
     &     )*za(jc,jd)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) + za(jb,jg)*za(jn,jd)*za(ju,jb
     &     )*zb(jb,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jn)*za(jc,jg)*za(je,jd)*zb(jc,je
     &     )*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju
     &     )*zb(jg,jb) - za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(jg,jb
     &     ) - za(jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,
     &     jg)*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) - za(jb,je)*za(jb,
     &     jg)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*za(jc,jd)*za(jn,
     &     jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,
     &     ju)*zb(jd,je)*zb(jg,jc) + za(jc,jd)*za(jc,jg)*za(jn,jd)*zb(jc,ju)*zb(jd,
     &     je)*zb(jg,jc) + za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(je,jb)*zb(jg,
     &     jc) - za(jb,jn)*za(jc,jd)*za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(j
     &     b,jg)*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + za(jb,jg)*za(j
     &     n,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(jc,jd)*za(j
     &     u,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jc,jg)*za(je,jd)*za(jn,jc)*zb(j
     &     c,je)*zb(je,ju)*zb(jg,jc) + za(jc,jg)*za(je,jc)*za(jn,jd)*zb(jc,je)*zb(j
     &     e,ju)*zb(jg,jc) + za(jc,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jn)*zb(je,ju)*zb(j
     &     g,jc) - za(jc,jg)*za(jn,jd)*za(ju,jc)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za
     &     (jc,jd)*za(jc,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za
     &     (jc,jd)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(jc,jd)*za(jc,jg)*za
     &     (jn,jd)*zb(jd,jc)*zb(je,ju)*zb(jg,jc) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb
     &     (je,jb)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(je,jb)*zb
     &     (je,ju)*zb(jg,jc) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb
     &     (jg,jd) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,jd) + 
     &     za(jc,jg)*za(jd,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jg,jd) - za(jb,jn)*
     &     za(jc,jd)*za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jg)*za(jc,jd)*
     &     za(jn,jc)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) + za(jb,jn)*za(jc,jg)*za(je,jd)*
     &     zb(jc,jb)*zb(je,ju)*zb(jg,je) - za(jb,jg)*za(je,jd)*za(jn,jc)*zb(jc,jb)*
     &     zb(je,ju)*zb(jg,je) - za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*
     &     zb(jg,je) + za(jb,jg)*za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,je) 
     &     + za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg
     &     )*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jg,je) - za(jb,jg)*za(je,jd
     &     )*za(jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,je) + za(jb,jg)*za(jc,jd)*za(jn,jg
     &     )*zb(jb,ju)*zb(jg,jc)*zb(jg,je) - 2*za(jc,jd)*za(jc,jg)*za(jn,jg)*zb(jc,
     &     ju)*zb(jg,jc)*zb(jg,je) - za(jc,jg)*za(je,jg)*za(jn,jd)*zb(je,ju)*zb(jg,
     &     jc)*zb(jg,je) + za(jc,jg)*za(je,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,
     &     je) - za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(j
     &     b,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jg,jn) + za(jb,jg)*za(j
     &     n,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,jn) - za(jc,jg)*za(jn,jd)*za(j
     &     n,jg)*zb(je,ju)*zb(jg,jc)*zb(jg,jn) - za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(j
     &     c,jb)*zb(jc,je)*zb(jg,ju) + za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jc,jb)*zb(j
     &     d,je)*zb(jg,ju) - za(jb,jg)*za(jb,jn)*za(jc,jd)*zb(jc,jb)*zb(je,jb)*zb(j
     &     g,ju) - za(jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) - za
     &     (jb,jg)*za(jn,jd)*za(ju,jc)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jg)*za
     &     (jc,jd)*za(ju,jn)*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jb,jn)*za(jc,jg)*za
     &     (jd,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jb,jg)*za(jd,jg)*za(jn,jc)*zb
     &     (jc,je)*zb(jg,jb)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jd)*zb(jd,je)*zb
     &     (jg,jb)*zb(jg,ju) + 2*za(jb,jg)*za(jb,jn)*za(jd,jg)*zb(je,jb)*zb(jg,jb)*
     &     zb(jg,ju) - za(jb,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) 
     &     - za(jb,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) - 2*za(jc,
     &     jg)*za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jd,
     &     jg)*za(jn,jd)*zb(jd,je)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jc,jg)*za(jd,
     &     jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,
     &     jb)*zb(jg,jc)*zb(jg,ju) + za(jc,jg)*za(jn,jd)*za(ju,jg)*zb(je,ju)*zb(jg,
     &     jc)*zb(jg,ju) + za(jc,jg)*za(jd,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc)*zb(jg,
     &     ju) + za(jb,jn)*za(jc,jg)*za(jd,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(j
     &     b,jg)*za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jg)*za(j
     &     c,jd)*za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(jb,jd)*(-(za(jb,jg)*(
     &     za(jn,jc)*zb(jb,ju)*(zb(jc,je)*zb(jg,jb) + zb(je,jb)*zb(jg,jc) + zb(jc,j
     &     b)*zb(jg,je)) + zb(jg,jb)*(za(jn,jd)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)
     &     *zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*zb(je,ju)) + 2*za(jn,jg)*zb(jg,je))
     &     ))) + za(jb,jn)*za(jc,jg)*(zb(jb,ju)*(zb(jc,je)*zb(jg,jb) + zb(jc,jb)*zb
     &     (jg,je)) + zb(je,jb)*(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*zb(jg,ju))) + za(j
     &     c,jg)*(-(za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb)) + za(jn,jg)*zb(jc,ju)*
     &     zb(jg,jb)*zb(jg,je) + za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) - za(ju,jn
     &     )*zb(jc,jb)*zb(je,ju)*zb(jg,ju) + za(jn,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju
     &     ) + za(jn,jc)*zb(jc,je)*(zb(jc,ju)*zb(jg,jb) - zb(jb,ju)*zb(jg,jc) + zb(
     &     jc,jb)*zb(jg,ju)) - za(jn,jd)*(zb(jc,ju)*zb(jd,je)*zb(jg,jb) - zb(jd,jc)
     &     *zb(je,ju)*zb(jg,jb) + zb(jc,jb)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,
     &     ju))))) - za(jb,jg)*za(jb,jn)*za(jn,jd)*zb(je,ju)*zb(jg,jb)*zb(jn,jb) + 
     &     za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + za(jb,jn)*
     &     za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jb,jg)*za(jb,jn)*
     &     za(jn,jd)*zb(je,jb)*zb(jg,jb)*zb(jn,ju) - za(jc,jg)*za(jn,jc)*za(jn,jd)*
     &     zb(jc,je)*zb(jg,jc)*zb(jn,ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*
     &     zb(jg,jc)*zb(jn,ju) + za(jb,jn)*za(jc,jg)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*
     &     zb(jn,ju) - za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jg,je)*zb(jn,ju) 
     &     - za(jb,jg)*za(jn,jd)*za(jn,jg)*zb(jg,jb)*zb(jg,je)*zb(jn,ju) + za(jc,jg
     &     )*za(jn,jd)*za(jn,jg)*zb(jg,jc)*zb(jg,je)*zb(jn,ju) + za(jb,jc)*(za(jn,j
     &     d)*za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,jc)*zb(jc,j
     &     b)*zb(jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(je,j
     &     u) + za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(jn,jd)*za(ju
     &     ,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je
     &     ,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju) + za(
     &     jd,jg)*za(jn,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za(jn,jd)*zb(
     &     jc,ju)*zb(jd,je)*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb(je,jb)*zb(
     &     jg,jb) + za(je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - za(je,jd)*z
     &     a(jn,jg)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*zb(jc,jn)*z
     &     b(je,ju)*zb(jg,jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) -
     &      za(jd,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jd)
     &     *zb(jd,jc)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju)*zb(jc,je)
     &     *zb(jg,jc) - za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + za(jb,j
     &     n)*za(jd,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,jg)*zb(jb,j
     &     u)*zb(je,ju)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,j
     &     c) + za(jd,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(je,jg)*za(jn
     &     ,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(je,jb)*zb(je
     &     ,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(
     &     jd,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(za(jn,jd)*zb
     &     (jc,jb)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,ju)*(-2*za(j
     &     b,jn)*zb(jc,jb)*zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + za(jn,jg)*(z
     &     b(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je)))) + za(j
     &     b,jd)*(2*za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*zb(jc,jb)*(
     &     -(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(za(ju,jn)*z
     &     b(jc,jb)*zb(je,ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,j
     &     c) + zb(jc,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb
     &     (jg,ju) + za(jb,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) - za(jd,jg)*
     &     za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)*zb(je,jb)*
     &     zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)*zb(jn,jb) 
     &     - za(jn,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,jc)*za(jn,jd
     &     )*zb(jc,jb)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb
     &     )*zb(jn,ju) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,ju) + za(jn,
     &     jd)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jn,ju))) + (im*imag(anomc7) - real(
     &     anomc7))*(2*gb**2*za(jb,je)*za(jc,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(j
     &     g,jb) + 2*gb**2*za(jc,jg)*za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,
     &     jb) + gb**2*za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) 
     &     + 3*gw**2*za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - 
     &     gb**2*za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) - 3*gw
     &     **2*za(jb,jg)*za(jc,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc) + gb**2*
     &     za(jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + 3*gw**2*za
     &     (jb,jg)*za(jc,jd)*za(jn,jc)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) + gb**2*za(jb,
     &     jg)*za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + 3*gw**2*za(jb,jg
     &     )*za(jn,jd)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - gb**2*za(jb,jg)*za
     &     (jc,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - 3*gw**2*za(jb,jg)*za(j
     &     c,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + gb**2*za(jb,jg)*za(jc,jd
     &     )*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) + 3*gw**2*za(jb,jg)*za(jc,jd)*
     &     za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - gb**2*za(jb,jg)*za(je,jd)*za(j
     &     n,jc)*zb(je,jb)*zb(je,ju)*zb(jg,jc) - 3*gw**2*za(jb,jg)*za(je,jd)*za(jn,
     &     jc)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + gb**2*za(jb,jg)*za(je,jc)*za(jn,jd)*
     &     zb(je,jb)*zb(je,ju)*zb(jg,jc) + 3*gw**2*za(jb,jg)*za(je,jc)*za(jn,jd)*zb
     &     (je,jb)*zb(je,ju)*zb(jg,jc) + gb**2*za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,
     &     ju)*zb(jg,jc)*zb(jg,je) + 3*gw**2*za(jb,jg)*za(jc,jd)*za(jn,jg)*zb(jb,ju
     &     )*zb(jg,jc)*zb(jg,je) + za(jb,jd)*(2*gb**2*za(jb,jn)*za(jc,jg)*(zb(jb,ju
     &     )*zb(jc,je) + zb(jc,ju)*zb(je,jb))*zb(jg,jb) + (gb**2 + 3*gw**2)*za(jb,j
     &     g)*za(jn,jc)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) + 2*gb**2*za(jc,jg)*zb(jg,jb)
     &     *(za(jn,jc)*zb(jc,je)*zb(jc,ju) + za(jn,jd)*(-(zb(jc,ju)*zb(jd,je)) + zb
     &     (jd,jc)*zb(je,ju)) + zb(jc,ju)*(-(za(ju,jn)*zb(je,ju)) + za(jn,jg)*zb(jg
     &     ,je)))) - gb**2*za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,
     &     ju) - 3*gw**2*za(jb,jg)*za(jd,jg)*za(jn,jc)*zb(je,jb)*zb(jg,jc)*zb(jg,ju
     &     ) + gb**2*za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) + 
     &     3*gw**2*za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - gb
     &     **2*za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju) - 3*gw**
     &     2*za(jb,jg)*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jg,jc)*zb(jn,ju) + za(jb,jn
     &     )*(za(jc,jd)*(2*gb**2*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) + (gb**2 +
     &      3*gw**2)*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) - 2*gb**2*za(jc,jg)*z
     &     b(jg,jb)*(za(je,jd)*zb(jc,je)*zb(je,ju) + za(jd,jg)*zb(jc,je)*zb(jg,ju) 
     &     + za(jn,jd)*(-(zb(jc,jn)*zb(je,ju)) + zb(jc,je)*zb(jn,ju)))))))/(3._dp*ecos
     &     sin**2*za(jb,jg)*za(jc,jg))
       end function streal_heavyWWZ_PMMM_P_L2

       function streal_heavyWWZ_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyWWZ_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propW16, propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW16  = 1._dp / (s(ju,jd) - wmass**2)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyWWZ_PMMM_M_L2 =
     &     (gw**2*propW16*propW34*propZ257*(im*imag(anomc4) - real(anomc4))*zb(jc,j
     &     b)*(za(jn,jd)*za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(je,jd)*za(jn,
     &     jc)*zb(jc,jb)*zb(jc,je)*zb(je,ju) + za(je,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,
     &     je)*zb(je,ju) + za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,jn)*zb(je,ju) - za(j
     &     n,jd)*za(ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jn)*za(je,jd)*zb(j
     &     c,jb)*zb(je,jb)*zb(je,ju) - za(jb,je)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(j
     &     e,ju) + za(jd,jg)*za(jn,jc)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jd,jg)*za
     &     (jn,jd)*zb(jc,ju)*zb(jd,je)*zb(jg,jb) + za(jb,jn)*za(jd,jg)*zb(jc,ju)*zb
     &     (je,jb)*zb(jg,jb) + za(je,jg)*za(jn,jd)*zb(jc,je)*zb(je,ju)*zb(jg,jb) - 
     &     za(je,jd)*za(jn,jg)*zb(jc,je)*zb(je,ju)*zb(jg,jb) + za(jn,jd)*za(jn,jg)*
     &     zb(jc,jn)*zb(je,ju)*zb(jg,jb) - za(jn,jd)*za(ju,jg)*zb(jc,ju)*zb(je,ju)*
     &     zb(jg,jb) - za(jd,jg)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jd,jg
     &     )*za(jn,jd)*zb(jd,jc)*zb(je,ju)*zb(jg,jb) + za(jd,jg)*za(jn,jc)*zb(jb,ju
     &     )*zb(jc,je)*zb(jg,jc) - za(jd,jg)*za(jn,jd)*zb(jb,ju)*zb(jd,je)*zb(jg,jc
     &     ) + za(jb,jn)*za(jd,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jn,jd)*za(ju,
     &     jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) - za(jd,jg)*za(ju,jn)*zb(jb,ju)*zb(je,
     &     ju)*zb(jg,jc) + za(jd,jg)*za(jn,jd)*zb(jd,jb)*zb(je,ju)*zb(jg,jc) - za(j
     &     e,jg)*za(jn,jd)*zb(je,jb)*zb(je,ju)*zb(jg,jc) + za(je,jd)*za(jn,jg)*zb(j
     &     e,jb)*zb(je,ju)*zb(jg,jc) + za(jd,jg)*za(jn,jg)*zb(jc,ju)*zb(jg,jb)*zb(j
     &     g,je) + za(jd,jg)*za(jn,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jd)*(z
     &     a(jn,jd)*zb(jc,jb)*(zb(jc,ju)*zb(jd,je) - zb(jd,jc)*zb(je,ju)) + zb(jc,j
     &     u)*(-2*za(jb,jn)*zb(jc,jb)*zb(je,jb) + za(ju,jn)*zb(jc,jb)*zb(je,ju) + z
     &     a(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,jb)*zb(jg,jc) - zb(jc,jb)*zb(jg,je
     &     )))) + za(jb,jd)*(2*za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je) + za(jn,jd)*
     &     zb(jc,jb)*(-(zb(jb,ju)*zb(jd,je)) + zb(jd,jb)*zb(je,ju)) + zb(jb,ju)*(-(
     &     za(ju,jn)*zb(jc,jb)*zb(je,ju)) + za(jn,jg)*(zb(jc,je)*zb(jg,jb) - zb(je,
     &     jb)*zb(jg,jc) + zb(jc,jb)*zb(jg,je)))) - za(jd,jg)*za(jn,jc)*zb(jc,jb)*z
     &     b(jc,je)*zb(jg,ju) + za(jb,jn)*za(jd,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju) -
     &      za(jd,jg)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) + za(jd,jg)*za(jn,jg)
     &     *zb(je,jb)*zb(jg,jc)*zb(jg,ju) - za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,ju)
     &     *zb(jn,jb) - za(jn,jd)*za(jn,jg)*zb(je,ju)*zb(jg,jc)*zb(jn,jb) - za(jn,j
     &     c)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(jn,ju) + za(jb,jn)*za(jn,jd)*zb(jc,j
     &     b)*zb(je,jb)*zb(jn,ju) - za(jn,jd)*za(jn,jg)*zb(jc,je)*zb(jg,jb)*zb(jn,j
     &     u) + za(jn,jd)*za(jn,jg)*zb(je,jb)*zb(jg,jc)*zb(jn,ju)))/(ecossin**2*zb(
     &     jg,jb)*zb(jg,jc))
       end function streal_heavyWWZ_PMMM_M_L2

       function streal_lightGR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_MPMM_P_L2 =
     &     (4*gb**2*propW34*(im*imag(anomc4) + real(anomc4))*za(jn,jd)*(za(ju,jd)*(
     &     -(za(jb,jd)*zb(jb,ju)*zb(jd,je)) + za(jc,jd)*zb(jc,ju)*zb(jd,je) + (za(j
     &     b,jn)*zb(jb,ju) + za(jn,jc)*zb(jc,ju))*zb(je,jn)) + za(jd,jg)*(-(za(jb,j
     &     d)*zb(jd,je)*zb(jg,jb)) + za(jb,jn)*zb(je,jn)*zb(jg,jb) + (za(jc,jd)*zb(
     &     jd,je) + za(jn,jc)*zb(je,jn))*zb(jg,jc)) + (za(jn,jd)*za(ju,jg)*zb(je,jn
     &     )*(za(jb,jc)*(zb(jc,ju)*zb(jg,jb) + zb(jb,ju)*zb(jg,jc)) + (-(za(ju,jb)*
     &     zb(jb,ju)) + za(ju,jc)*zb(jc,ju))*zb(jg,ju)))/(s(jb,jc) + s(jb,ju) + s(j
     &     c,ju))))/(3._dp*ecossin**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(jd,jg)*za(ju
     &     ,jg)*zb(jc,jb))
       end function streal_lightGR_MPMM_P_L2

       function streal_lightGR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_MPMM_M_L2 =
     &     (4*gb**2*propW34*(im*imag(anomc4) + real(anomc4))*((za(ju,jb)*zb(jb,ju) 
     &     - za(ju,jc)*zb(jc,ju))*zb(je,ju)*(za(jn,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,
     &     ju)) + za(jb,jc)*(za(jn,jd)*(-(zb(jc,ju)*zb(jd,ju)*zb(je,jb)) + (zb(jb,j
     &     u)*((s(jd,je) + s(jd,jn) + s(je,jn))*zb(jc,je)*zb(jd,ju) + 2*za(jd,jg)*z
     &     b(jc,ju)*zb(jd,je)*zb(jg,jd) - 2*za(jn,jg)*zb(jc,ju)*zb(je,jn)*zb(jg,jd)
     &     ))/(s(jd,je) + s(jd,jn) + s(je,jn))) + za(jn,jg)*(zb(jb,ju)*zb(jc,je) - 
     &     zb(jc,ju)*zb(je,jb))*zb(jg,ju))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,ju) + 
     &     s(jc,ju))*zb(jc,jb)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightGR_MPMM_M_L2

       function streal_lightGR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_PMMM_P_L2 =
     &     (4*gb**2*propW34*(im*imag(anomc4) - real(anomc4))*za(jn,jd)*(za(ju,jd)*(
     &     -(za(jb,jd)*zb(jb,ju)*zb(jd,je)) + za(jc,jd)*zb(jc,ju)*zb(jd,je) + (za(j
     &     b,jn)*zb(jb,ju) + za(jn,jc)*zb(jc,ju))*zb(je,jn)) + za(jd,jg)*(-(za(jb,j
     &     d)*zb(jd,je)*zb(jg,jb)) + za(jb,jn)*zb(je,jn)*zb(jg,jb) + (za(jc,jd)*zb(
     &     jd,je) + za(jn,jc)*zb(je,jn))*zb(jg,jc)) + (za(jn,jd)*za(ju,jg)*zb(je,jn
     &     )*(za(jb,jc)*(zb(jc,ju)*zb(jg,jb) + zb(jb,ju)*zb(jg,jc)) + (-(za(ju,jb)*
     &     zb(jb,ju)) + za(ju,jc)*zb(jc,ju))*zb(jg,ju)))/(s(jb,jc) + s(jb,ju) + s(j
     &     c,ju))))/(3._dp*ecossin**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(jb,jc)*za(jd
     &     ,jg)*za(ju,jg))
       end function streal_lightGR_PMMM_P_L2

       function streal_lightGR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGR_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGR_PMMM_M_L2 =
     &     (4*gb**2*propW34*(im*imag(anomc4) - real(anomc4))*((za(ju,jb)*zb(jb,ju) 
     &     - za(ju,jc)*zb(jc,ju))*zb(je,ju)*(za(jn,jd)*zb(jd,ju) + za(jn,jg)*zb(jg,
     &     ju)) + za(jb,jc)*(za(jn,jd)*(-(zb(jc,ju)*zb(jd,ju)*zb(je,jb)) + (zb(jb,j
     &     u)*((s(jd,je) + s(jd,jn) + s(je,jn))*zb(jc,je)*zb(jd,ju) + 2*za(jd,jg)*z
     &     b(jc,ju)*zb(jd,je)*zb(jg,jd) - 2*za(jn,jg)*zb(jc,ju)*zb(je,jn)*zb(jg,jd)
     &     ))/(s(jd,je) + s(jd,jn) + s(je,jn))) + za(jn,jg)*(zb(jb,ju)*zb(jc,je) - 
     &     zb(jc,ju)*zb(je,jb))*zb(jg,ju))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,ju) + 
     &     s(jc,ju))*za(jb,jc)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightGR_PMMM_M_L2

       function streal_lightGL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_MPMM_P_L2 =
     &     (2*gb**2*propW34*(im*imag(anomc4) + real(anomc4))*(za(ju,jd)*(za(jb,jd)*
     &     (za(jn,jc)*zb(jc,jb) + za(jn,jd)*zb(jd,jb)) - za(jc,jd)*(za(jb,jn)*zb(jc
     &     ,jb) + za(jn,jd)*zb(jd,jc)))*zb(je,ju) + za(jd,jg)*(za(jb,jd)*(za(jn,jc)
     &     *zb(jc,jb) + za(jn,jd)*zb(jd,jb)) - za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(
     &     jn,jd)*zb(jd,jc)))*zb(jg,je) + (2*za(jb,jd)*za(jc,jd)*za(ju,jg)*zb(jc,jb
     &     )*zb(je,ju)*(za(jn,je)*zb(jg,je) - za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(j
     &     e,ju) + s(jn,ju))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,jd) + s(jc,jd))*za(j
     &     d,jg)*za(ju,jg)*zb(jc,jb))
       end function streal_lightGL_MPMM_P_L2

       function streal_lightGL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_MPMM_M_L2 =
     &     (2*gb**2*propW34*(-(im*imag(anomc4)) - real(anomc4))*zb(je,ju)*(za(ju,jn
     &     )*(za(jb,jd)*zb(jb,ju)*zb(jd,ju) - za(jc,jd)*zb(jc,ju)*zb(jd,ju) + (za(j
     &     b,jg)*zb(jb,ju) - za(jc,jg)*zb(jc,ju))*zb(jg,ju)) + za(jn,je)*(za(jb,jd)
     &     *(zb(jd,ju)*zb(je,jb) + ((za(jc,jg)*zb(jc,jb) + za(jd,jg)*zb(jd,jb))*zb(
     &     je,ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) + za(jc,jd)*(zb(jc,j
     &     e)*zb(jd,ju) + ((za(jb,jg)*zb(jc,jb) - za(jd,jg)*zb(jd,jc))*zb(je,ju)*zb
     &     (jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) + (za(jc,jg)*zb(jc,je) + za(j
     &     b,jg)*zb(je,jb))*zb(jg,ju))))/(3._dp*ecossin**2*(s(je,jn) + s(je,ju) + s(jn
     &     ,ju))*zb(jc,jb)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightGL_MPMM_M_L2

       function streal_lightGL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_PMMM_P_L2 =
     &     (2*gb**2*propW34*(im*imag(anomc4)*(za(ju,jd)*(za(jb,jd)*(za(jn,jc)*zb(jc
     &     ,jb) + za(jn,jd)*zb(jd,jb)) - za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(jn,jd)
     &     *zb(jd,jc)))*zb(je,ju) + za(jd,jg)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(
     &     jn,jd)*zb(jd,jb)) - za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)
     &     ))*zb(jg,je) + (2*za(jb,jd)*za(jc,jd)*za(ju,jg)*zb(jc,jb)*zb(je,ju)*(za(
     &     jn,je)*zb(jg,je) - za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(je,ju) + s(jn,ju)
     &     )) + real(anomc4)*(za(ju,jd)*(-(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(jn,j
     &     d)*zb(jd,jb))) + za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)))*
     &     zb(je,ju) + za(jd,jg)*(-(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(jn,jd)*zb(j
     &     d,jb))) + za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)))*zb(jg,j
     &     e) + (2*za(jb,jd)*za(jc,jd)*za(ju,jg)*zb(jc,jb)*zb(je,ju)*(-(za(jn,je)*z
     &     b(jg,je)) + za(ju,jn)*zb(jg,ju)))/(s(je,jn) + s(je,ju) + s(jn,ju)))))/(3
     &     ._dp*ecossin**2*(s(jb,jc) + s(jb,jd) + s(jc,jd))*za(jb,jc)*za(jd,jg)*za(ju,
     &     jg))
       end function streal_lightGL_PMMM_P_L2

       function streal_lightGL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightGL_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_lightGL_PMMM_M_L2 =
     &     (2*gb**2*propW34*(-(im*imag(anomc4)) + real(anomc4))*zb(je,ju)*(za(ju,jn
     &     )*(za(jb,jd)*zb(jb,ju)*zb(jd,ju) - za(jc,jd)*zb(jc,ju)*zb(jd,ju) + (za(j
     &     b,jg)*zb(jb,ju) - za(jc,jg)*zb(jc,ju))*zb(jg,ju)) + za(jn,je)*(za(jb,jd)
     &     *(zb(jd,ju)*zb(je,jb) + ((za(jc,jg)*zb(jc,jb) + za(jd,jg)*zb(jd,jb))*zb(
     &     je,ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) + za(jc,jd)*(zb(jc,j
     &     e)*zb(jd,ju) + ((za(jb,jg)*zb(jc,jb) - za(jd,jg)*zb(jd,jc))*zb(je,ju)*zb
     &     (jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) + (za(jc,jg)*zb(jc,je) + za(j
     &     b,jg)*zb(je,jb))*zb(jg,ju))))/(3._dp*ecossin**2*(s(je,jn) + s(je,ju) + s(jn
     &     ,ju))*za(jb,jc)*zb(jg,jd)*zb(jg,ju))
       end function streal_lightGL_PMMM_M_L2

       function streal_heavyGR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_MPMM_P_L2 =
     &     (4*gb**2*propW34*(im*imag(anomc4) + real(anomc4))*za(jb,jc)*za(jn,jd)*(z
     &     a(jb,jc)**2*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) - (za(jb,jg)*za(
     &     ju,jc)*zb(je,ju) + za(jc,jg)*(za(ju,jb)*zb(je,ju) + 2*za(jb,jg)*zb(jg,je
     &     )))*zb(jg,ju) + za(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb(je,ju) - za(ju,jc)*zb(
     &     jc,ju)*zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jg,je) - za(jc,jg)*zb(jc,ju)*z
     &     b(jg,je) - za(jc,jg)*zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(je,jb)*zb(jg,ju)
     &     )))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je) + s(jd,jn)
     &      + s(je,jn))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGR_MPMM_P_L2

       function streal_heavyGR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_MPMM_M_L2 =
     &     (4*gb**2*propW34*za(jn,jd)*(2*(im*imag(anomc7) + real(anomc7))*(-(za(jb,
     &     jg)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc)) + za(jb,jc)*(za(jc,jg)*zb(j
     &     c,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) + z
     &     a(jc,jg)*(za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*(zb(jc,ju)
     &     *zb(jg,jb) - zb(jb,ju)*zb(jg,jc))*zb(jg,je))) + 3*im*imag(anomc4)*(za(jb
     &     ,jc)**2*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jb,jg
     &     )**2*zb(jg,jb)*(zb(jb,ju)*zb(jg,je) - zb(je,jb)*zb(jg,ju)) + za(jb,jc)*(
     &     za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc,jb)*zb(jc,ju)*
     &     zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jb) + 2*za(jc,jg)*zb(jc,
     &     je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) - za(j
     &     u,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jb,ju)*zb(jc,je)*zb(j
     &     g,jc) - 2*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jc,jg)*zb(jc,ju)*
     &     zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg
     &     )*zb(jb,ju)*zb(jc,jb)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je
     &     ) - za(jc,jg)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(jc,jb)*zb(je,
     &     jb)*zb(jg,ju)) + za(jb,jg)*(za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(jg,jb) + za
     &     (jc,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) - za(jc,jg)*zb(jb,ju)*zb(jg,jc)*zb
     &     (jg,je) + za(jc,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) - za(ju,jg)*zb(je,ju)*
     &     zb(jg,jb)*zb(jg,ju) + za(jc,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - 2*za(jc,
     &     jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(ju,jc)*zb(je,ju)*(zb(jb,ju)*zb(jg
     &     ,jc) + zb(jc,jb)*zb(jg,ju))) - za(jc,jg)*(za(ju,jb)*zb(je,ju)*(-(zb(jc,j
     &     u)*zb(jg,jb)) + zb(jc,jb)*zb(jg,ju)) + zb(jg,jc)*(za(ju,jc)*zb(jc,ju)*zb
     &     (je,ju) - za(ju,jg)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*(zb(jc,ju)*zb(jg,je)
     &      + zb(jc,je)*zb(jg,ju))))) + 3*real(anomc4)*(za(jb,jc)**2*zb(jc,jb)*(zb(
     &     jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jb,jg)**2*zb(jg,jb)*(zb(jb,
     &     ju)*zb(jg,je) - zb(je,jb)*zb(jg,ju)) + za(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb
     &     (jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jg)*
     &     zb(jb,ju)*zb(jc,je)*zb(jg,jb) + 2*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb
     &     ) - za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*zb(je,
     &     ju)*zb(jg,jb) + za(jc,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - 2*za(jb,jg)*zb
     &     (jb,ju)*zb(je,jb)*zb(jg,jc) - za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - 
     &     za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*
     &     zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) - za(jc,jg)*zb(jc,jb
     &     )*zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju)) + za(jb
     &     ,jg)*(za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jc,ju)*zb(j
     &     g,jb)*zb(jg,je) - za(jc,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jg)*zb
     &     (jc,je)*zb(jg,jb)*zb(jg,ju) - za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) + 
     &     za(jc,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - 2*za(jc,jg)*zb(jc,jb)*zb(jg,je
     &     )*zb(jg,ju) - za(ju,jc)*zb(je,ju)*(zb(jb,ju)*zb(jg,jc) + zb(jc,jb)*zb(jg
     &     ,ju))) - za(jc,jg)*(za(ju,jb)*zb(je,ju)*(-(zb(jc,ju)*zb(jg,jb)) + zb(jc,
     &     jb)*zb(jg,ju)) + zb(jg,jc)*(za(ju,jc)*zb(jc,ju)*zb(je,ju) - za(ju,jg)*zb
     &     (je,ju)*zb(jg,ju) + za(jc,jg)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju)
     &     ))))))/(9._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je) + s(jd,
     &     jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyGR_MPMM_M_L2

       function streal_heavyGR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_PMMM_P_L2 =
     &     (4*gb**2*propW34*za(jn,jd)*(-2*(im*imag(anomc7) - real(anomc7))*(-(za(jb
     &     ,jg)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg,jc)) + za(jb,jc)*(za(jc,jg)*zb(
     &     jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) + 
     &     za(jc,jg)*(za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*(zb(jc,ju
     &     )*zb(jg,jb) - zb(jb,ju)*zb(jg,jc))*zb(jg,je))) + 3*im*imag(anomc4)*(za(j
     &     b,jc)**2*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jb,j
     &     g)**2*zb(jg,jb)*(zb(jb,ju)*zb(jg,je) - zb(je,jb)*zb(jg,ju)) + za(jb,jc)*
     &     (za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc,jb)*zb(jc,ju)
     &     *zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jb) - za(jb,jg)*zb(jc,j
     &     u)*zb(je,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jc
     &     ,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg
     &     ,jc) - za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*zb(jb,ju)*zb(
     &     jc,jb)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) - za(jc,jg)*z
     &     b(jc,jb)*zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju)) 
     &     + za(jb,jg)*(za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*zb(jc,j
     &     u)*zb(jg,jb)*zb(jg,je) + za(jc,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc
     &     ,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) - za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg
     &     ,ju) + za(jc,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - 2*za(jc,jg)*zb(jc,jb)*z
     &     b(jg,je)*zb(jg,ju) + za(ju,jc)*zb(je,ju)*(zb(jb,ju)*zb(jg,jc) - zb(jc,jb
     &     )*zb(jg,ju))) - za(jc,jg)*(za(ju,jb)*zb(je,ju)*(zb(jc,ju)*zb(jg,jb) + zb
     &     (jc,jb)*zb(jg,ju)) + zb(jg,jc)*(za(ju,jc)*zb(jc,ju)*zb(je,ju) - za(ju,jg
     &     )*zb(je,ju)*zb(jg,ju) + za(jc,jg)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg
     &     ,ju))))) - 3*real(anomc4)*(za(jb,jc)**2*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) -
     &      zb(jc,ju)*zb(je,jb)) + za(jb,jg)**2*zb(jg,jb)*(zb(jb,ju)*zb(jg,je) - zb
     &     (je,jb)*zb(jg,ju)) + za(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) 
     &     - za(ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je
     &     )*zb(jg,jb) - za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,
     &     ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(j
     &     c,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(j
     &     g,jc) + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb
     &     (jc,ju)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) - za(jb,jg)*
     &     zb(jc,jb)*zb(je,jb)*zb(jg,ju)) + za(jb,jg)*(za(ju,jb)*zb(jb,ju)*zb(je,ju
     &     )*zb(jg,jb) - za(jc,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jc,jg)*zb(jb,
     &     ju)*zb(jg,jc)*zb(jg,je) + za(jc,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) - za(j
     &     u,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) + za(jc,jg)*zb(je,jb)*zb(jg,jc)*zb(j
     &     g,ju) - 2*za(jc,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(ju,jc)*zb(je,ju)*
     &     (zb(jb,ju)*zb(jg,jc) - zb(jc,jb)*zb(jg,ju))) - za(jc,jg)*(za(ju,jb)*zb(j
     &     e,ju)*(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*zb(jg,ju)) + zb(jg,jc)*(za(ju,jc)
     &     *zb(jc,ju)*zb(je,ju) - za(ju,jg)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*(zb(jc,
     &     ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju)))))))/(9._dp*ecossin**2*(s(jb,jc) + s(j
     &     b,jg) + s(jc,jg))*(s(jd,je) + s(jd,jn) + s(je,jn))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGR_PMMM_P_L2

       function streal_heavyGR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGR_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGR_PMMM_M_L2 =
     &     (4*gb**2*propW34*(im*imag(anomc4) - real(anomc4))*za(jn,jd)*zb(jc,jb)*(z
     &     a(jb,jc)*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(ju,j
     &     b)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,j
     &     u) + za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*zb(jc,ju)*zb(je
     &     ,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(
     &     jb,ju)*zb(jc,je)*zb(jg,jc) - za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - z
     &     a(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*z
     &     b(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je)))/(3._dp*ecossin**2*(s(j
     &     b,jc) + s(jb,jg) + s(jc,jg))*(s(jd,je) + s(jd,jn) + s(je,jn))*zb(jg,jb)*
     &     zb(jg,jc))
       end function streal_heavyGR_PMMM_M_L2

       function streal_heavyGL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_MPMM_P_L2 =
     &     (2*gb**2*propW34*(im*imag(anomc4) + real(anomc4))*za(jb,jc)*zb(je,ju)*(-
     &     (za(jc,jd)*((-(za(jn,je)*za(ju,jb)) + za(jb,je)*za(ju,jn))*zb(je,ju) + z
     &     a(jb,jn)*(za(jn,je)*zb(je,jn) + za(ju,jn)*zb(jn,ju)))) + za(jb,jd)*(2*za
     &     (jc,jd)*(za(jn,je)*zb(jd,je) - za(ju,jn)*zb(jd,ju)) + (za(jn,je)*za(ju,j
     &     c) + za(je,jc)*za(ju,jn))*zb(je,ju) + za(jn,jc)*(za(jn,je)*zb(je,jn) + z
     &     a(ju,jn)*zb(jn,ju)))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(
     &     s(je,jn) + s(je,ju) + s(jn,ju))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGL_MPMM_P_L2

       function streal_heavyGL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_MPMM_M_L2 =
     &     (2*gb**2*propW34*zb(je,ju)*(2*(im*imag(anomc7) + real(anomc7))*(za(jb,jg
     &     )*za(jc,jd)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb))*zb(jg,jc) + za(j
     &     b,jd)*zb(jg,jb)*(za(jd,jg)*(-(za(jn,je)*zb(jd,je)) + za(ju,jn)*zb(jd,ju)
     &     ) + za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb)) + za(jn,je)*za
     &     (jn,jg)*zb(je,jn) + za(jn,je)*za(ju,jg)*zb(je,ju) + za(je,jg)*za(ju,jn)*
     &     zb(je,ju) + za(jn,jg)*za(ju,jn)*zb(jn,ju))) - 3*im*imag(anomc4)*(-(za(jn
     &     ,je)*za(ju,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju)) + za(jb,je)*za(ju,jd)*za(j
     &     u,jn)*zb(jb,ju)*zb(je,ju) + za(je,jd)*za(jn,je)*za(ju,jb)*zb(je,jb)*zb(j
     &     e,ju) - za(jb,je)*za(je,jd)*za(ju,jn)*zb(je,jb)*zb(je,ju) - 2*za(jc,jd)*
     &     za(jd,jg)*za(jn,je)*zb(jd,je)*zb(jg,jc) + 2*za(jc,jd)*za(jd,jg)*za(ju,jn
     &     )*zb(jd,ju)*zb(jg,jc) - za(jd,jg)*za(jn,jc)*za(jn,je)*zb(je,jn)*zb(jg,jc
     &     ) + za(jc,jd)*za(jn,je)*za(jn,jg)*zb(je,jn)*zb(jg,jc) - za(jd,jg)*za(jn,
     &     je)*za(ju,jc)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jn,je)*za(ju,jg)*zb(je,
     &     ju)*zb(jg,jc) - za(jd,jg)*za(je,jc)*za(ju,jn)*zb(je,ju)*zb(jg,jc) + za(j
     &     c,jd)*za(je,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc) + za(jn,jd)*za(jn,je)*za(j
     &     u,jb)*zb(je,ju)*zb(jn,jb) - za(jb,je)*za(jn,jd)*za(ju,jn)*zb(je,ju)*zb(j
     &     n,jb) - za(jd,jg)*za(jn,jc)*za(ju,jn)*zb(jg,jc)*zb(jn,ju) + za(jc,jd)*za
     &     (jn,jg)*za(ju,jn)*zb(jg,jc)*zb(jn,ju) + za(jb,jn)*(za(ju,jd)*zb(jb,ju) -
     &      za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb))*(za(jn,je)*zb(je,jn) + za(ju
     &     ,jn)*zb(jn,ju)) + za(jb,jd)*(za(jn,je)**2*zb(je,jb)*zb(je,jn) + za(ju,jn
     &     )*(2*za(ju,jd)*zb(jb,ju)*zb(jd,ju) + za(ju,je)*zb(jb,ju)*zb(je,ju) + za(
     &     je,jd)*(-2*zb(jd,ju)*zb(je,jb) + zb(jd,jb)*zb(je,ju)) - 2*za(jn,jd)*zb(j
     &     d,ju)*zb(jn,jb) + za(ju,jn)*zb(jb,ju)*zb(jn,ju) + za(jn,jd)*zb(jd,jb)*zb
     &     (jn,ju)) + za(jn,je)*(2*za(je,jd)*zb(jd,je)*zb(je,jb) + za(ju,jn)*zb(jb,
     &     ju)*zb(je,jn) + za(jn,jd)*zb(jd,jb)*zb(je,jn) + za(ju,je)*zb(je,jb)*zb(j
     &     e,ju) + za(ju,jd)*(-2*zb(jb,ju)*zb(jd,je) + zb(jd,jb)*zb(je,ju)) + 2*za(
     &     jn,jd)*zb(jd,je)*zb(jn,jb) + za(ju,jn)*zb(je,jb)*zb(jn,ju)))) - 3*real(a
     &     nomc4)*(-(za(jn,je)*za(ju,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju)) + za(jb,je)
     &     *za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju) + za(je,jd)*za(jn,je)*za(ju,jb)
     &     *zb(je,jb)*zb(je,ju) - za(jb,je)*za(je,jd)*za(ju,jn)*zb(je,jb)*zb(je,ju)
     &      - 2*za(jc,jd)*za(jd,jg)*za(jn,je)*zb(jd,je)*zb(jg,jc) + 2*za(jc,jd)*za(
     &     jd,jg)*za(ju,jn)*zb(jd,ju)*zb(jg,jc) - za(jd,jg)*za(jn,jc)*za(jn,je)*zb(
     &     je,jn)*zb(jg,jc) + za(jc,jd)*za(jn,je)*za(jn,jg)*zb(je,jn)*zb(jg,jc) - z
     &     a(jd,jg)*za(jn,je)*za(ju,jc)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jn,je)*z
     &     a(ju,jg)*zb(je,ju)*zb(jg,jc) - za(jd,jg)*za(je,jc)*za(ju,jn)*zb(je,ju)*z
     &     b(jg,jc) + za(jc,jd)*za(je,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc) + za(jn,jd)
     &     *za(jn,je)*za(ju,jb)*zb(je,ju)*zb(jn,jb) - za(jb,je)*za(jn,jd)*za(ju,jn)
     &     *zb(je,ju)*zb(jn,jb) - za(jd,jg)*za(jn,jc)*za(ju,jn)*zb(jg,jc)*zb(jn,ju)
     &      + za(jc,jd)*za(jn,jg)*za(ju,jn)*zb(jg,jc)*zb(jn,ju) + za(jb,jn)*(za(ju,
     &     jd)*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb))*(za(jn,je)*zb
     &     (je,jn) + za(ju,jn)*zb(jn,ju)) + za(jb,jd)*(za(jn,je)**2*zb(je,jb)*zb(je
     &     ,jn) + za(ju,jn)*(2*za(ju,jd)*zb(jb,ju)*zb(jd,ju) + za(ju,je)*zb(jb,ju)*
     &     zb(je,ju) + za(je,jd)*(-2*zb(jd,ju)*zb(je,jb) + zb(jd,jb)*zb(je,ju)) - 2
     &     *za(jn,jd)*zb(jd,ju)*zb(jn,jb) + za(ju,jn)*zb(jb,ju)*zb(jn,ju) + za(jn,j
     &     d)*zb(jd,jb)*zb(jn,ju)) + za(jn,je)*(2*za(je,jd)*zb(jd,je)*zb(je,jb) + z
     &     a(ju,jn)*zb(jb,ju)*zb(je,jn) + za(jn,jd)*zb(jd,jb)*zb(je,jn) + za(ju,je)
     &     *zb(je,jb)*zb(je,ju) + za(ju,jd)*(-2*zb(jb,ju)*zb(jd,je) + zb(jd,jb)*zb(
     &     je,ju)) + 2*za(jn,jd)*zb(jd,je)*zb(jn,jb) + za(ju,jn)*zb(je,jb)*zb(jn,ju
     &     ))))))/(9._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(je,jn) + s(je,
     &     ju) + s(jn,ju))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyGL_MPMM_M_L2

       function streal_heavyGL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_PMMM_P_L2 =
     &     (2*gb**2*propW34*zb(je,ju)*(2*im*imag(anomc7)*(za(jb,jd)*(za(jc,jg)*(-(z
     &     a(jn,je)*zb(jc,je)) + za(ju,jn)*zb(jc,ju)) + za(jb,jg)*(za(ju,jn)*zb(jb,
     &     ju) + za(jn,je)*zb(je,jb)))*zb(jg,jb) + za(jb,jg)*(za(ju,jn)*zb(jb,ju) +
     &      za(jn,je)*zb(je,jb))*(za(je,jd)*zb(jg,je) + za(jn,jd)*zb(jg,jn) + za(ju
     &     ,jd)*zb(jg,ju))) - 2*real(anomc7)*(za(jb,jd)*(za(jc,jg)*(-(za(jn,je)*zb(
     &     jc,je)) + za(ju,jn)*zb(jc,ju)) + za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,
     &     je)*zb(je,jb)))*zb(jg,jb) + za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*z
     &     b(je,jb))*(za(je,jd)*zb(jg,je) + za(jn,jd)*zb(jg,jn) + za(ju,jd)*zb(jg,j
     &     u))) + (3*im*imag(anomc4)*(za(jc,jg)*za(jd,jg)*(za(ju,jd)*(-2*za(ju,jn)*
     &     zb(jc,ju)*zb(jg,ju) + za(jn,je)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,j
     &     u))) + za(je,jd)*(2*za(jn,je)*zb(jc,je)*zb(jg,je) - za(ju,jn)*(zb(jc,ju)
     &     *zb(jg,je) + zb(jc,je)*zb(jg,ju))) + za(jn,jd)*(za(jn,je)*(zb(jc,jn)*zb(
     &     jg,je) + zb(jc,je)*zb(jg,jn)) - za(ju,jn)*(zb(jc,ju)*zb(jg,jn) + zb(jc,j
     &     n)*zb(jg,ju)))) + za(jb,jd)*(2*za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je
     &     )*zb(je,jb))*(za(ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(j
     &     n,jb)) + za(jc,jg)*(za(jn,je)*(2*za(je,jd)*zb(jc,je)*zb(je,jb) + za(ju,j
     &     d)*(-(zb(jb,ju)*zb(jc,je)) + zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jc,jn)
     &     *zb(je,jb) + zb(jc,je)*zb(jn,jb))) + za(ju,jn)*(2*za(ju,jd)*zb(jb,ju)*zb
     &     (jc,ju) + za(je,jd)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jn,
     &     jd)*(zb(jb,ju)*zb(jc,jn) - zb(jc,ju)*zb(jn,jb))))) + za(jb,jg)*(za(ju,jd
     &     )*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb))*((za(jn,je)*za(
     &     ju,jd) + za(je,jd)*za(ju,jn))*zb(je,ju) + za(jn,jd)*(za(jn,je)*zb(je,jn)
     &      + za(ju,jn)*zb(jn,ju)))))/za(jd,jg) + (3*real(anomc4)*(za(jc,jg)*za(jd,
     &     jg)*(-(za(ju,jd)*(-2*za(ju,jn)*zb(jc,ju)*zb(jg,ju) + za(jn,je)*(zb(jc,ju
     &     )*zb(jg,je) + zb(jc,je)*zb(jg,ju)))) + za(je,jd)*(-2*za(jn,je)*zb(jc,je)
     &     *zb(jg,je) + za(ju,jn)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju))) + za
     &     (jn,jd)*(-(za(jn,je)*(zb(jc,jn)*zb(jg,je) + zb(jc,je)*zb(jg,jn))) + za(j
     &     u,jn)*(zb(jc,ju)*zb(jg,jn) + zb(jc,jn)*zb(jg,ju)))) - za(jb,jd)*(2*za(jb
     &     ,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb))*(za(ju,jd)*zb(jb,ju) - 
     &     za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb)) + za(jc,jg)*(za(jn,je)*(2*za(
     &     je,jd)*zb(jc,je)*zb(je,jb) + za(ju,jd)*(-(zb(jb,ju)*zb(jc,je)) + zb(jc,j
     &     u)*zb(je,jb)) + za(jn,jd)*(zb(jc,jn)*zb(je,jb) + zb(jc,je)*zb(jn,jb))) +
     &      za(ju,jn)*(2*za(ju,jd)*zb(jb,ju)*zb(jc,ju) + za(je,jd)*(zb(jb,ju)*zb(jc
     &     ,je) - zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jb,ju)*zb(jc,jn) - zb(jc,ju)
     &     *zb(jn,jb))))) + za(jb,jg)*(-(za(ju,jd)*zb(jb,ju)) + za(je,jd)*zb(je,jb)
     &      + za(jn,jd)*zb(jn,jb))*((za(jn,je)*za(ju,jd) + za(je,jd)*za(ju,jn))*zb(
     &     je,ju) + za(jn,jd)*(za(jn,je)*zb(je,jn) + za(ju,jn)*zb(jn,ju)))))/za(jd,
     &     jg)))/(9._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc,jg))*(s(je,jn) + s(je,j
     &     u) + s(jn,ju))*za(jb,jg)*za(jc,jg))
       end function streal_heavyGL_PMMM_P_L2

       function streal_heavyGL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyGL_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)

           streal_heavyGL_PMMM_M_L2 =
     &     (2*gb**2*propW34*(im*imag(anomc4) - real(anomc4))*zb(jc,jb)*zb(je,ju)*(z
     &     a(jn,je)*(2*za(je,jd)*zb(jc,je)*zb(je,jb) + za(ju,jd)*(-(zb(jb,ju)*zb(jc
     &     ,je)) + zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jc,jn)*zb(je,jb) + zb(jc,je
     &     )*zb(jn,jb))) + za(ju,jn)*(2*za(ju,jd)*zb(jb,ju)*zb(jc,ju) + za(je,jd)*(
     &     zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jb,ju)*zb(jc,
     &     jn) - zb(jc,ju)*zb(jn,jb)))))/(3._dp*ecossin**2*(s(jb,jc) + s(jb,jg) + s(jc
     &     ,jg))*(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyGL_PMMM_M_L2

       function streal_lightZR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_MPMM_P_L2 =
     &     ((gb**2 - 3*gw**2)*propW34*propZ25*(-(im*imag(anomc4)) - real(anomc4))*z
     &     a(jb,jc)*za(jn,jd)*(za(ju,jd)*(-(za(jb,jd)*zb(jb,ju)*zb(jd,je)) + za(jc,
     &     jd)*zb(jc,ju)*zb(jd,je) + (za(jb,jn)*zb(jb,ju) + za(jn,jc)*zb(jc,ju))*zb
     &     (je,jn)) + za(jd,jg)*(-(za(jb,jd)*zb(jd,je)*zb(jg,jb)) + za(jb,jn)*zb(je
     &     ,jn)*zb(jg,jb) + (za(jc,jd)*zb(jd,je) + za(jn,jc)*zb(je,jn))*zb(jg,jc)) 
     &     + (za(jn,jd)*za(ju,jg)*zb(je,jn)*(za(jb,jc)*(zb(jc,ju)*zb(jg,jb) + zb(jb
     &     ,ju)*zb(jg,jc)) + (-(za(ju,jb)*zb(jb,ju)) + za(ju,jc)*zb(jc,ju))*zb(jg,j
     &     u)))/(s(jb,jc) + s(jb,ju) + s(jc,ju))))/(3._dp*ecossin**2*(s(jd,je) + s(jd,
     &     jn) + s(je,jn))*za(jd,jg)*za(ju,jg))
       end function streal_lightZR_MPMM_P_L2

       function streal_lightZR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_MPMM_M_L2 =
     &     ((gb**2 - 3*gw**2)*propW34*propZ25*(-(im*imag(anomc4)) - real(anomc4))*z
     &     a(jb,jc)*((za(ju,jb)*zb(jb,ju) - za(ju,jc)*zb(jc,ju))*zb(je,ju)*(za(jn,j
     &     d)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) + za(jb,jc)*(za(jn,jd)*(-(zb(jc,ju)*
     &     zb(jd,ju)*zb(je,jb)) + (zb(jb,ju)*((s(jd,je) + s(jd,jn) + s(je,jn))*zb(j
     &     c,je)*zb(jd,ju) + 2*za(jd,jg)*zb(jc,ju)*zb(jd,je)*zb(jg,jd) - 2*za(jn,jg
     &     )*zb(jc,ju)*zb(je,jn)*zb(jg,jd)))/(s(jd,je) + s(jd,jn) + s(je,jn))) + za
     &     (jn,jg)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb))*zb(jg,ju))))/(3._dp*eco
     &     ssin**2*(s(jb,jc) + s(jb,ju) + s(jc,ju))*zb(jg,jd)*zb(jg,ju))
       end function streal_lightZR_MPMM_M_L2

       function streal_lightZR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_PMMM_P_L2 =
     &     ((gb**2 - 3*gw**2)*propW34*propZ25*(-(im*imag(anomc4)) + real(anomc4))*z
     &     a(jn,jd)*zb(jc,jb)*(za(ju,jd)*(-(za(jb,jd)*zb(jb,ju)*zb(jd,je)) + za(jc,
     &     jd)*zb(jc,ju)*zb(jd,je) + (za(jb,jn)*zb(jb,ju) + za(jn,jc)*zb(jc,ju))*zb
     &     (je,jn)) + za(jd,jg)*(-(za(jb,jd)*zb(jd,je)*zb(jg,jb)) + za(jb,jn)*zb(je
     &     ,jn)*zb(jg,jb) + (za(jc,jd)*zb(jd,je) + za(jn,jc)*zb(je,jn))*zb(jg,jc)) 
     &     + (za(jn,jd)*za(ju,jg)*zb(je,jn)*(za(jb,jc)*(zb(jc,ju)*zb(jg,jb) + zb(jb
     &     ,ju)*zb(jg,jc)) + (-(za(ju,jb)*zb(jb,ju)) + za(ju,jc)*zb(jc,ju))*zb(jg,j
     &     u)))/(s(jb,jc) + s(jb,ju) + s(jc,ju))))/(3._dp*ecossin**2*(s(jd,je) + s(jd,
     &     jn) + s(je,jn))*za(jd,jg)*za(ju,jg))
       end function streal_lightZR_PMMM_P_L2

       function streal_lightZR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZR_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZR_PMMM_M_L2 =
     &     ((gb**2 - 3*gw**2)*propW34*propZ25*(-(im*imag(anomc4)) + real(anomc4))*z
     &     b(jc,jb)*((za(ju,jb)*zb(jb,ju) - za(ju,jc)*zb(jc,ju))*zb(je,ju)*(za(jn,j
     &     d)*zb(jd,ju) + za(jn,jg)*zb(jg,ju)) + za(jb,jc)*(za(jn,jd)*(-(zb(jc,ju)*
     &     zb(jd,ju)*zb(je,jb)) + (zb(jb,ju)*((s(jd,je) + s(jd,jn) + s(je,jn))*zb(j
     &     c,je)*zb(jd,ju) + 2*za(jd,jg)*zb(jc,ju)*zb(jd,je)*zb(jg,jd) - 2*za(jn,jg
     &     )*zb(jc,ju)*zb(je,jn)*zb(jg,jd)))/(s(jd,je) + s(jd,jn) + s(je,jn))) + za
     &     (jn,jg)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb))*zb(jg,ju))))/(3._dp*eco
     &     ssin**2*(s(jb,jc) + s(jb,ju) + s(jc,ju))*zb(jg,jd)*zb(jg,ju))
       end function streal_lightZR_PMMM_M_L2

       function streal_lightZL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_MPMM_P_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ25*(im*imag(anomc4) + real(anomc4))*za(j
     &     b,jc)*(za(ju,jd)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(jn,jd)*zb(jd,jb)) 
     &     - za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)))*zb(je,ju) + za(
     &     jd,jg)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(jn,jd)*zb(jd,jb)) - za(jc,jd
     &     )*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)))*zb(jg,je) + (2*za(jb,jd)*
     &     za(jc,jd)*za(ju,jg)*zb(jc,jb)*zb(je,ju)*(za(jn,je)*zb(jg,je) - za(ju,jn)
     &     *zb(jg,ju)))/(s(je,jn) + s(je,ju) + s(jn,ju))))/(3._dp*ecossin**2*(s(jb,jc)
     &      + s(jb,jd) + s(jc,jd))*za(jd,jg)*za(ju,jg))
       end function streal_lightZL_MPMM_P_L2

       function streal_lightZL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_MPMM_M_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ25*(-(im*imag(anomc4)) - real(anomc4))*z
     &     a(jb,jc)*zb(je,ju)*(za(ju,jn)*(za(jb,jd)*zb(jb,ju)*zb(jd,ju) - za(jc,jd)
     &     *zb(jc,ju)*zb(jd,ju) + (za(jb,jg)*zb(jb,ju) - za(jc,jg)*zb(jc,ju))*zb(jg
     &     ,ju)) + za(jn,je)*(za(jb,jd)*(zb(jd,ju)*zb(je,jb) + ((za(jc,jg)*zb(jc,jb
     &     ) + za(jd,jg)*zb(jd,jb))*zb(je,ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(j
     &     c,jd))) + za(jc,jd)*(zb(jc,je)*zb(jd,ju) + ((za(jb,jg)*zb(jc,jb) - za(jd
     &     ,jg)*zb(jd,jc))*zb(je,ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) +
     &      (za(jc,jg)*zb(jc,je) + za(jb,jg)*zb(je,jb))*zb(jg,ju))))/(3._dp*ecossin**2
     &     *(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jd)*zb(jg,ju))
       end function streal_lightZL_MPMM_M_L2

       function streal_lightZL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_PMMM_P_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ25*(im*imag(anomc4) - real(anomc4))*zb(j
     &     c,jb)*(za(ju,jd)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(jn,jd)*zb(jd,jb)) 
     &     - za(jc,jd)*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)))*zb(je,ju) + za(
     &     jd,jg)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb) + za(jn,jd)*zb(jd,jb)) - za(jc,jd
     &     )*(za(jb,jn)*zb(jc,jb) + za(jn,jd)*zb(jd,jc)))*zb(jg,je) + (2*za(jb,jd)*
     &     za(jc,jd)*za(ju,jg)*zb(jc,jb)*zb(je,ju)*(za(jn,je)*zb(jg,je) - za(ju,jn)
     &     *zb(jg,ju)))/(s(je,jn) + s(je,ju) + s(jn,ju))))/(3._dp*ecossin**2*(s(jb,jc)
     &      + s(jb,jd) + s(jc,jd))*za(jd,jg)*za(ju,jg))
       end function streal_lightZL_PMMM_P_L2

       function streal_lightZL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightZL_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ25

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ25  = 1._dp / (s(jb,jc) - zmass**2)

           streal_lightZL_PMMM_M_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ25*(-(im*imag(anomc4)) + real(anomc4))*z
     &     b(jc,jb)*zb(je,ju)*(za(ju,jn)*(za(jb,jd)*zb(jb,ju)*zb(jd,ju) - za(jc,jd)
     &     *zb(jc,ju)*zb(jd,ju) + (za(jb,jg)*zb(jb,ju) - za(jc,jg)*zb(jc,ju))*zb(jg
     &     ,ju)) + za(jn,je)*(za(jb,jd)*(zb(jd,ju)*zb(je,jb) + ((za(jc,jg)*zb(jc,jb
     &     ) + za(jd,jg)*zb(jd,jb))*zb(je,ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(j
     &     c,jd))) + za(jc,jd)*(zb(jc,je)*zb(jd,ju) + ((za(jb,jg)*zb(jc,jb) - za(jd
     &     ,jg)*zb(jd,jc))*zb(je,ju)*zb(jg,jd))/(s(jb,jc) + s(jb,jd) + s(jc,jd))) +
     &      (za(jc,jg)*zb(jc,je) + za(jb,jg)*zb(je,jb))*zb(jg,ju))))/(3._dp*ecossin**2
     &     *(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jd)*zb(jg,ju))
       end function streal_lightZL_PMMM_M_L2

       function streal_heavyZR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_MPMM_P_L2 =
     &     -((gb**2 - 3*gw**2)*propW34*propZ257*(im*imag(anomc4) + real(anomc4))*za
     &     (jb,jc)*za(jn,jd)*(za(jb,jc)**2*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,j
     &     b)) - (za(jb,jg)*za(ju,jc)*zb(je,ju) + za(jc,jg)*(za(ju,jb)*zb(je,ju) + 
     &     2*za(jb,jg)*zb(jg,je)))*zb(jg,ju) + za(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb(je
     &     ,ju) - za(ju,jc)*zb(jc,ju)*zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jg,je) - z
     &     a(jc,jg)*zb(jc,ju)*zb(jg,je) - za(jc,jg)*zb(jc,je)*zb(jg,ju) - za(jb,jg)
     &     *zb(je,jb)*zb(jg,ju))))/(3._dp*ecossin**2*(s(jd,je) + s(jd,jn) + s(je,jn))*
     &     za(jb,jg)*za(jc,jg))
       end function streal_heavyZR_MPMM_P_L2

       function streal_heavyZR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_MPMM_M_L2 =
     &     ((gb**2 - 3*gw**2)*propW34*propZ257*za(jn,jd)*((im*imag(anomc7) + real(a
     &     nomc7))*((gb**2 + 3*gw**2)*za(jb,jg)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(jg
     &     ,jc) + za(jb,jc)*(2*gb**2*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) + (gb*
     &     *2 + 3*gw**2)*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) + za(jc,jg)*(2*gb
     &     **2*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*(2*gb**2*zb(jc,j
     &     u)*zb(jg,jb) + (gb**2 + 3*gw**2)*zb(jb,ju)*zb(jg,jc))*zb(jg,je))) - 3*gw
     &     **2*im*imag(anomc4)*(za(jb,jc)**2*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc
     &     ,ju)*zb(je,jb)) + za(jb,jg)**2*zb(jg,jb)*(zb(jb,ju)*zb(jg,je) - zb(je,jb
     &     )*zb(jg,ju)) + za(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(j
     &     u,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(j
     &     g,jb) + 2*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*zb(jc,ju)*
     &     zb(je,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg
     &     )*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - 2*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,
     &     jc) - za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(jb,ju)*zb(j
     &     e,ju)*zb(jg,jc) + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*zb(jg,je) - za(jc,jg)*zb
     &     (jc,jb)*zb(jc,ju)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) - 
     &     za(jb,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju)) + za(jb,jg)*(za(ju,jb)*zb(jb,ju
     &     )*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) - za(jc,
     &     jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,
     &     ju) - za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) + za(jc,jg)*zb(je,jb)*zb(j
     &     g,jc)*zb(jg,ju) - 2*za(jc,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(ju,jc)*
     &     zb(je,ju)*(zb(jb,ju)*zb(jg,jc) + zb(jc,jb)*zb(jg,ju))) - za(jc,jg)*(za(j
     &     u,jb)*zb(je,ju)*(-(zb(jc,ju)*zb(jg,jb)) + zb(jc,jb)*zb(jg,ju)) + zb(jg,j
     &     c)*(za(ju,jc)*zb(jc,ju)*zb(je,ju) - za(ju,jg)*zb(je,ju)*zb(jg,ju) + za(j
     &     c,jg)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju))))) - 3*gw**2*real(anom
     &     c4)*(za(jb,jc)**2*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) 
     &     + za(jb,jg)**2*zb(jg,jb)*(zb(jb,ju)*zb(jg,je) - zb(je,jb)*zb(jg,ju)) + z
     &     a(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc,jb)*
     &     zb(jc,ju)*zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jb) + 2*za(jc,
     &     jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,
     &     jb) - za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jb,ju)*zb(j
     &     c,je)*zb(jg,jc) - 2*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc) - za(jc,jg)*
     &     zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) 
     &     + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju
     &     )*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(jc,
     &     jb)*zb(je,jb)*zb(jg,ju)) + za(jb,jg)*(za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(j
     &     g,jb) + za(jc,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) - za(jc,jg)*zb(jb,ju)*zb
     &     (jg,jc)*zb(jg,je) + za(jc,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) - za(ju,jg)*
     &     zb(je,ju)*zb(jg,jb)*zb(jg,ju) + za(jc,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) 
     &     - 2*za(jc,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) - za(ju,jc)*zb(je,ju)*(zb(jb
     &     ,ju)*zb(jg,jc) + zb(jc,jb)*zb(jg,ju))) - za(jc,jg)*(za(ju,jb)*zb(je,ju)*
     &     (-(zb(jc,ju)*zb(jg,jb)) + zb(jc,jb)*zb(jg,ju)) + zb(jg,jc)*(za(ju,jc)*zb
     &     (jc,ju)*zb(je,ju) - za(ju,jg)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*(zb(jc,ju)
     &     *zb(jg,je) + zb(jc,je)*zb(jg,ju)))))))/(9._dp*ecossin**2*gw**2*(s(jd,je) + 
     &     s(jd,jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZR_MPMM_M_L2

       function streal_heavyZR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_PMMM_P_L2 =
     &     -((gb**2 - 3*gw**2)*propW34*propZ257*za(jn,jd)*((im*imag(anomc7) - real(
     &     anomc7))*((gb**2 + 3*gw**2)*za(jb,jg)*za(ju,jc)*zb(jb,ju)*zb(je,ju)*zb(j
     &     g,jc) + za(jb,jc)*(2*gb**2*za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) + (gb
     &     **2 + 3*gw**2)*za(jb,jg)*zb(jb,ju)*zb(je,jb)*zb(jg,jc)) + za(jc,jg)*(2*g
     &     b**2*za(ju,jb)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jb,jg)*(2*gb**2*zb(jc,
     &     ju)*zb(jg,jb) + (gb**2 + 3*gw**2)*zb(jb,ju)*zb(jg,jc))*zb(jg,je))) + 3*g
     &     w**2*im*imag(anomc4)*(za(jb,jc)**2*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(j
     &     c,ju)*zb(je,jb)) + za(jb,jg)**2*zb(jg,jb)*(zb(jb,ju)*zb(jg,je) - zb(je,j
     &     b)*zb(jg,ju)) + za(jb,jc)*(za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(
     &     ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(
     &     jg,jb) - za(jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*z
     &     b(je,ju)*zb(jg,jb) + za(jc,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jc,jg)
     &     *zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc)
     &      + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,j
     &     u)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(jc
     &     ,jb)*zb(je,jb)*zb(jg,ju)) + za(jb,jg)*(za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(
     &     jg,jb) - za(jc,jg)*zb(jc,ju)*zb(jg,jb)*zb(jg,je) + za(jc,jg)*zb(jb,ju)*z
     &     b(jg,jc)*zb(jg,je) + za(jc,jg)*zb(jc,je)*zb(jg,jb)*zb(jg,ju) - za(ju,jg)
     &     *zb(je,ju)*zb(jg,jb)*zb(jg,ju) + za(jc,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju)
     &      - 2*za(jc,jg)*zb(jc,jb)*zb(jg,je)*zb(jg,ju) + za(ju,jc)*zb(je,ju)*(zb(j
     &     b,ju)*zb(jg,jc) - zb(jc,jb)*zb(jg,ju))) - za(jc,jg)*(za(ju,jb)*zb(je,ju)
     &     *(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*zb(jg,ju)) + zb(jg,jc)*(za(ju,jc)*zb(j
     &     c,ju)*zb(je,ju) - za(ju,jg)*zb(je,ju)*zb(jg,ju) + za(jc,jg)*(zb(jc,ju)*z
     &     b(jg,je) + zb(jc,je)*zb(jg,ju))))) - 3*gw**2*real(anomc4)*(za(jb,jc)**2*
     &     zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jb,jg)**2*zb(
     &     jg,jb)*(zb(jb,ju)*zb(jg,je) - zb(je,jb)*zb(jg,ju)) + za(jb,jc)*(za(ju,jb
     &     )*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc,jb)*zb(jc,ju)*zb(je,ju
     &     ) + za(jb,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jb) - za(jb,jg)*zb(jc,ju)*zb(je,
     &     jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(jg,jb) + za(jc,jg)*zb(j
     &     b,ju)*zb(jc,je)*zb(jg,jc) - za(jc,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jc) - za
     &     (ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)*zb(jb,ju)*zb(jc,jb)*zb
     &     (jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*
     &     zb(jc,je)*zb(jg,ju) - za(jb,jg)*zb(jc,jb)*zb(je,jb)*zb(jg,ju)) + za(jb,j
     &     g)*(za(ju,jb)*zb(jb,ju)*zb(je,ju)*zb(jg,jb) - za(jc,jg)*zb(jc,ju)*zb(jg,
     &     jb)*zb(jg,je) + za(jc,jg)*zb(jb,ju)*zb(jg,jc)*zb(jg,je) + za(jc,jg)*zb(j
     &     c,je)*zb(jg,jb)*zb(jg,ju) - za(ju,jg)*zb(je,ju)*zb(jg,jb)*zb(jg,ju) + za
     &     (jc,jg)*zb(je,jb)*zb(jg,jc)*zb(jg,ju) - 2*za(jc,jg)*zb(jc,jb)*zb(jg,je)*
     &     zb(jg,ju) + za(ju,jc)*zb(je,ju)*(zb(jb,ju)*zb(jg,jc) - zb(jc,jb)*zb(jg,j
     &     u))) - za(jc,jg)*(za(ju,jb)*zb(je,ju)*(zb(jc,ju)*zb(jg,jb) + zb(jc,jb)*z
     &     b(jg,ju)) + zb(jg,jc)*(za(ju,jc)*zb(jc,ju)*zb(je,ju) - za(ju,jg)*zb(je,j
     &     u)*zb(jg,ju) + za(jc,jg)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju))))))
     &     )/(9._dp*ecossin**2*gw**2*(s(jd,je) + s(jd,jn) + s(je,jn))*za(jb,jg)*za(jc,
     &     jg))
       end function streal_heavyZR_PMMM_P_L2

       function streal_heavyZR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZR_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZR_PMMM_M_L2 =
     &     ((gb**2 - 3*gw**2)*propW34*propZ257*(-(im*imag(anomc4)) + real(anomc4))*
     &     za(jn,jd)*zb(jc,jb)*(za(jb,jc)*zb(jc,jb)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju
     &     )*zb(je,jb)) + za(ju,jb)*zb(jb,ju)*zb(jc,jb)*zb(je,ju) - za(ju,jc)*zb(jc
     &     ,jb)*zb(jc,ju)*zb(je,ju) + za(jc,jg)*zb(jc,je)*zb(jc,ju)*zb(jg,jb) - za(
     &     jb,jg)*zb(jc,ju)*zb(je,jb)*zb(jg,jb) - za(ju,jg)*zb(jc,ju)*zb(je,ju)*zb(
     &     jg,jb) + za(jc,jg)*zb(jb,ju)*zb(jc,je)*zb(jg,jc) - za(jb,jg)*zb(jb,ju)*z
     &     b(je,jb)*zb(jg,jc) - za(ju,jg)*zb(jb,ju)*zb(je,ju)*zb(jg,jc) + za(jb,jg)
     &     *zb(jb,ju)*zb(jc,jb)*zb(jg,je) - za(jc,jg)*zb(jc,jb)*zb(jc,ju)*zb(jg,je)
     &     ))/(3._dp*ecossin**2*(s(jd,je) + s(jd,jn) + s(je,jn))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZR_PMMM_M_L2

       function streal_heavyZL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_MPMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_MPMM_P_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ257*(im*imag(anomc4) + real(anomc4))*za(
     &     jb,jc)*zb(je,ju)*(za(jn,je)*(za(jc,jd)*(-(za(jb,jn)*zb(je,jn)) + za(ju,j
     &     b)*zb(je,ju)) + za(jb,jd)*(2*za(jc,jd)*zb(jd,je) + za(jn,jc)*zb(je,jn) +
     &      za(ju,jc)*zb(je,ju))) + za(ju,jn)*(-(za(jc,jd)*(za(jb,je)*zb(je,ju) + z
     &     a(jb,jn)*zb(jn,ju))) + za(jb,jd)*(-2*za(jc,jd)*zb(jd,ju) + za(je,jc)*zb(
     &     je,ju) + za(jn,jc)*zb(jn,ju)))))/(3._dp*ecossin**2*(s(je,jn) + s(je,ju) + s
     &     (jn,ju))*za(jb,jg)*za(jc,jg))
       end function streal_heavyZL_MPMM_P_L2

       function streal_heavyZL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_MPMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_MPMM_M_L2 =
     &     -((gb**2 + 3*gw**2)*propW34*propZ257*zb(je,ju)*((im*imag(anomc7) + real(
     &     anomc7))*(-((gb**2 + 3*gw**2)*za(jb,jg)*za(jc,jd)*(za(ju,jn)*zb(jb,ju) +
     &      za(jn,je)*zb(je,jb))*zb(jg,jc)) + 2*gb**2*za(jb,jd)*zb(jg,jb)*(za(jd,jg
     &     )*(-(za(jn,je)*zb(jd,je)) + za(ju,jn)*zb(jd,ju)) + za(jb,jg)*(za(ju,jn)*
     &     zb(jb,ju) + za(jn,je)*zb(je,jb)) + za(jn,je)*za(jn,jg)*zb(je,jn) + za(jn
     &     ,je)*za(ju,jg)*zb(je,ju) + za(je,jg)*za(ju,jn)*zb(je,ju) + za(jn,jg)*za(
     &     ju,jn)*zb(jn,ju))) + 3*gw**2*im*imag(anomc4)*(-(za(jn,je)*za(ju,jb)*za(j
     &     u,jd)*zb(jb,ju)*zb(je,ju)) + za(jb,je)*za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(
     &     je,ju) + za(je,jd)*za(jn,je)*za(ju,jb)*zb(je,jb)*zb(je,ju) - za(jb,je)*z
     &     a(je,jd)*za(ju,jn)*zb(je,jb)*zb(je,ju) - 2*za(jc,jd)*za(jd,jg)*za(jn,je)
     &     *zb(jd,je)*zb(jg,jc) + 2*za(jc,jd)*za(jd,jg)*za(ju,jn)*zb(jd,ju)*zb(jg,j
     &     c) - za(jd,jg)*za(jn,jc)*za(jn,je)*zb(je,jn)*zb(jg,jc) + za(jc,jd)*za(jn
     &     ,je)*za(jn,jg)*zb(je,jn)*zb(jg,jc) - za(jd,jg)*za(jn,je)*za(ju,jc)*zb(je
     &     ,ju)*zb(jg,jc) + za(jc,jd)*za(jn,je)*za(ju,jg)*zb(je,ju)*zb(jg,jc) - za(
     &     jd,jg)*za(je,jc)*za(ju,jn)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(je,jg)*za(
     &     ju,jn)*zb(je,ju)*zb(jg,jc) + za(jn,jd)*za(jn,je)*za(ju,jb)*zb(je,ju)*zb(
     &     jn,jb) - za(jb,je)*za(jn,jd)*za(ju,jn)*zb(je,ju)*zb(jn,jb) - za(jd,jg)*z
     &     a(jn,jc)*za(ju,jn)*zb(jg,jc)*zb(jn,ju) + za(jc,jd)*za(jn,jg)*za(ju,jn)*z
     &     b(jg,jc)*zb(jn,ju) + za(jb,jn)*(za(ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,jb
     &     ) - za(jn,jd)*zb(jn,jb))*(za(jn,je)*zb(je,jn) + za(ju,jn)*zb(jn,ju)) + z
     &     a(jb,jd)*(za(jn,je)**2*zb(je,jb)*zb(je,jn) + za(ju,jn)*(2*za(ju,jd)*zb(j
     &     b,ju)*zb(jd,ju) + za(ju,je)*zb(jb,ju)*zb(je,ju) + za(je,jd)*(-2*zb(jd,ju
     &     )*zb(je,jb) + zb(jd,jb)*zb(je,ju)) - 2*za(jn,jd)*zb(jd,ju)*zb(jn,jb) + z
     &     a(ju,jn)*zb(jb,ju)*zb(jn,ju) + za(jn,jd)*zb(jd,jb)*zb(jn,ju)) + za(jn,je
     &     )*(2*za(je,jd)*zb(jd,je)*zb(je,jb) + za(ju,jn)*zb(jb,ju)*zb(je,jn) + za(
     &     jn,jd)*zb(jd,jb)*zb(je,jn) + za(ju,je)*zb(je,jb)*zb(je,ju) + za(ju,jd)*(
     &     -2*zb(jb,ju)*zb(jd,je) + zb(jd,jb)*zb(je,ju)) + 2*za(jn,jd)*zb(jd,je)*zb
     &     (jn,jb) + za(ju,jn)*zb(je,jb)*zb(jn,ju)))) + 3*gw**2*real(anomc4)*(-(za(
     &     jn,je)*za(ju,jb)*za(ju,jd)*zb(jb,ju)*zb(je,ju)) + za(jb,je)*za(ju,jd)*za
     &     (ju,jn)*zb(jb,ju)*zb(je,ju) + za(je,jd)*za(jn,je)*za(ju,jb)*zb(je,jb)*zb
     &     (je,ju) - za(jb,je)*za(je,jd)*za(ju,jn)*zb(je,jb)*zb(je,ju) - 2*za(jc,jd
     &     )*za(jd,jg)*za(jn,je)*zb(jd,je)*zb(jg,jc) + 2*za(jc,jd)*za(jd,jg)*za(ju,
     &     jn)*zb(jd,ju)*zb(jg,jc) - za(jd,jg)*za(jn,jc)*za(jn,je)*zb(je,jn)*zb(jg,
     &     jc) + za(jc,jd)*za(jn,je)*za(jn,jg)*zb(je,jn)*zb(jg,jc) - za(jd,jg)*za(j
     &     n,je)*za(ju,jc)*zb(je,ju)*zb(jg,jc) + za(jc,jd)*za(jn,je)*za(ju,jg)*zb(j
     &     e,ju)*zb(jg,jc) - za(jd,jg)*za(je,jc)*za(ju,jn)*zb(je,ju)*zb(jg,jc) + za
     &     (jc,jd)*za(je,jg)*za(ju,jn)*zb(je,ju)*zb(jg,jc) + za(jn,jd)*za(jn,je)*za
     &     (ju,jb)*zb(je,ju)*zb(jn,jb) - za(jb,je)*za(jn,jd)*za(ju,jn)*zb(je,ju)*zb
     &     (jn,jb) - za(jd,jg)*za(jn,jc)*za(ju,jn)*zb(jg,jc)*zb(jn,ju) + za(jc,jd)*
     &     za(jn,jg)*za(ju,jn)*zb(jg,jc)*zb(jn,ju) + za(jb,jn)*(za(ju,jd)*zb(jb,ju)
     &      - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb))*(za(jn,je)*zb(je,jn) + za(
     &     ju,jn)*zb(jn,ju)) + za(jb,jd)*(za(jn,je)**2*zb(je,jb)*zb(je,jn) + za(ju,
     &     jn)*(2*za(ju,jd)*zb(jb,ju)*zb(jd,ju) + za(ju,je)*zb(jb,ju)*zb(je,ju) + z
     &     a(je,jd)*(-2*zb(jd,ju)*zb(je,jb) + zb(jd,jb)*zb(je,ju)) - 2*za(jn,jd)*zb
     &     (jd,ju)*zb(jn,jb) + za(ju,jn)*zb(jb,ju)*zb(jn,ju) + za(jn,jd)*zb(jd,jb)*
     &     zb(jn,ju)) + za(jn,je)*(2*za(je,jd)*zb(jd,je)*zb(je,jb) + za(ju,jn)*zb(j
     &     b,ju)*zb(je,jn) + za(jn,jd)*zb(jd,jb)*zb(je,jn) + za(ju,je)*zb(je,jb)*zb
     &     (je,ju) + za(ju,jd)*(-2*zb(jb,ju)*zb(jd,je) + zb(jd,jb)*zb(je,ju)) + 2*z
     &     a(jn,jd)*zb(jd,je)*zb(jn,jb) + za(ju,jn)*zb(je,jb)*zb(jn,ju))))))/(9._dp*ec
     &     ossin**2*gw**2*(s(je,jn) + s(je,ju) + s(jn,ju))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZL_MPMM_M_L2

       function streal_heavyZL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_PMMM_P_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_PMMM_P_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ257*zb(je,ju)*((im*imag(anomc7)*(za(jb,j
     &     d)*(2*gb**2*za(jc,jg)*(za(jn,je)*zb(jc,je) - za(ju,jn)*zb(jc,ju)) + (gb*
     &     *2 + 3*gw**2)*za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb)))*zb(
     &     jg,jb) + (gb**2 + 3*gw**2)*za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb
     &     (je,jb))*(za(je,jd)*zb(jg,je) + za(jn,jd)*zb(jg,jn) + za(ju,jd)*zb(jg,ju
     &     ))))/gw**2 - (real(anomc7)*(za(jb,jd)*(2*gb**2*za(jc,jg)*(za(jn,je)*zb(j
     &     c,je) - za(ju,jn)*zb(jc,ju)) + (gb**2 + 3*gw**2)*za(jb,jg)*(za(ju,jn)*zb
     &     (jb,ju) + za(jn,je)*zb(je,jb)))*zb(jg,jb) + (gb**2 + 3*gw**2)*za(jb,jg)*
     &     (za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb))*(za(je,jd)*zb(jg,je) + za(jn
     &     ,jd)*zb(jg,jn) + za(ju,jd)*zb(jg,ju))))/gw**2 + (3*im*imag(anomc4)*(za(j
     &     c,jg)*za(jd,jg)*(za(ju,jd)*(-2*za(ju,jn)*zb(jc,ju)*zb(jg,ju) + za(jn,je)
     &     *(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju))) + za(je,jd)*(2*za(jn,je)*z
     &     b(jc,je)*zb(jg,je) - za(ju,jn)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,ju
     &     ))) + za(jn,jd)*(za(jn,je)*(zb(jc,jn)*zb(jg,je) + zb(jc,je)*zb(jg,jn)) -
     &      za(ju,jn)*(zb(jc,ju)*zb(jg,jn) + zb(jc,jn)*zb(jg,ju)))) + za(jb,jd)*(2*
     &     za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(jn,je)*zb(je,jb))*(za(ju,jd)*zb(jb,j
     &     u) - za(je,jd)*zb(je,jb) - za(jn,jd)*zb(jn,jb)) + za(jc,jg)*(za(jn,je)*(
     &     2*za(je,jd)*zb(jc,je)*zb(je,jb) + za(ju,jd)*(-(zb(jb,ju)*zb(jc,je)) + zb
     &     (jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jc,jn)*zb(je,jb) + zb(jc,je)*zb(jn,jb
     &     ))) + za(ju,jn)*(2*za(ju,jd)*zb(jb,ju)*zb(jc,ju) + za(je,jd)*(zb(jb,ju)*
     &     zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jb,ju)*zb(jc,jn) - zb(j
     &     c,ju)*zb(jn,jb))))) + za(jb,jg)*(za(ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,j
     &     b) - za(jn,jd)*zb(jn,jb))*((za(jn,je)*za(ju,jd) + za(je,jd)*za(ju,jn))*z
     &     b(je,ju) + za(jn,jd)*(za(jn,je)*zb(je,jn) + za(ju,jn)*zb(jn,ju)))))/za(j
     &     d,jg) + (3*real(anomc4)*(za(jc,jg)*za(jd,jg)*(-(za(ju,jd)*(-2*za(ju,jn)*
     &     zb(jc,ju)*zb(jg,ju) + za(jn,je)*(zb(jc,ju)*zb(jg,je) + zb(jc,je)*zb(jg,j
     &     u)))) + za(je,jd)*(-2*za(jn,je)*zb(jc,je)*zb(jg,je) + za(ju,jn)*(zb(jc,j
     &     u)*zb(jg,je) + zb(jc,je)*zb(jg,ju))) + za(jn,jd)*(-(za(jn,je)*(zb(jc,jn)
     &     *zb(jg,je) + zb(jc,je)*zb(jg,jn))) + za(ju,jn)*(zb(jc,ju)*zb(jg,jn) + zb
     &     (jc,jn)*zb(jg,ju)))) - za(jb,jd)*(2*za(jb,jg)*(za(ju,jn)*zb(jb,ju) + za(
     &     jn,je)*zb(je,jb))*(za(ju,jd)*zb(jb,ju) - za(je,jd)*zb(je,jb) - za(jn,jd)
     &     *zb(jn,jb)) + za(jc,jg)*(za(jn,je)*(2*za(je,jd)*zb(jc,je)*zb(je,jb) + za
     &     (ju,jd)*(-(zb(jb,ju)*zb(jc,je)) + zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(j
     &     c,jn)*zb(je,jb) + zb(jc,je)*zb(jn,jb))) + za(ju,jn)*(2*za(ju,jd)*zb(jb,j
     &     u)*zb(jc,ju) + za(je,jd)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + z
     &     a(jn,jd)*(zb(jb,ju)*zb(jc,jn) - zb(jc,ju)*zb(jn,jb))))) + za(jb,jg)*(-(z
     &     a(ju,jd)*zb(jb,ju)) + za(je,jd)*zb(je,jb) + za(jn,jd)*zb(jn,jb))*((za(jn
     &     ,je)*za(ju,jd) + za(je,jd)*za(ju,jn))*zb(je,ju) + za(jn,jd)*(za(jn,je)*z
     &     b(je,jn) + za(ju,jn)*zb(jn,ju)))))/za(jd,jg)))/(9._dp*ecossin**2*(s(je,jn) 
     &     + s(je,ju) + s(jn,ju))*za(jb,jg)*za(jc,jg))
       end function streal_heavyZL_PMMM_P_L2

       function streal_heavyZL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg, za,zb)
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_heavyZL_PMMM_M_L2
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propW34
           real(dp) :: propZ257

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propZ257 = 1._dp / (s(jb,jc)+s(jb,jg)+s(jc,jg) - zmass**2)

           streal_heavyZL_PMMM_M_L2 =
     &     ((gb**2 + 3*gw**2)*propW34*propZ257*(im*imag(anomc4) - real(anomc4))*zb(
     &     jc,jb)*zb(je,ju)*(za(jn,je)*(2*za(je,jd)*zb(jc,je)*zb(je,jb) + za(ju,jd)
     &     *(-(zb(jb,ju)*zb(jc,je)) + zb(jc,ju)*zb(je,jb)) + za(jn,jd)*(zb(jc,jn)*z
     &     b(je,jb) + zb(jc,je)*zb(jn,jb))) + za(ju,jn)*(2*za(ju,jd)*zb(jb,ju)*zb(j
     &     c,ju) + za(je,jd)*(zb(jb,ju)*zb(jc,je) - zb(jc,ju)*zb(je,jb)) + za(jn,jd
     &     )*(zb(jb,ju)*zb(jc,jn) - zb(jc,ju)*zb(jn,jb)))))/(3._dp*ecossin**2*(s(je,jn
     &     ) + s(je,ju) + s(jn,ju))*zb(jg,jb)*zb(jg,jc))
       end function streal_heavyZL_PMMM_M_L2

