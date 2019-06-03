      module zaj_virtamps_m
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'
        include 'scale.f' !impure
        include 'epinv.f' !impure
        include 'epinv2.f' !impure

      contains

c     virt helicity amplitude (q;1), couplings stripped,
c       receiving contributions from different leading color
c       and subleading color primitive amps, each of it
c       having a pole, and cc, sc part (see BDK paper)
      subroutine zaj_virtamp_q1(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)
        complex(dp) :: amps_lc(2,2,2,2), amps_slc(2,2,2,2)

        call zaj_virt_a6_lc(i1,i2,i3,i4,i5,i6,za,zb,amps_lc)

        call zaj_virt_a6_slc(i1,i2,i3,i4,i5,i6,za,zb,amps_slc)

        amps = -Nc*amps_lc - 1._dp/Nc*amps_slc

      end subroutine

c     virt helicity amplitude (q;2), couplings stripped
      subroutine zaj_virtamp_q2(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)
        complex(dp) :: amps_fvs(2,2,2,2), amps_fvf(2,2,2,2)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a64v_fvs_pp, zaj_virt_a64v_fvs_pm,
     &     zaj_virt_a64v_fvs_mp, amps_fvs)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a64v_fvf_pp, zaj_virt_a64v_fvf_pm,
     &     zaj_virt_a64v_fvf_mp, amps_fvf)

        ! apparently there is a 3<->4 symmetry,
        ! such that a factor 2 is sufficient
        amps = -2*(amps_fvs + amps_fvf)

      end subroutine 

      subroutine zaj_virtamp_q3(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)
        complex(dp) :: amps_s1(2,2,2,2), amps_s2(2,2,2,2)

       call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &    zaj_virt_a64ax_pp, zaj_virt_a64ax_pm,
     &    zaj_virt_a64ax_mp, amps_s1)

       call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &    zaj_virt_a64ax_pp, zaj_virt_a64ax_mp,
     &    zaj_virt_a64ax_pm, amps_s2)

        amps = amps_s1 + amps_s2

      end subroutine

      ! assembly of q;1 leading color (lc) piece
      subroutine zaj_virt_a6_lc(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)
        complex(dp) :: amps_pole(2,2,2,2), amps_fcc(2,2,2,2), amps_fsc(2,2,2,2)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_lc_pole_pp, zaj_virt_a6_lc_pole_pm,
     &     zaj_virt_a6_lc_pole_mp, amps_pole)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_lc_fcc_pp, zaj_virt_a6_lc_fcc_pm,
     &     zaj_virt_a6_lc_fcc_mp, amps_fcc)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_lc_fsc_pp, zaj_virt_a6_lc_fsc_pm,
     &     zaj_virt_a6_lc_fsc_mp, amps_fsc)

        amps = amps_pole + amps_fcc + amps_fsc

      end subroutine

      subroutine zaj_virt_a6_slc(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)
        complex(dp) :: amps_pole(2,2,2,2), amps_fcc(2,2,2,2), amps_fsc(2,2,2,2)
        complex(dp) :: amps_pole_fl(2,2,2,2), amps_fcc_fl(2,2,2,2), amps_fsc_fl(2,2,2,2)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_slc_pole_pp, zaj_virt_a6_slc_pole_pm,
     &     zaj_virt_a6_slc_pole_mp, amps_pole)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fcc_pp, zaj_virt_a6_slc_fcc_pm,
     &     zaj_virt_a6_slc_fcc_mp, amps_fcc)

        call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fsc_pp, zaj_virt_a6_slc_fsc_pm,
     &     zaj_virt_a6_slc_fsc_mp, amps_fsc)

c       symmetrization configuration
        call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &     zaj_virt_a6_slc_pole_pp, zaj_virt_a6_slc_pole_mp,
     &     zaj_virt_a6_slc_pole_pm, amps_pole_fl)

        call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fcc_pp, zaj_virt_a6_slc_fcc_mp,
     &     zaj_virt_a6_slc_fcc_pm, amps_fcc_fl)

        call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fsc_pp, zaj_virt_a6_slc_fsc_mp,
     &     zaj_virt_a6_slc_fsc_pm, amps_fsc_fl)

        amps = amps_pole + amps_fcc + amps_fsc
     &        + amps_pole_fl + amps_fcc_fl + amps_fsc_fl

      end subroutine
        
 


c     primitive amplitudes for q;1 subleading color piece
      function zaj_virt_a6_slc_pole_pp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_pole_pp
        complex(dp) :: zaj_virt_a6_slc_tree_pp
        complex(dp) :: Vcc, Vsc

        complex(dp) :: xl12, xl56
        complex(dp) :: lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        xl12=Lnrat(musq,real(-s(i1,i2),dp))
        xl56=Lnrat(musq,real(-s(i5,i6),dp))

        Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
        Vsc=half*(one+epinv+xl56)

        zaj_virt_a6_slc_tree_pp =
     &    (-za(i2,i5)**2)/(za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))

        zaj_virt_a6_slc_pole_pp = zaj_virt_a6_slc_tree_pp*(Vcc+Vsc)

      end function

      function zaj_virt_a6_slc_pole_pm(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_pole_pm
        complex(dp) :: zaj_virt_a6_slc_tree_pm
        complex(dp) :: Vcc, Vsc

        complex(dp) :: xl12, xl56
        complex(dp) :: lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        xl12=Lnrat(musq,real(-s(i1,i2),dp))
        xl56=Lnrat(musq,real(-s(i5,i6),dp))

        zaj_virt_a6_slc_tree_pm =
     .    -(-((za(i2,i4)*za(i2,i5)*zb(i1,i3)*zb(i1,i6))/
     .    (s(i3,i4)*s(i5,i6)*za(i2,i3)*zb(i1,i4)))+
     .    (za(i2,i5)*zb(i1,i3)**2*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6
     .    )))/
     .    (s(i3,i4)*s(i5,i6)*zb(i1,i4)*t(i1,i3,i4))+
     .    (za(i2,i4)**2*zb(i1,i6)*(-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4
     .    )))/
     .    (s(i3,i4)*s(i5,i6)*za(i2,i3)*t(i2,i3,i4)))

        Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
        Vsc=half*(one+epinv+xl56)

        zaj_virt_a6_slc_pole_pm = zaj_virt_a6_slc_tree_pm*(Vcc+Vsc)

      end function

      function zaj_virt_a6_slc_pole_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_pole_mp
        complex(dp) :: zaj_virt_a6_slc_tree_mp
        complex(dp) :: Vcc, Vsc

        complex(dp) :: xl12, xl56
        complex(dp) :: lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        xl12=Lnrat(musq,real(-s(i1,i2),dp))
        xl56=Lnrat(musq,real(-s(i5,i6),dp))

        zaj_virt_a6_slc_tree_mp =
     &    (((-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i3,i4)*s(i5,i6)*za(i1,i4)*zb(i2,i3))-
     &    (za(i1,i3)*za(i2,i5)*zb(i1,i4)*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i3,i4)*s(i5,i6)*za(i1,i4)*t(i1,i3,i4))-
     &    (za(i2,i3)*zb(i1,i6)*zb(i2,i4)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    (s(i3,i4)*s(i5,i6)*zb(i2,i3)*t(i2,i3,i4)))

        Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
        Vsc=half*(one+epinv+xl56)

        zaj_virt_a6_slc_pole_mp = zaj_virt_a6_slc_tree_mp*(Vcc+Vsc)

      end function

      function zaj_virt_a6_slc_fcc_pp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_fcc_pp

        complex(dp) :: L0,Lsm1,Lsm1_2me

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_virt_a6_slc_fcc_pp =
     &    ((Lsm1(-s(i1,i4),-t(i1,i2,i4),-s(i1,i2),-t(i1,i2,i4))+
     &    Lsm1_2me(t(i1,i3,i4),t(i1,i2,i4),s(i1,i4),s(i5,i6)))*za(i2,i5)**2)
     &    /
     &    (za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))-
     &    (Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))*za(i2,i5)*
     &    (za(i1,i5)*za(i2,i3)-za(i1,i2)*za(i3,i5)))/
     &    (za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1_2me(t(i1,i2,i3),t(i2,i3,i4),s(i2,i3),s(i5,i6))*za(i2,i5)*
     &    (-(za(i1,i5)*za(i2,i4))+za(i1,i2)*za(i4,i5)))/
     &    (za(i1,i4)**2*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1_2me(t(i1,i2,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*za(i2,i5)*
     &    (za(i2,i4)*za(i3,i5)+za(i2,i3)*za(i4,i5)))/
     &    (za(i1,i4)*za(i2,i3)*za(i3,i4)**2*za(i5,i6))-
     &    (2._dp*L0(-t(i2,i3,i4),-s(i5,i6))*za(i1,i5)*za(i2,i5)*
     &    (-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4)))/
     &    (s(i5,i6)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))-
     &    (2._dp*L0(-t(i2,i3,i4),-s(i2,i3))*za(i2,i5)*za(i4,i5)*zb(i3,i4))/
     &    (s(i2,i3)*za(i1,i4)*za(i3,i4)*za(i5,i6))

      end function

      function zaj_virt_a6_slc_fcc_pm(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_fcc_pm

        complex(dp)  :: L0,Lsm1,Lsm1_2mh,Lnrat,I3m

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_virt_a6_slc_fcc_pm=
     &    Lsm1(-s(i1,i2),-t(i1,i2,i4),-s(i1,i4),-t(i1,i2,i4))*
     &    (-((zb(i1,i2)**2*(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))**2)/
     &    (za(i5,i6)*zb(i1,i4)*zb(i2,i4)**2*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))+
     &    (zb(i1,i4)*(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))**2)/
     &    (za(i5,i6)*zb(i2,i4)**2*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4
     &    ))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))-
     &    (2._dp*L0(-t(i2,i3,i4),-s(i5,i6))*za(i1,i5)*zb(i1,i3)*
     &    (za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4))*
     &    (-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4)))/
     &    (s(i5,i6)*za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     &    Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))*
     &    ((za(i1,i2)**2*(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))**2)/
     &    (za(i1,i3)**2*za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-
     &    za(i2,i3)*zb(i2,i4))*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    zb(i5,i6))-(za(i2,i3)*
     &    (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     &    (za(i1,i3)**2*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6)))-
     &    (2._dp*L0(-t(i1,i3,i4),-s(i5,i6))*za(i2,i4)*
     &    (-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*zb(i2,i6)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6)))/
     &    (s(i5,i6)*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6))
     
          zaj_virt_a6_slc_fcc_pm=zaj_virt_a6_slc_fcc_pm+
     &    I3m(s(i1,i4),s(i2,i3),s(i5,i6))*
     &    ((za(i1,i5)*za(i2,i4)*zb(i1,i3)*zb(i2,i6)+
     &    ((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    (za(i1,i3)*za(i2,i5)*zb(i1,i6)*zb(i2,i4)+
     &    za(i1,i3)*zb(i1,i4)*
     &    (-(za(i2,i5)*zb(i2,i6))+za(i3,i5)*zb(i3,i6))))/
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))+
     &    ((-(za(i1,i4)*zb(i1,i3))+za(i2,i4)*zb(i2,i3))*
     &    (-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6)))/2)/
     &    ((-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))-
     &    (za(i4,i5)*zb(i1,i3)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))
     &    /
     &    ((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*t(i1,i2,i3)))+
     &    Lsm1_2mh(s(i1,i4),t(i1,i2,i3),s(i2,i3),s(i5,i6))*
     &    (-(((-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))**2*
     &    (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     &    (za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**3*zb(i5,i6)))+
     &    (za(i1,i2)**2*zb(i4,i6)**2*t(i1,i2,i3)**2)/
     &    (za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**3*zb(i5,i6)))+
     &    Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*
     &    (((-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))**2*
     &    (-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))**2)/
     &    (za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**3*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))-
     &    (za(i2,i3)*zb(i4,i6)**2*t(i1,i2,i3)**2)/
     &    ((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**3*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6)))+
     &    (2._dp*Lnrat(-s(i2,i3),-s(i5,i6))*(-2._dp*(s(i1,i2)-s(i3,i4))*za(i2,i4
     &    )*zb(i2,i3)*
     &    (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6))+
     &    (za(i2,i5)*zb(i1,i3)*((s(i1,i2)-s(i3,i4))*
     &    ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i3,i5)*zb(i2,i3)+
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*zb(i2,i6))+
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i5)*zb(i2,i3))+
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*zb(i3,i6))))/
     &    (za(i5,i6)*zb(i1,i4))-
     &    (za(i2,i4)*(zb(i1,i3)*
     &    (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))-
     &    zb(i3,i4)*(-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6))-
     &    zb(i2,i3)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))*
     &    (za(i3,i5)*zb(i2,i3)*zb(i5,i6)+zb(i2,i6)*t(i1,i2,i3)))/zb(i5,i6)
     &    ))/((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(
     &    i5,i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
          
          zaj_virt_a6_slc_fcc_pm=zaj_virt_a6_slc_fcc_pm+
     &    I3m(s(i2,i3),s(i1,i4),s(i5,i6))*
     &    ((za(i1,i5)*za(i2,i4)*zb(i1,i3)*zb(i2,i6)+
     &    ((za(i1,i4)*zb(i1,i3)-za(i2,i4)*zb(i2,i3))*
     &    (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/2+
     &    ((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    (za(i1,i3)*za(i2,i5)*zb(i1,i6)*zb(i2,i4)+
     &    za(i2,i3)*zb(i2,i4)*
     &    (-(za(i1,i5)*zb(i1,i6))+za(i4,i5)*zb(i4,i6))))/
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))/
     &    ((-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))-
     &    (za(i2,i4)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*zb(i3,i6))/
     &    ((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*t(i1,i2,i4)))+
     &    Lsm1_2mh(s(i2,i3),t(i1,i2,i4),s(i1,i4),s(i5,i6))*
     &    (((za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))**2*
     &    (-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))**2)/
     &    (za(i5,i6)*zb(i1,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))**3)-
     &    (za(i3,i5)**2*zb(i1,i2)**2*t(i1,i2,i4)**2)/
     &    (za(i5,i6)*zb(i1,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))**3))+
     &    Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),s(i5,i6))*
     &    (-(((za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))**2*
     &    (-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))**2)/
     &    (za(i5,i6)*zb(i1,i4)*(-(za(i1,i3)*zb(i1,i4))-
     &    za(i2,i3)*zb(i2,i4))**3*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))+
     &    (za(i3,i5)**2*zb(i1,i4)*t(i1,i2,i4)**2)/
     &    (za(i5,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**3*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))+
     &    (2._dp*Lnrat(-s(i1,i4),-s(i5,i6))*(-2._dp*(s(i1,i2)-s(i3,i4))*za(i1,i4
     &    )*zb(i1,i3)*
     &    (-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6))-
     &    (za(i2,i4)*zb(i1,i6)*((s(i1,i2)-s(i3,i4))*
     &    ((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i4,i6)-
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*zb(i5,i6))+
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    (-((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i1,i6))-
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i4,i5)*zb(i5,i6))))/
     &    (za(i2,i3)*zb(i5,i6))+
     &    (zb(i1,i3)*(-(za(i1,i4)*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3)))+
     &    za(i2,i4)*(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))+
     &    za(i3,i4)*(-(za(i1,i5)*zb(i1,i3))+za(i4,i5)*zb(i3,i4)))*
     &    (-(za(i1,i4)*za(i5,i6)*zb(i4,i6))+za(i1,i5)*t(i1,i2,i4)))/
     &    za(i5,i6)))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
     
          zaj_virt_a6_slc_fcc_pm=zaj_virt_a6_slc_fcc_pm+
     &    Lsm1(-s(i1,i4),-t(i1,i3,i4),-s(i3,i4),-t(i1,i3,i4))*
     &    ((za(i3,i4)**2*(za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))**2)/
     &    (za(i1,i3)**3*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    zb(i5,i6)*t(i1,i3,i4))-
     &    (za(i1,i4)**2*(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6))**2)/
     &    (za(i1,i3)**3*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i
     &    6)*
     &    t(i1,i3,i4)))+(2._dp*L0(-t(i1,i3,i4),-s(i3,i4))*za(i1,i4)*zb(i1,i3)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6))*
     &    (za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6)))/
     &    (s(i3,i4)*za(i1,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    zb(i5,i6)*t(i1,i3,i4))+(2._dp*L0(-t(i1,i3,i4),-s(i1,i4))*za(i1,i4)*z
     &    b(i1,i3)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6))*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i1,i4)*za(i1,i3)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    zb(i5,i6)*t(i1,i3,i4))+Lsm1_2mh(s(i1,i2),t(i1,i3,i4),s(i3,i4),
     &    s(i5,i6))*(-((za(i2,i5)**2*zb(i1,i3)**3)/
     &    (za(i5,i6)*zb(i1,i4)*(-(za(i2,i3)*zb(i1,i3))-
     &    za(i2,i4)*zb(i1,i4))*zb(i3,i4)*t(i1,i3,i4)))-
     &    ((-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))**3*
     &    (za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))**2)/
     &    (za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**3*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6)*t(i1,i3,i4)
     &    )+(za(i1,i4)**2*(-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*zb(i2
     &    ,i6)**2*
     &    t(i1,i3,i4))/
     &    (za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**3*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6)))+
     &    (I3m(s(i1,i4),s(i2,i3),s(i5,i6))*zb(i1,i3)*
     &    (-(s(i5,i6)*za(i2,i4)*((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*
     &    zb(i2,i6)-2._dp*za(i1,i4)*za(i3,i5)*zb(i2,i3)*zb(i4,i6)))+
     &    2._dp*za(i1,i4)*za(i1,i5)*za(i2,i5)*
     &    ((-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i1,i2)-
     &    2._dp*za(i3,i4)*zb(i1,i4)*zb(i2,i3))*zb(i5,i6)-
     &    2._dp*za(i1,i4)*za(i4,i5)*zb(i4,i6)*
     &    (s(i5,i6)*(-s(i1,i4)-s(i2,i3)+s(i5,i6))+
     &    (s(i1,i4)-s(i2,i3)-s(i5,i6))*t(i1,i3,i4))))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
     
          zaj_virt_a6_slc_fcc_pm=zaj_virt_a6_slc_fcc_pm+
     &    Lsm1(-s(i2,i3),-t(i2,i3,i4),-s(i3,i4),-t(i2,i3,i4))*
     &    (-(((za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2*zb(i3,i4)**2)/
     &    (za(i5,i6)*zb(i2,i4)**3*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    t(i2,i3,i4)))+(zb(i2,i3)**2*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))**2)/
     &    (za(i5,i6)*zb(i2,i4)**3*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    t(i2,i3,i4)))-(2._dp*L0(-t(i2,i3,i4),-s(i3,i4))*za(i2,i4)*zb(i2,i3)*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))*
     &    (-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4)))/
     &    (s(i3,i4)*za(i5,i6)*zb(i2,i4)*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*t(i2,i3,i4))-
     &    (2._dp*L0(-t(i2,i3,i4),-s(i2,i3))*za(i2,i4)*zb(i2,i3)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))*
     &    (-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4)))/
     &    (s(i2,i3)*za(i5,i6)*zb(i2,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3
     &    ,i4))*
     &    t(i2,i3,i4))+Lsm1_2mh(s(i1,i2),t(i2,i3,i4),s(i3,i4),s(i5,i6))*
     &    (-(((za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2*
     &    (za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4))**3)/
     &    (za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**3*
     &    zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*t(i2,i3,i4))
     &    )-(za(i2,i4)**3*zb(i1,i6)**2)/
     &    (za(i2,i3)*za(i3,i4)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     &    zb(i5,i6)*t(i2,i3,i4))+
     &    (za(i1,i5)**2*zb(i2,i3)**2*(za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4
     &    ))*
     &    t(i2,i3,i4))/
     &    (za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**3*zb(i3,i
     &    4)*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))+
     &    (I3m(s(i2,i3),s(i1,i4),s(i5,i6))*za(i2,i4)*
     &    (-2._dp*za(i5,i6)*zb(i1,i6)*zb(i2,i3)*zb(i2,i6)*
     &    (-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i2))+
     &    2._dp*za(i1,i4)*za(i2,i3)*zb(i3,i4))-
     &    s(i5,i6)*zb(i1,i3)*((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*
     &    zb(i2,i6)-2._dp*za(i1,i4)*za(i3,i5)*zb(i2,i3)*zb(i4,i6))-
     &    2._dp*za(i3,i5)*zb(i2,i3)*zb(i3,i6)*
     &    (s(i5,i6)*(-s(i1,i4)-s(i2,i3)+s(i5,i6))+
     &    (-s(i1,i4)+s(i2,i3)-s(i5,i6))*t(i2,i3,i4))))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
     
          zaj_virt_a6_slc_fcc_pm=zaj_virt_a6_slc_fcc_pm+
     &    I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    ((za(i3,i5)*zb(i2,i3)*(s(i1,i4)*za(i1,i4)*zb(i1,i6)-
     &    s(i2,i4)*za(i2,i4)*zb(i2,i6))-
     &    5._dp*za(i1,i3)*za(i2,i4)*za(i3,i5)*zb(i1,i2)*zb(i2,i3)*zb(i3,i6)-
     &    (-s(i1,i3)+s(i2,i4))*za(i3,i4)*
     &    (-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))*zb(i3,i6)-
     &    za(i1,i4)*zb(i1,i2)*(-(s(i3,i4)*za(i1,i5)*zb(i1,i6))-
     &    2._dp*(-((s(i1,i3)+s(i2,i3))*za(i2,i5)*zb(i2,i6))-
     &    za(i2,i3)*za(i4,i5)*zb(i2,i6)*zb(i3,i4))-
     &    s(i2,i4)*za(i3,i5)*zb(i3,i6)+
     &    (2._dp*za(i2,i3)*zb(i1,i2)*zb(i3,i4)*
     &    (-(za(i2,i5)*zb(i2,i6))-za(i4,i5)*zb(i4,i6)))/zb(i1,i4)+
     &    3._dp*(za(i2,i5)*za(i3,i4)*zb(i2,i4)*zb(i3,i6)+
     &    s(i1,i3)*za(i4,i5)*zb(i4,i6))))/
     &    (2._dp*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))-
     &    (za(i2,i4)*zb(i1,i2)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &    zb(i3,i6))/
     &    (zb(i1,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*t(i1,i2,i4
     &    ))
     &    -(za(i2,i5)*zb(i1,i3)**2*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i
     &    6)))/
     &    (zb(i1,i4)*t(i1,i3,i4)**2)+
     &    (za(i2,i5)*zb(i1,i3)**2*(-(za(i2,i5)*zb(i2,i3))+
     &    za(i4,i5)*zb(i3,i4)))/
     &    (2._dp*za(i5,i6)*zb(i1,i4)*zb(i3,i4)*t(i1,i3,i4))+
     &    (za(i2,i4)*zb(i1,i3)*zb(i1,i6)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6)))/
     &    (2._dp*za(i3,i4)*zb(i1,i4)*zb(i5,i6)*t(i1,i3,i4))+
     &    ((za(i1,i5)**2*zb(i1,i3)**3*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))/
     &    (za(i5,i6)*zb(i1,i4)*zb(i3,i4))+
     &    (zb(i1,i2)*((-s(i1,i2)+2._dp*s(i1,i3)+s(i3,i4)-s(i5,i6))*
     &    za(i4,i5)-2._dp*za(i2,i4)*za(i5,i6)*zb(i2,i6))*zb(i3,i6))/
     &    zb(i1,i4)+(zb(i1,i3)*
     &    (-(za(i1,i3)*za(i1,i5)*zb(i1,i2)*zb(i1,i3))+
     &    za(i1,i3)*za(i2,i5)*zb(i1,i2)*zb(i2,i3)-
     &    za(i1,i5)*za(i3,i4)*zb(i1,i3)*zb(i2,i4)+
     &    za(i2,i5)*za(i3,i4)*zb(i2,i3)*zb(i2,i4))*zb(i3,i6))/
     &    (zb(i1,i4)*zb(i3,i4))+
     &    (za(i2,i4)**2*(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i2))+
     &    za(i1,i3)*zb(i1,i2)*zb(i1,i3)+
     &    za(i3,i4)*zb(i1,i3)*zb(i2,i4))*zb(i2,i6)**2)/
     &    (za(i3,i4)*zb(i1,i4)*zb(i5,i6))+
     &    (4._dp*za(i1,i4)*za(i2,i5)*zb(i1,i3)*
     &    (-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*zb(i2,i6))/
     &    t(i1,i3,i4)+((-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6))**2*
     &    (2._dp*s(i3,i4)*s(i5,i6)+(s(i1,i2)-s(i3,i4)-s(i5,i6))*t(i1,i3,i4))
     &    )/(za(i3,i4)*zb(i5,i6)*t(i1,i3,i4)**2)+
     &    (za(i4,i5)*zb(i1,i2)*(za(i1,i4)*zb(i1,i6)-za(i2,i4)*zb(i2,i6))*
     &    (s(i1,i4)+t(i2,i3,i4)))/(za(i3,i4)*zb(i1,i4)))/
     &    (2._dp*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))
     
          zaj_virt_a6_slc_fcc_pm=zaj_virt_a6_slc_fcc_pm+
     &    I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    (((s(i1,i3)-s(i2,i4))*za(i4,i5)*zb(i3,i4)*
     &    (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))+
     &    za(i1,i4)*(-(s(i1,i3)*za(i1,i5)*zb(i1,i3))+
     &    s(i2,i3)*za(i2,i5)*zb(i2,i3))*zb(i4,i6)+
     &    5._dp*za(i1,i2)*za(i1,i4)*za(i4,i5)*zb(i1,i3)*zb(i2,i4)*zb(i4,i6)+
     &    za(i1,i2)*zb(i2,i3)*(-(s(i3,i4)*za(i2,i5)*zb(i2,i6))+
     &    (2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i4)*
     &    (-(za(i1,i5)*zb(i1,i6))-za(i3,i5)*zb(i3,i6)))/za(i2,i3)+
     &    3._dp*(-(za(i1,i3)*za(i4,i5)*zb(i1,i6)*zb(i3,i4))+
     &    s(i2,i4)*za(i3,i5)*zb(i3,i6))-
     &    2._dp*(-((s(i1,i4)+s(i2,i4))*za(i1,i5)*zb(i1,i6))+
     &    za(i1,i5)*za(i3,i4)*zb(i1,i4)*zb(i3,i6))-
     &    s(i1,i3)*za(i4,i5)*zb(i4,i6)))/
     &    (2._dp*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     &    (za(i1,i2)*za(i4,i5)*zb(i1,i3)*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))/
     &    (za(i2,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*t(i1,i2,i3))-
     &    (za(i2,i4)**2*zb(i1,i6)*(-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4
     &    )))/
     &    (za(i2,i3)*t(i2,i3,i4)**2)+
     &    (za(i2,i4)*za(i2,i5)*zb(i1,i3)*
     &    (-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4)))/
     &    (2._dp*za(i2,i3)*za(i5,i6)*zb(i3,i4)*t(i2,i3,i4))+
     &    (za(i2,i4)**2*zb(i1,i6)*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6
     &    )))/
     &    (2._dp*za(i2,i3)*za(i3,i4)*zb(i5,i6)*t(i2,i3,i4))+
     &    ((za(i1,i5)**2*zb(i1,i3)**2*
     &    ((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)-
     &    za(i1,i2)*za(i2,i4)*zb(i2,i4)-za(i1,i3)*za(i2,i4)*zb(i3,i4)
     &    ))/(za(i2,i3)*za(i5,i6)*zb(i3,i4))-
     &    (za(i2,i4)*za(i4,i5)*(-(za(i1,i2)*za(i1,i4)*zb(i1,i6)*
     &    zb(i2,i4))+za(i1,i2)*za(i2,i4)*zb(i2,i4)*zb(i2,i6)-
     &    za(i1,i3)*za(i1,i4)*zb(i1,i6)*zb(i3,i4)+
     &    za(i1,i3)*za(i2,i4)*zb(i2,i6)*zb(i3,i4)))/
     &    (za(i2,i3)*za(i3,i4))+
     &    (za(i2,i4)**3*zb(i2,i6)**2*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))/
     &    (za(i2,i3)*za(i3,i4)*zb(i5,i6))-
     &    (za(i1,i2)*za(i4,i5)*((-s(i1,i2)+2._dp*s(i2,i4)+s(i3,i4)-s(i5,i6))*
     &    zb(i3,i6)+2._dp*za(i1,i5)*zb(i1,i3)*zb(i5,i6)))/za(i2,i3)+
     &    (za(i1,i2)*(-(za(i1,i5)*zb(i1,i3))+za(i2,i5)*zb(i2,i3))*
     &    zb(i3,i6)*(s(i2,i3)+t(i1,i3,i4)))/(za(i2,i3)*zb(i3,i4))+
     &    (4._dp*za(i1,i5)*za(i2,i4)*zb(i1,i6)*zb(i2,i3)*
     &    (za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4)))/t(i2,i3,i4)+
     &    ((za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4))*
     &    (-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4))**2*
     &    (2._dp*s(i3,i4)*s(i5,i6)+(s(i1,i2)-s(i3,i4)-s(i5,i6))*t(i2,i3,i4)))/
     &    (za(i5,i6)*zb(i3,i4)*t(i2,i3,i4)**2))/
     &    (2._dp*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))

      end function

      function zaj_virt_a6_slc_fcc_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_fcc_mp

        complex(dp)  :: L0,Lsm1,Lsm1_2mh,Lsm1_2mht,I3m

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)


        zaj_virt_a6_slc_fcc_mp=
     &    ((Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))+
     &    Lsm1_2mht(s(i1,i4),t(i1,i2,i3),s(i2,i3),s(i5,i6)))*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))**2)/
     &    (za(i5,i6)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i2,i3)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (2._dp*L0(-t(i2,i3,i4),-s(i5,i6))*za(i1,i5)*zb(i1,i4)*zb(i2,i4)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    (s(i5,i6)*za(i5,i6)*zb(i2,i3)*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i3,i4))-
     &    (I3m(s(i1,i4),s(i2,i3),s(i5,i6))*zb(i1,i4)*
     &    (-(s(i1,i4)*s(i2,i3)*za(i2,i5)**2)+
     &    (-(za(i1,i2)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3)))+
     &    za(i2,i4)*(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))**2))/
     &    (2._dp*za(i1,i4)*za(i5,i6)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i
     &    4))*
     &    zb(i2,i3)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))-
     &    ((Lsm1(-s(i1,i2),-t(i1,i2,i4),-s(i1,i4),-t(i1,i2,i4))+
     &    Lsm1_2mht(s(i2,i3),t(i1,i2,i4),s(i1,i4),s(i5,i6)))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))**2)/
     &    (za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i5,i6))-
     &    (2._dp*L0(-t(i1,i3,i4),-s(i5,i6))*za(i1,i3)*za(i2,i3)*zb(i2,i6)*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i5,i6)*za(i1,i4)*za(i3,i4)*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i5,i6))+
     &    (I3m(s(i2,i3),s(i1,i4),s(i5,i6))*za(i2,i3)*
     &    (-(s(i1,i4)*s(i2,i3)*zb(i1,i6)**2)+
     &    (zb(i1,i2)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))+
     &    zb(i1,i3)*(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))**2))/
     &    (2._dp*za(i1,i4)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i2,i3)*zb(i5,i6))+
     &    Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*
     &    (-(((za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))**2)/
     &    (za(i5,i6)*zb(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-
     &    za(i2,i4)*zb(i2,i3))**3))+
     &    (za(i4,i5)**2*zb(i1,i3)**2*t(i1,i2,i3)**2)/
     &    (za(i5,i6)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i2,i3)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**3))
     
          zaj_virt_a6_slc_fcc_mp=zaj_virt_a6_slc_fcc_mp+
     &    Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),s(i5,i6))*
     &    (((-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))**2*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))/
     &    (za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**3*
     &    zb(i5,i6))-(za(i2,i4)**2*zb(i3,i6)**2*t(i1,i2,i4)**2)/
     &    (za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**3*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i5,i6)))+
     &    (Lsm1(-s(i1,i4),-t(i1,i3,i4),-s(i3,i4),-t(i1,i3,i4))*za(i2,i5)**2*
     &    zb(i1,i4)**2)/
     &    (za(i5,i6)*zb(i1,i3)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    t(i1,i3,i4))-(2._dp*L0(-t(i1,i3,i4),-s(i3,i4))*za(i1,i3)*zb(i1,i4)*
     &    (za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i3,i4)*za(i1,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    zb(i5,i6)*t(i1,i3,i4))+Lsm1_2mh(s(i1,i2),t(i1,i3,i4),s(i3,i4),
     &    s(i5,i6))*(-((za(i2,i5)**2*zb(i1,i4)**2*
     &    (-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4)))/
     &    (za(i5,i6)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     &    zb(i3,i4)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    t(i1,i3,i4)))-(za(i1,i3)*
     &    (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))**2*
     &    (za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))**2)/
     &    (za(i1,i4)*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    *
     &    3._dp*zb(i5,i6)*t(i1,i3,i4))+
     &    (za(i1,i3)**3*zb(i2,i6)**2*t(i1,i3,i4))/
     &    (za(i1,i4)*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    *3*
     &    zb(i5,i6)))
     
          zaj_virt_a6_slc_fcc_mp=zaj_virt_a6_slc_fcc_mp+
     &    I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    ((za(i1,i2)*za(i3,i5)*zb(i1,i4)*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    (za(i1,i4)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*t(i1,i2,i4
     &    ))
     &    +(4._dp*za(i1,i2)*za(i2,i4)*za(i3,i5)*zb(i1,i4)*zb(i2,i6)*zb(i3,i4)-
     &    s(i2,i3)*za(i1,i5)*zb(i1,i6)*
     &    (-3._dp*za(i1,i2)*zb(i1,i4)+za(i2,i3)*zb(i3,i4))-
     &    s(i3,i4)*za(i1,i2)*za(i3,i5)*zb(i1,i4)*zb(i3,i6)+
     &    s(i2,i4)*za(i2,i3)*(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))*
     &    zb(i3,i6)+(za(i1,i5)*za(i2,i4)*zb(i1,i4)-
     &    za(i2,i4)*za(i2,i5)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i4)*zb(i3,i6))-
     &    za(i2,i3)*zb(i2,i4)*zb(i3,i6))-
     &    (za(i1,i2)*za(i3,i4)*zb(i1,i4)*zb(i2,i3)*
     &    (-2._dp*za(i1,i2)*za(i3,i5)*zb(i3,i6)+
     &    za(i2,i5)*(za(i1,i2)*zb(i2,i6)-za(i1,i3)*zb(i3,i6))))/
     &    za(i1,i4)-za(i2,i4)*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i4))-za(i2,i3)*zb(i3,i4))*zb(i4,i6)-
     &    za(i4,i5)*(-(s(i2,i4)*za(i1,i2)*zb(i1,i4))+
     &    s(i1,i2)*za(i2,i3)*zb(i3,i4))*zb(i4,i6)-
     &    s(i1,i4)*za(i3,i5)*zb(i3,i4)*
     &    (za(i2,i3)*zb(i3,i6)+za(i2,i4)*zb(i4,i6))+
     &    (-(za(i1,i3)*za(i2,i3)*zb(i1,i3))-
     &    za(i1,i3)*za(i2,i4)*zb(i1,i4))*
     &    (-(za(i1,i5)*zb(i1,i4)*zb(i3,i6))-
     &    za(i2,i5)*zb(i2,i3)*zb(i4,i6))-
     &    zb(i1,i4)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    (za(i1,i3)*za(i2,i5)*zb(i2,i6)+
     &    2._dp*za(i1,i5)*(-(za(i1,i3)*zb(i1,i6))-za(i3,i4)*zb(i4,i6)))+
     &    2._dp*za(i2,i3)*za(i2,i5)*zb(i3,i4)*
     &    (-(za(i1,i4)*zb(i1,i6)*zb(i2,i4))-
     &    zb(i2,i6)*(s(i1,i2)+t(i1,i2,i4))))/
     &    (2._dp*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))+
     &    (za(i2,i5)*zb(i1,i4)*(-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    ((-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*t(i1,i3,i4)**2)-
     &    (2._dp*za(i1,i3)**2*za(i2,i5)*zb(i1,i4)*zb(i2,i6))/
     &    (za(i1,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*t(i1,i3,i4
     &    ))
     &    +(za(i1,i3)*za(i2,i3)*za(i2,i5)*zb(i1,i4)*zb(i1,i6))/
     &    (2._dp*za(i3,i4)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    t(i1,i3,i4))+(za(i2,i5)*zb(i1,i4)*
     &    (-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    (2._dp*za(i5,i6)*zb(i3,i4)*(-(za(i1,i2)*zb(i1,i3))-
     &    za(i2,i4)*zb(i3,i4))*t(i1,i3,i4))+
     &    (za(i1,i3)*(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6))**2*
     &    (2._dp*s(i3,i4)*s(i5,i6)+(s(i1,i2)-s(i3,i4)-s(i5,i6))*t(i1,i3,i4)))/
     &    (2._dp*za(i1,i4)*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-
     &    za(i1,i4)*zb(i2,i4))*zb(i5,i6)*t(i1,i3,i4)**2)+
     &    (-((za(i1,i3)*za(i2,i3)*za(i2,i5)*zb(i1,i6)*zb(i2,i4))/za(i3,i4))+
     &    za(i2,i3)*(-(za(i1,i5)*zb(i1,i4))-2._dp*za(i3,i5)*zb(i3,i4))*
     &    zb(i4,i6)+(2._dp*za(i1,i2)*za(i5,i6)*zb(i4,i6)*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/za(i1,i4)-
     &    (za(i2,i5)*zb(i1,i4)*
     &    (-(za(i1,i5)*zb(i1,i4)*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))-
     &    za(i1,i2)*za(i5,i6)*zb(i2,i4)*zb(i4,i6)))/
     &    (za(i5,i6)*zb(i3,i4))+
     &    ((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i3)*za(i2,i3)*zb(i3,i6)*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (za(i1,i4)*za(i3,i4)*zb(i5,i6))-
     &    (za(i1,i3)*za(i2,i5)*(-(za(i1,i3)*zb(i1,i6)*
     &    (s(i5,i6)-t(i1,i2,i4)))+
     &    za(i3,i4)*zb(i4,i6)*(s(i5,i6)+t(i1,i3,i4))))/
     &    (za(i1,i4)*za(i3,i4)))/
     &    (2._dp*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))))
     
          zaj_virt_a6_slc_fcc_mp=zaj_virt_a6_slc_fcc_mp+
     &    (2._dp*L0(-t(i2,i3,i4),-s(i3,i4))*za(i2,i3)*zb(i2,i4)*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    (s(i3,i4)*za(i5,i6)*zb(i2,i3)*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*t(i2,i3,i4))-
     &    (Lsm1(-s(i2,i3),-t(i2,i3,i4),-s(i3,i4),-t(i2,i3,i4))*za(i2,i3)**2*
     &    zb(i1,i6)**2)/
     &    (za(i2,i4)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i5,i6)*
     &    t(i2,i3,i4))+Lsm1_2mh(s(i1,i2),t(i2,i3,i4),s(i3,i4),s(i5,i6))*
     &    (-((zb(i2,i4)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**2)/
     &    (za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-
     &    za(i1,i4)*zb(i2,i4))**3*zb(i3,i4)*t(i2,i3,i4)))-
     &    (za(i2,i3)**2*(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*zb(i1,i6)*
     &    *2)/
     &    (za(i3,i4)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i5,i6)*t(i2,i3,i4)
     &    )+(za(i1,i5)**2*zb(i2,i4)**3*t(i2,i3,i4))/
     &    (za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    *3*
     &    zb(i3,i4)))
     
          zaj_virt_a6_slc_fcc_mp=zaj_virt_a6_slc_fcc_mp+
     &    I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    (-((za(i2,i3)*zb(i1,i2)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &    zb(i4,i6))/
     &    ((za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i2,i3)*t(i1,i2,i3)))
     &    +(-(s(i1,i4)*za(i2,i5)*(3._dp*za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4
     &    ))*
     &    zb(i2,i6))+s(i1,i3)*za(i4,i5)*zb(i1,i4)*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))+
     &    (-(za(i2,i3)*zb(i1,i3)*zb(i2,i4))-
     &    za(i2,i4)*zb(i1,i4)*zb(i2,i4))*
     &    (-(za(i1,i4)*za(i3,i5)*zb(i1,i6))-
     &    za(i2,i3)*za(i4,i5)*zb(i2,i6))+
     &    (-(za(i1,i3)*za(i4,i5)*zb(i1,i4))-
     &    za(i2,i3)*za(i4,i5)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i3)*zb(i1,i6))+
     &    za(i2,i3)*zb(i1,i3)*zb(i2,i6))-
     &    za(i2,i3)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (za(i1,i5)*zb(i1,i6)*zb(i2,i4)+
     &    2._dp*zb(i2,i6)*(-(za(i2,i5)*zb(i2,i4))+za(i3,i5)*zb(i3,i4)))-
     &    za(i3,i5)*(s(i1,i3)*za(i2,i3)*zb(i1,i2)-
     &    s(i1,i2)*za(i3,i4)*zb(i1,i4))*zb(i3,i6)+
     &    s(i3,i4)*za(i2,i3)*za(i4,i5)*zb(i1,i2)*zb(i4,i6)+
     &    4._dp*za(i1,i5)*za(i2,i3)*za(i3,i4)*zb(i1,i2)*zb(i1,i3)*zb(i4,i6)+
     &    s(i2,i3)*za(i3,i4)*(za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4))*
     &    zb(i4,i6)-za(i3,i5)*zb(i1,i3)*
     &    (za(i2,i3)*zb(i1,i2)+za(i3,i4)*zb(i1,i4))*
     &    (za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))-
     &    (za(i1,i4)*za(i2,i3)*zb(i1,i2)*zb(i3,i4)*
     &    (zb(i1,i6)*(-(za(i1,i5)*zb(i1,i2))-za(i4,i5)*zb(i2,i4))+
     &    2._dp*za(i4,i5)*zb(i1,i2)*zb(i4,i6)))/zb(i2,i3)-
     &    2._dp*za(i3,i4)*zb(i1,i4)*zb(i1,i6)*
     &    (-(za(i1,i3)*za(i2,i5)*zb(i2,i3))-
     &    za(i1,i5)*(s(i1,i2)+t(i1,i2,i3))))/
     &    (2._dp*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     &    (za(i2,i3)*(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*zb(i1,i6)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    ((za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*t(i2,i3,i4)**2)-
     &    (2._dp*za(i1,i5)*za(i2,i3)*zb(i1,i6)*zb(i2,i4)**2)/
     &    (zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*t(i2,i3,i4
     &    ))
     &    -(za(i2,i3)*za(i2,i5)*zb(i1,i4)*zb(i1,i6)*zb(i2,i4))/
     &    (2._dp*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i3,i4)*t(i2,i3,i
     &    4))+
     &    (za(i2,i3)*(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*zb(i1,i6)*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (2._dp*za(i3,i4)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i5,i6)
     &    *
     &    t(i2,i3,i4))+(zb(i2,i4)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))**2*
     &    (2._dp*s(i3,i4)*s(i5,i6)+(s(i1,i2)-s(i3,i4)-s(i5,i6))*t(i2,i3,i4)))/
     &    (2._dp*za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i
     &    4))*
     &    zb(i3,i4)*t(i2,i3,i4)**2)+
     &    ((za(i1,i3)*za(i2,i5)*zb(i1,i4)*zb(i1,i6)*zb(i2,i4))/zb(i3,i4)+
     &    ((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i4,i5)*zb(i1,i4)*zb(i2,i4)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    (za(i5,i6)*zb(i2,i3)*zb(i3,i4))+
     &    za(i3,i5)*zb(i1,i4)*(-(za(i2,i3)*zb(i2,i6))+
     &    2._dp*za(i3,i4)*zb(i4,i6))+
     &    (2._dp*za(i3,i5)*zb(i1,i2)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))*zb(i5,i6))/
     &    zb(i2,i3)-(za(i2,i3)*zb(i1,i6)*
     &    (-(za(i2,i3)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &    zb(i2,i6))-za(i1,i3)*za(i3,i5)*zb(i1,i2)*zb(i5,i6)))/
     &    (za(i3,i4)*zb(i5,i6))+
     &    (zb(i1,i6)*zb(i2,i4)*(-(za(i2,i5)*zb(i2,i4)*
     &    (s(i5,i6)-t(i1,i2,i3)))-
     &    za(i3,i5)*zb(i3,i4)*(s(i5,i6)+t(i2,i3,i4))))/
     &    (zb(i2,i3)*zb(i3,i4)))/
     &    (2._dp*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))

      end function

      function zaj_virt_a6_slc_fsc_pp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_fsc_pp

        complex(dp)  :: L0,L1,Lsm1,Lsm1_2me,Lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_virt_a6_slc_fsc_pp =
     &    -((za(i1,i5)*za(i2,i5)*((Lnrat(-s(i2,i3),-s(i5,i6))*za(i1,i2))/
     &    (za(i1,i4)*za(i2,i3))+Lnrat(-s(i1,i2),-s(i5,i6))/za(i3,i4)))/
     &    (za(i1,i3)*za(i1,i4)*za(i5,i6)))+
     &    (Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))*za(i1,i5)**2*
     &    za(i2,i3))/
     &    (za(i1,i3)**2*za(i1,i4)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1_2me(t(i1,i2,i3),t(i2,i3,i4),s(i2,i3),s(i5,i6))*za(i1,i5)**2*
     &    za(i2,i4)**2)/(za(i1,i4)**3*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lnrat(-t(i1,i2,i3),-s(i5,i6))*za(i2,i5)**2)/
     &    (za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lnrat(-t(i2,i3,i4),-s(i5,i6))*za(i2,i5)**2)/
     &    (2._dp*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1_2me(t(i1,i2,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*za(i2,i3)*
     &    za(i4,i5)**2)/(za(i1,i4)*za(i3,i4)**3*za(i5,i6))-
     &    (L0(-t(i1,i2,i3),-s(i2,i3))*za(i1,i2)**2*
     &    (za(i1,i5)*za(i2,i5)*zb(i1,i2)+za(i1,i5)*za(i3,i5)*zb(i1,i3)))/
     &    (s(i2,i3)*za(i1,i3)*za(i1,i4)**2*za(i2,i3)*za(i5,i6))+
     &    (L1(-t(i2,i3,i4),-s(i5,i6))*za(i1,i2)**2*za(i5,i6)*zb(i1,i6)**2)/
     &    (2._dp*s(i5,i6)**2*za(i1,i4)*za(i2,i3)*za(i3,i4))+
     &    (L0(-t(i2,i3,i4),-s(i5,i6))*((za(i1,i5)**2*za(i2,i4)*
     &    (-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4)))/
     &    (za(i1,i4)**2*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (za(i1,i2)*za(i2,i5)*zb(i1,i6))/(za(i1,i4)*za(i2,i3)*za(i3,i4))))/
     &    s(i5,i6)+(L0(-t(i1,i2,i3),-s(i1,i2))*za(i2,i3)*
     &    (-(za(i1,i5)*za(i3,i5)*zb(i1,i3))-za(i2,i5)*za(i3,i5)*zb(i2,i3)))/
     &    (s(i1,i2)*za(i1,i3)*za(i3,i4)**2*za(i5,i6))-
     &    (L0(-t(i1,i2,i4),-s(i1,i2))*za(i2,i4)*za(i4,i5)*
     &    (-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))/
     &    (s(i1,i2)*za(i1,i4)*za(i3,i4)**2*za(i5,i6))+
     &    (L0(-t(i2,i3,i4),-s(i2,i3))*(za(i1,i5)*za(i2,i4)+za(i1,i4)*za(i2,i
     &    5))*
     &    za(i4,i5)*zb(i3,i4))/(s(i2,i3)*za(i1,i4)**2*za(i3,i4)*za(i5,i6))-
     &    (L0(-t(i1,i2,i4),-s(i5,i6))*za(i2,i4)*za(i3,i5)*zb(i3,i6))/
     &    (s(i5,i6)*za(i1,i4)*za(i3,i4)**2)+
     &    (L1(-t(i1,i2,i4),-s(i5,i6))*za(i2,i3)*za(i5,i6)*zb(i3,i6)**2)/
     &    (s(i5,i6)**2*za(i1,i4)*za(i3,i4))+
     &    (L1(-t(i1,i2,i3),-s(i5,i6))*za(i2,i4)**2*za(i5,i6)*zb(i4,i6)**2)/
     &    (s(i5,i6)**2*za(i1,i4)*za(i2,i3)*za(i3,i4))+
     &    (L0(-t(i1,i2,i3),-s(i5,i6))*(-((za(i2,i5)*za(i4,i5)*
     &    (-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4)))/
     &    (za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6)))-
     &    (za(i2,i4)*(-(za(i1,i5)*za(i2,i4))+2._dp*za(i1,i4)*za(i2,i5))*
     &    zb(i4,i6))/(za(i1,i4)**2*za(i2,i3)*za(i3,i4))+
     &    (za(i2,i4)*za(i4,i5)*zb(i4,i6))/(za(i1,i4)*za(i3,i4)**2)))/s(i5,i6
     &    )+
     &    (-(za(i2,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6)))+
     &    (za(i2,i5)*za(i4,i5)*zb(i3,i4))/
     &    (za(i1,i4)*za(i3,i4)*za(i5,i6)*t(i2,i3,i4))-
     &    (za(i1,i5)*za(i2,i5)*zb(i1,i3)*zb(i3,i4))/
     &    (za(i1,i4)*za(i5,i6)*t(i1,i3,i4)*t(i2,i3,i4))-
     &    (zb(i3,i4)*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6))*
     &    (za(i2,i3)*zb(i3,i6)+za(i2,i4)*zb(i4,i6)))/
     &    (za(i1,i4)*za(i3,i4)*zb(i5,i6)*t(i1,i3,i4)*t(i2,i3,i4)))/2-
     &    (L1(-s(i2,i3),-t(i2,i3,i4))*za(i2,i3)*za(i4,i5)**2*zb(i3,i4)**2)/
     &    (2._dp*za(i1,i4)*za(i3,i4)*za(i5,i6)*t(i2,i3,i4)**2)

      end function

      ! what an awful amplitude ..
      function zaj_virt_a6_slc_fsc_pm(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_slc_fsc_pm

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        complex(dp) :: L0,L1,Lsm1,Lsm1_2mh,I3m,Lnrat,Ls1
        complex(dp) :: tmp

      zaj_virt_a6_slc_fsc_pm= 
     .-((Lsm1(-s(i1,i2),-t(i1,i2,i4),-s(i1,i4),-t(i1,i2,i4))*zb(i1,i4)*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))**2)/
     .(za(i5,i6)*zb(i2,i4)**2*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4
     .))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))-
     .(za(i2,i4)*(-((L0(-t(i1,i2,i4),-s(i1,i2))*zb(i1,i2)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))**2)/
     .(s(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))-
     .(L0(-t(i1,i2,i4),-s(i1,i4))*zb(i1,i4)*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))**2)/
     .(s(i1,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))))/
     .(za(i5,i6)*zb(i2,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))
     .+
     .(za(i1,i5)*za(i4,i5)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))
     ./
     .(2._dp*za(i3,i4)*za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i
     .4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     .(za(i1,i5)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .((za(i1,i4)*((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*zb(i1,i2)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6)))/za(i3,i4)
     .-((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)+
     .2._dp*za(i1,i2)*(za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4)))*
     .zb(i2,i3)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))-
     .(Lnrat(-s(i1,i2),-s(i5,i6))*za(i3,i5)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4)))/
     .(2._dp*s(i3,i4)*za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4
     .))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     .(za(i3,i5)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .((za(i3,i5)*zb(i1,i2)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i4)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i3)))/za(i3,i4)-
     .((-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i3)+
     .2._dp*za(i2,i4)*zb(i1,i2)*zb(i3,i4)))/zb(i3,i4)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .((-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i5)*za(i4,i5))/
     .(za(i2,i3)*za(i5,i6))+
     .((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i5)*zb(i1,i6))/
     .(2._dp*za(i2,i3)*zb(i1,i4))+za(i4,i5)*zb(i3,i6)))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s
     .(i5,i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     .((-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i5)*zb(i1,i6))/
     .(2._dp*za(i2,i3)*zb(i1,i4))+za(i4,i5)*zb(i3,i6)-
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i1,i6)*zb(i3,i6))/
     .(zb(i1,i4)*zb(i5,i6))))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s
     .(i5,i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     .((za(i2,i4)**2*za(i3,i5)**2*zb(i1,i2))/(za(i2,i3)*za(i3,i4)*za(i5,
     .i6))-
     .(za(i3,i5)*zb(i1,i3)*zb(i3,i6))/zb(i3,i4)+
     .((za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*zb(i1,i6)*zb(i3,i6))/
     .(zb(i1,i4)*zb(i5,i6)))/
     .(2._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .(-((za(i2,i5)*za(i4,i5)*(-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4
     .)))/
     .(za(i2,i3)*za(i5,i6)))+
     .(za(i2,i4)*za(i4,i5)*zb(i4,i6))/za(i3,i4)-
     .(za(i1,i2)*zb(i1,i3)**2*zb(i4,i6)**2)/(zb(i1,i4)*zb(i3,i4)*zb(i5,i
     .6)))/
     .(2._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     .((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i2,i6)*zb(i3,i6))
     ./
     .(2._dp*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i3,i4)*zb(i5,i6))+
     .(Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))*za(i2,i3)*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     .(za(i1,i3)**2*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))+
     .(zb(i1,i3)*((L0(-t(i1,i2,i3),-s(i1,i2))*za(i1,i2)*
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))**2)/
     .(s(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     .(L0(-t(i1,i2,i3),-s(i2,i3))*za(i2,i3)*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     .(s(i2,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))))/
     .(za(i1,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*zb(i5,i6))

      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+
     .(Lnrat(-s(i1,i2),-s(i5,i6))*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i
     .2,i3))*
     .(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6))*zb(i4,i6))/
     .(2._dp*s(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*zb(i5,i6))+
     .((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i4,i6)*
     .((((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i4)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)))/za(i3,i4)+
     .(za(i1,i2)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i3)+
     .2._dp*za(i2,i4)*zb(i1,i2)*zb(i3,i4))*zb(i4,i6))/zb(i3,i4)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))-
     .((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i2,i6)*
     .(-(za(i1,i4)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i6)-
     .2._dp*zb(i1,i2)*(za(i2,i3)*zb(i3,i6)+za(i2,i4)*zb(i4,i6))))-
     .(zb(i2,i3)*(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6))-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6)))/zb(i3,i4))
     .)/(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)
     .*s(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6))+
     .I3m(s(i1,i4),s(i2,i3),s(i5,i6))*
     .(-(((-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(-(za(i1,i4)*za(i2,i5)*zb(i1,i6)*zb(i2,i3))+
     .(za(i4,i5)*(-2._dp*za(i2,i5)*zb(i2,i3)-za(i5,i6)*zb(i3,i6))*
     .zb(i5,i6))/2))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     .i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))+
     .(3._dp*za(i1,i4)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i2,
     .i3)*
     .(((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i5)*za(i5,i6)*zb(i1,i6)*
     .zb(i5,i6))/2+za(i4,i5)*zb(i1,i4)*
     .(-(s(i5,i6)*za(i2,i3)*zb(i3,i6))+
     .(s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     .i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm+
     .I3m(s(i2,i3),s(i1,i4),s(i5,i6))*
     .((3._dp*za(i1,i4)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i2
     .,i3)*
     .(za(i2,i3)*(-(s(i5,i6)*za(i4,i5)*zb(i1,i4))-
     .(-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i5,i6)*zb(i1,i6))*
     .zb(i3,i6)+((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i5)*
     .za(i5,i6)*zb(i1,i6)*zb(i5,i6))/2))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     .i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     .((-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(-(za(i1,i4)*za(i2,i5)*zb(i1,i6)*zb(i2,i3))-
     .(za(i5,i6)*zb(i3,i6)*
     .(-2._dp*za(i1,i4)*zb(i1,i6)+za(i4,i5)*zb(i5,i6)))/2))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     .i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))+
     .Lnrat(-s(i1,i4),-s(i5,i6))*((4._dp*za(i2,i4)*za(i2,i5)*zb(i1,i6)*zb(
     .i2,i3)+
     .(2._dp*za(i4,i5)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i5)+
     .2._dp*za(i2,i3)*za(i5,i6)*zb(i3,i6)))/(za(i2,i3)*za(i5,i6))+
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i4)+
     .za(i1,i4)*za(i2,i5)*zb(i1,i5)+za(i1,i4)*za(i2,i6)*zb(i1,i6))*
     .(-((za(i2,i5)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3)))/
     .(za(i2,i3)*za(i5,i6)))+(zb(i1,i6)*zb(i3,i6))/zb(i5,i6)))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s
     .(i5,i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     .(3._dp*za(i1,i4)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i4,i5)*zb(i1,i4)+
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*zb(i1,i6))*zb(i3,i6)
     .-za(i2,i5)*zb(i2,i3)*((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i1,i6)-
     .2._dp*za(i4,i5)*zb(i1,i4)*zb(i5,i6))))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     .i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))+
     .Lnrat(-s(i2,i3),-s(i5,i6))*((3._dp*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)
     .*zb(i1,i4))*
     .zb(i2,i3)*(-(za(i1,i4)*zb(i1,i6)*
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i5)+
     .2._dp*za(i2,i3)*za(i5,i6)*zb(i3,i6)))+
     .za(i4,i5)*((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*zb(i3,i6)-
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     .i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     .(4._dp*za(i1,i4)*za(i2,i5)*zb(i1,i3)*zb(i1,i6)+
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i1,i3)+
     .za(i2,i5)*zb(i1,i5)*zb(i2,i3)+za(i2,i6)*zb(i1,i6)*zb(i2,i3))*
     .(-((za(i2,i5)*za(i4,i5))/za(i5,i6))+
     .(zb(i1,i6)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     .(zb(i1,i4)*zb(i5,i6)))-
     .(2._dp*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i3,i6)*
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i1,i6)-
     .2._dp*za(i4,i5)*zb(i1,i4)*zb(i5,i6)))/(zb(i1,i4)*zb(i5,i6)))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s
     .(i5,i6)-
     .2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))

      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+(Lnrat(-s(i1,i2),-s(i5,i6))*za(i1,i2)*
     .((-6._dp*zb(i1,i2)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i3)-
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i4))*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2-
     .((-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*
     .(zb(i1,i6)*zb(i3,i4)-zb(i1,i3)*zb(i4,i6)))/
     .((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*zb(i3,i4)*
     .zb(i5,i6))+(-((za(i4,i5)*
     .(2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i5)*
     .zb(i1,i3)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))
     .+za(i4,i5)*(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*zb(i3,i4)))/
     .za(i5,i6))-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i3)*
     .zb(i1,i2)*(-(za(i3,i5)*zb(i3,i6))+za(i4,i5)*zb(i4,i6)))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))-
     .(za(i2,i4)*zb(i1,i2)*
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*zb(i3,i6))/
     .zb(i5,i6)-(zb(i1,i2)*zb(i4,i6)*
     .(-((2._dp*(s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)+
     .za(i1,i2)*za(i3,i4)*zb(i1,i3))*zb(i3,i6))+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .zb(i4,i6))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))/
     .(zb(i3,i4)*zb(i5,i6))+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i4)*zb(i3,i6)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6)-
     .za(i4,i5)*zb(i5,i6)))/(zb(i3,i4)*zb(i5,i6))+
     .za(i4,i5)*(-(za(i2,i3)*zb(i1,i2)*zb(i3,i6))+
     .za(i4,i5)*zb(i1,i4)*zb(i5,i6))+
     .4._dp*za(i3,i5)*zb(i1,i2)*
     .(za(i2,i4)*zb(i3,i6)-
     .(za(i2,i3)*za(i4,i5)*zb(i3,i4)*zb(i5,i6))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*za(i1,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .-(Lnrat(-s(i1,i2),-s(i5,i6))*zb(i1,i2)*
     .(-(((-(za(i2,i5)*za(i3,i4))-za(i2,i4)*za(i3,i5))*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))/
     .(za(i3,i4)*za(i5,i6)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))+
     .(6._dp*za(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i4)-
     .2._dp*za(i2,i3)*zb(i1,i2)*zb(i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2+
     .(-((za(i1,i2)*za(i4,i5)*zb(i1,i3)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))/za(i5,i6))
     .+(za(i1,i2)*za(i3,i5)*(((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i3,i5)*
     .zb(i1,i4)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))-
     .za(i4,i5)*(2._dp*(s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i3)+
     .za(i2,i4)*zb(i1,i2)*zb(i3,i4))))/(za(i3,i4)*za(i5,i6))
     .+((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i3)*za(i4,i5)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)+
     .za(i5,i6)*zb(i3,i6)))/(za(i3,i4)*za(i5,i6))+
     .zb(i3,i6)*(za(i1,i2)*za(i4,i5)*zb(i1,i4)-
     .za(i2,i3)*za(i5,i6)*zb(i3,i6))-
     .4._dp*za(i1,i2)*(za(i4,i5)*zb(i1,i3)-
     .(za(i3,i4)*za(i5,i6)*zb(i1,i4)*zb(i3,i6))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))*zb(i4,i6)+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i4)*
     .(za(i3,i5)*zb(i3,i6)-za(i4,i5)*zb(i4,i6)))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))+
     .(zb(i3,i6)*(-(za(i3,i4)*
     .(-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))*
     .zb(i3,i6))+
     .2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i4,i6)-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .zb(i5,i6))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*zb(i2,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm+
     .(Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*za(i2,i3)*zb(i4,
     .i6)**2*
     .t(i1,i2,i3)**2)/
     .((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**3*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))+
     .(L0(-t(i1,i2,i3),-s(i5,i6))*za(i2,i4)*zb(i4,i6)**2*t(i1,i2,i3)**2)
     ./
     .(s(i5,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))+
     .(Lnrat(-s(i3,i4),-s(i5,i6))*((-6._dp*za(i1,i2)*za(i3,i4)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i4)-
     .2._dp*za(i2,i3)*zb(i1,i2)*zb(i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     .(za(i4,i5)*(za(i1,i2)*za(i3,i4)*zb(i1,i4)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))+
     .3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i3)*za(i4,i5)*
     .zb(i3,i4)-za(i2,i3)*za(i3,i4)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))*zb(i3,i4)-
     .2._dp*za(i1,i2)*zb(i1,i3)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i3,i5)+
     .2._dp*za(i3,i4)*za(i5,i6)*zb(i4,i6))))/za(i5,i6)+
     .(za(i3,i4)*zb(i3,i6)*((-(za(i1,i4)*zb(i1,i6))-
     .za(i2,i4)*zb(i2,i6))*
     .(-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))+
     .2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i4,i6)-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .zb(i5,i6)+(za(i2,i3)*
     .((za(i3,i5)*za(i4,i5)*zb(i3,i4))/za(i5,i6)+
     .(za(i3,i4)*zb(i3,i6)*zb(i4,i6))/zb(i5,i6))*
     .(-s(i1,i2)**2+2._dp*s(i1,i2)*s(i3,i4)-s(i3,i4)**2+
     .2._dp*s(i1,i2)*s(i5,i6)+2._dp*s(i3,i4)*s(i5,i6)-s(i5,i6)**2+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*(t(i1,i2,i3)-t(i1,i2,i4))))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i3)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     .(Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),s(i5,i6))*za(i3,i5)**2*zb(
     .i1,i4)*
     .t(i1,i2,i4)**2)/
     .(za(i5,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**3*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))-
     .(L0(-t(i1,i2,i4),-s(i5,i6))*za(i3,i5)**2*zb(i1,i3)*t(i1,i2,i4)**2)
     ./
     .(s(i5,i6)*za(i5,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**
     .2*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .(Lnrat(-s(i3,i4),-s(i5,i6))*((za(i4,i5)*
     .(2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i5)*zb(i1,i3)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))+
     .(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)))*zb(i3,i4))/
     .za(i5,i6)-(6._dp*zb(i1,i2)*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i3)-
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i4))*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*zb(i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     .(zb(i3,i6)*(za(i2,i3)*zb(i1,i2)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*zb(i3,i4)-
     .za(i3,i4)*zb(i1,i4)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*zb(i3,i4)-
     .3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*zb(i1,i4)*
     .zb(i3,i6)+2._dp*za(i2,i4)*zb(i1,i2)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i4,i6)+
     .2._dp*za(i3,i5)*zb(i3,i4)*zb(i5,i6))))/zb(i5,i6)+
     .(zb(i1,i4)*((za(i3,i5)*za(i4,i5)*zb(i3,i4))/za(i5,i6)+
     .(za(i3,i4)*zb(i3,i6)*zb(i4,i6))/zb(i5,i6))*
     .(-s(i1,i2)**2+2._dp*s(i1,i2)*s(i3,i4)-s(i3,i4)**2+
     .2._dp*s(i1,i2)*s(i5,i6)+2._dp*s(i3,i4)*s(i5,i6)-s(i5,i6)**2+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*(-t(i1,i2,i3)+t(i1,i2,i4))))/
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i4)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     .(Ls1(-s(i1,i4),-t(i1,i3,i4),-s(i3,i4),-t(i1,i3,i4))*za(i1,i4)**2*z
     .b(i1,i3)**2*
     .(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6))**2)/
     .(za(i1,i3)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6)*
     .t(i1,i3,i4)**3)-(zb(i1,i3)**2*
     .(za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))*
     .(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     .(2._dp*za(i1,i3)*zb(i1,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i
     .4))*
     .zb(i3,i4)*zb(i5,i6)*t(i1,i3,i4))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .-(za(i1,i4)*zb(i1,i3)**2*((L1(-t(i1,i3,i4),-s(i3,i4))*za(i3,i4)*
     .(za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))**2)/
     .(s(i3,i4)**2*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     .(L1(-t(i1,i3,i4),-s(i1,i4))*za(i1,i4)*
     .(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6))**2)/
     .(s(i1,i4)**2*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))))/
     .(2._dp*za(i1,i3)*zb(i5,i6)*t(i1,i3,i4))-
     .(Lsm1_2mh(s(i1,i2),t(i1,i3,i4),s(i3,i4),s(i5,i6))*za(i1,i4)**2*
     .(-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*zb(i2,i6)**2*t(i1,i3,
     .i4))/
     .(za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**3*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6))+
     .(L0(-t(i1,i3,i4),-s(i3,i4))*za(i1,i4)**2*zb(i1,i3)*zb(i2,i6)**2*t(
     .i1,i3,i4))/
     .(s(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**2*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6))+
     .(za(i1,i4)*za(i2,i4)*(-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*
     .zb(i2,i6)**2*((L1(-s(i5,i6),-t(i1,i3,i4))*za(i2,i4))/
     .(za(i1,i4)*t(i1,i3,i4)**2)-
     .(2._dp*L0(-s(i5,i6),-t(i1,i3,i4)))/
     .((-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*t(i1,i3,i4)))*
     .t(i1,i3,i4))/
     .(2._dp*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*zb(i5,i6))+
     .(Lnrat(-s(i2,i3),-s(i5,i6))*zb(i1,i3)*zb(i4,i6)*
     .(za(i1,i3)*za(i4,i5)+((-(za(i1,i3)*zb(i1,i6))+
     .za(i3,i4)*zb(i4,i6))*t(i1,i2,i3)-
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*t(i1,i3,i4))/
     .(zb(i1,i4)*zb(i5,i6))))/
     .(2._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))

      tmp=
     .-((Lsm1_2mh(s(i2,i3),t(i1,i2,i4),s(i1,i4),s(i5,i6))*
     .(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))**2)/
     .(za(i5,i6)*zb(i2,i4)*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))**3))

      tmp=tmp
     .-(L0(-t(i1,i2,i4),-s(i1,i4))*za(i2,i4)*zb(i1,i4)*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))**2)/
     .(s(i1,i4)*za(i5,i6)*zb(i2,i4)*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))**2)+
     .(za(i2,i5)*((-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i1,i2)-
     .2._dp*za(i3,i4)*zb(i1,i4)*zb(i2,i3))*zb(i3,i6))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i4)
     .*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .((-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*zb(i1,i2))-
     .(s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i3,i4)*zb(i1,i4))*
     .((za(i2,i5)**2*zb(i2,i3))/za(i5,i6)+
     .(za(i2,i3)*zb(i3,i6)**2)/zb(i5,i6)))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i2,i3)
     .*
     .zb(i2,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))


      tmp=tmp
     .-(Lnrat(-s(i1,i4),-s(i5,i6))*zb(i1,i4)*
     .(((za(i2,i4)*za(i3,i5)+za(i2,i3)*za(i4,i5))*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4)))/
     .(za(i2,i3)*za(i5,i6)*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .(6._dp*za(i1,i4)*((-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i1,i2)-
     .2._dp*za(i3,i4)*zb(i1,i4)*zb(i2,i3))*
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2+
     .(-((za(i1,i4)*za(i2,i5)*zb(i1,i3)*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4)))/
     .za(i5,i6))-
     .(za(i1,i4)*za(i3,i5)*
     .(-(za(i2,i5)*
     .(2._dp*(s(i1,i4)-s(i2,i3)-s(i5,i6))*zb(i1,i3)+
     .za(i2,i4)*zb(i1,i4)*zb(i2,i3)))+
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i3,i5)*
     .zb(i1,i2)*
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))
     .)/(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))/(za(i2,i3)*za(i5,
     .i6))
     .+((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i1,i2)*
     .(-(za(i2,i5)*zb(i2,i6))+za(i3,i5)*zb(i3,i6)))/
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))+
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i5)*za(i3,i4)*
     .(-(za(i1,i5)*zb(i1,i3))+za(i4,i5)*zb(i3,i4)+
     .za(i5,i6)*zb(i3,i6)))/(za(i2,i3)*za(i5,i6))+
     .zb(i3,i6)*(za(i1,i4)*za(i2,i5)*zb(i1,i2)+
     .za(i3,i4)*za(i5,i6)*zb(i3,i6))-
     .4._dp*za(i1,i4)*zb(i2,i6)*
     .(za(i2,i5)*zb(i1,i3)+
     .(za(i2,i3)*za(i5,i6)*zb(i1,i2)*zb(i3,i6))/
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .(zb(i3,i6)*(za(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i2))+za(i3,i4)*zb(i2,i3))*
     .zb(i3,i6)+
     .2._dp*(-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i4)*
     .zb(i2,i6))-
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i4,i5)*
     .zb(i5,i6))))/zb(i5,i6))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*zb(i2,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))-
     .(L1(-s(i5,i6),-t(i1,i2,i4))*za(i3,i5)*zb(i1,i3)*zb(i3,i6))/
     .(zb(i2,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     .t(i1,i2,i4))+(L0(-t(i1,i2,i4),-s(i5,i6))*za(i3,i5)*zb(i1,i3)*
     .(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))*t(i1,i2,i4))/
     .(s(i5,i6)*za(i5,i6)*zb(i2,i4)*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))**2)


      tmp=tmp
     .+(Lnrat(-s(i2,i3),-s(i5,i6))*(-((za(i2,i5)*zb(i2,i3)*
     .(2._dp*((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i3,i5)*
     .zb(i1,i3)+
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*
     .zb(i1,i6))+
     .(za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*
     .(-(za(i1,i5)*zb(i1,i3))+za(i4,i5)*zb(i3,i4))))/
     .za(i5,i6))+(6._dp*
     .(-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i3,i4))+
     .2._dp*za(i1,i4)*za(i2,i3)*zb(i1,i2))*zb(i1,i4)*zb(i2,i3)*
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)-
     .(zb(i3,i6)*(3._dp*(-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*
     .zb(i1,i2)*zb(i3,i6)-
     .za(i2,i3)*zb(i1,i2)*zb(i2,i3)*
     .(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))+
     .za(i3,i4)*zb(i1,i4)*zb(i2,i3)*
     .(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))-
     .2._dp*za(i2,i4)*zb(i1,i4)*
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*zb(i2,i6)-
     .2._dp*za(i3,i5)*zb(i2,i3)*zb(i5,i6))))/zb(i5,i6)+
     .(zb(i1,i2)*(-((za(i2,i5)*za(i3,i5)*zb(i2,i3))/za(i5,i6))-
     .(za(i2,i3)*zb(i2,i6)*zb(i3,i6))/zb(i5,i6))*
     .(-s(i1,i4)**2+2._dp*s(i1,i4)*s(i2,i3)-s(i2,i3)**2+
     .2._dp*s(i1,i4)*s(i5,i6)+2._dp*s(i2,i3)*s(i5,i6)-s(i5,i6)**2+
     .(-s(i1,i4)+s(i2,i3)-s(i5,i6))*
     .(t(i1,i2,i4)-t(i1,i3,i4))))/
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i4)
     .*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))


      tmp=tmp
     .-(I3m(s(i1,i4),s(i2,i3),s(i5,i6))*zb(i1,i4)*
     .((3._dp*za(i1,i4)*(-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*
     .zb(i1,i2))-
     .(s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i3,i4)*zb(i1,i4))*
     .zb(i2,i3)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2-
     .(za(i2,i4)*(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))*
     .zb(i3,i6))/t(i1,i2,i4)-
     .(-(za(i2,i4)*((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i3,i5)*
     .zb(i2,i3)+
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*zb(i2,i6))*
     .zb(i3,i6))+
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     .(3._dp*za(i1,i4)*za(i3,i5)*zb(i1,i2)*zb(i3,i6)+
     .za(i3,i4)*zb(i2,i3)*
     .(-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6)))+
     .za(i1,i4)*za(i2,i5)*zb(i1,i3)*
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*zb(i2,i6)-
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     .zb(i3,i6)-2._dp*za(i3,i5)*zb(i2,i3)*zb(i5,i6))-
     .(za(i1,i4)*zb(i1,i2)*zb(i2,i6)*
     .(-(za(i1,i2)*za(i3,i5)*zb(i1,i3))-
     .za(i2,i4)*za(i3,i5)*zb(i3,i4))*
     .(-t(i1,i2,i4)+t(i1,i3,i4)))/
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     .(zb(i2,i4)*(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))


c --- DSW. Introduce tmp to split up the above expression
      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+((za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))*zb(i2,i4)*
     .(tmp))/
     .(zb(i1,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm+
     .I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     .((za(i2,i4)*za(i3,i5)*zb(i1,i2)*zb(i3,i6))/
     .((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))-
     .(3._dp*za(i5,i6)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .zb(i2,i6)*(-((za(i1,i4)*zb(i1,i2)*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i2,i6))+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i1,i5)*zb(i5,i6)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     .za(i2,i4)*zb(i1,i2)*
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6))-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))-
     .(3._dp*s(i1,i2)*za(i3,i5)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3
     .))*
     .zb(i5,i6)*(za(i5,i6)*zb(i1,i2)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))+
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*zb(i1,i2))-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))*
     .t(i1,i2,i3)))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))-
     .(za(i3,i5)*zb(i1,i2)*((za(i1,i2)*za(i3,i4)*zb(i1,i3)+
     .za(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))+
     .za(i2,i4)*(-(za(i1,i4)*zb(i1,i4))-za(i2,i4)*zb(i2,i4)))*
     .zb(i3,i6)*t(i1,i2,i4)-
     .((-(za(i1,i2)*za(i3,i4)*zb(i1,i4))+
     .za(i2,i3)*(-(za(i1,i3)*zb(i1,i3))-
     .za(i2,i3)*zb(i2,i3))+
     .za(i2,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))
     .)*(-(za(i1,i4)*zb(i1,i3)*zb(i4,i6))-za(i2,i4)*zb(i2,i3)*zb(i4,i6))
     .*
     .t(i1,i2,i4))/(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))+
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6))-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6)-
     .(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))*t(i1,i2,i3)+
     .(-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6))*t(i1,i2,i4))-
     .za(i1,i2)*zb(i5,i6)*
     .(3._dp*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     .za(i4,i5)*zb(i1,i3)*(-t(i1,i2,i3)+t(i1,i2,i4)))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))+
     .(zb(i2,i6)*(-(za(i2,i4)*(za(i1,i4)*za(i5,i6)*zb(i1,i3)*zb(i1,i6)-
     .(za(i2,i5)*zb(i1,i2)+za(i5,i6)*zb(i1,i6))*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     .za(i2,i4)*za(i5,i6)*zb(i1,i2)*zb(i3,i6)))+
     .(za(i1,i4)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .(za(i1,i5)*zb(i1,i2)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))+
     .3._dp*za(i1,i4)*za(i5,i6)*zb(i1,i3)*zb(i2,i6)-
     .(za(i1,i4)*za(i5,i6)*zb(i2,i3)*zb(i2,i6)*
     .(t(i1,i3,i4)-t(i2,i3,i4)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .-(Lnrat(-s(i3,i4),-s(i5,i6))*zb(i2,i3)*
     .(-2._dp*(-(za(i1,i4)*zb(i1,i3))+za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))+
     .4._dp*(-(za(i1,i2)*za(i4,i5)*zb(i1,i3)*zb(i2,i6))+
     .za(i1,i4)*za(i2,i5)*zb(i1,i2)*zb(i3,i6))+
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)+
     .(4._dp*s(i5,i6)*za(i1,i4)*zb(i2,i3))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))*
     .(-((za(i1,i5)*za(i2,i5)*zb(i1,i2))/za(i5,i6))+
     .(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6))-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     .((za(i1,i5)*za(i4,i5)*zb(i1,i3))/za(i5,i6)+
     .(2._dp*za(i1,i5)*zb(i2,i3)*
     .(za(i1,i5)*za(i2,i4)*zb(i1,i2)+
     .za(i1,i5)*za(i3,i4)*zb(i1,i3)-
     .za(i2,i4)*za(i5,i6)*zb(i2,i6)))/
     .(za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))-
     .(za(i2,i4)*zb(i2,i6)*zb(i3,i6))/zb(i5,i6))+
     .(8._dp*s(i1,i2)*za(i1,i5)*za(i4,i5)*zb(i2,i3)*zb(i5,i6))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .(3._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(t(i1,i3,i4)-t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i4)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm-
     .(Ls1(-s(i2,i3),-t(i2,i3,i4),-s(i3,i4),-t(i2,i3,i4))*za(i2,i4)**2*z
     .b(i2,i3)**2*
     .(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))**2)/
     .(za(i5,i6)*zb(i2,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     .t(i2,i3,i4)**3)-(za(i2,i4)**2*
     .(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))*
     .(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     .(2._dp*za(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i4)*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*t(i2,i3,i4))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+(za(i2,i4)**2*zb(i2,i3)*(-((L1(-t(i2,i3,i4),-s(i3,i4))*
     .(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2*zb(i3,i4))/
     .(s(i3,i4)**2*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))+
     .(L1(-t(i2,i3,i4),-s(i2,i3))*zb(i2,i3)*
     .(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))**2)/
     .(s(i2,i3)**2*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))))/
     .(2._dp*za(i5,i6)*zb(i2,i4)*t(i2,i3,i4))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm-
     .(L0(-t(i2,i3,i4),-s(i3,i4))*za(i1,i5)**2*za(i2,i4)*zb(i2,i3)**2*t(
     .i2,i3,i4))/
     .(s(i3,i4)*za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**
     .2*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm-
     .(Lsm1_2mh(s(i1,i2),t(i2,i3,i4),s(i3,i4),s(i5,i6))*za(i1,i5)**2*zb(
     .i2,i3)**2*
     .(za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4))*t(i2,i3,i4))/
     .(za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))**3*zb(i3,i
     .4)*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+(za(i1,i5)**2*zb(i1,i3)*zb(i2,i3)*
     .(za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4))*
     .((L1(-s(i5,i6),-t(i2,i3,i4))*zb(i1,i3))/(zb(i2,i3)*t(i2,i3,i4)**2)
     .-
     .(2._dp*L0(-s(i5,i6),-t(i2,i3,i4)))/
     .((-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*t(i2,i3,i4)))*
     .t(i2,i3,i4))/
     .(2._dp*za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i3,
     .i4)*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .-(Lnrat(-s(i1,i2),-s(i5,i6))*zb(i2,i3)*
     .(((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i2)*
     .(-(za(i2,i5)*(-(za(i1,i5)*za(i2,i4))+za(i1,i2)*za(i4,i5))*
     .zb(i2,i3))-za(i1,i5)*za(i2,i5)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))+
     .(za(i1,i4)*za(i1,i5)**2*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i2,i3))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)
     .)
     .-(2._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i1,i2)*za(i4,i5)*
     .(zb(i1,i6)*zb(i2,i3)+zb(i1,i3)*zb(i2,i6)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     .(za(i1,i4)*zb(i2,i3)*(-((za(i1,i5)*za(i2,i5)*zb(i1,i2))/
     .za(i5,i6))-(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))-
     .(za(i2,i4)*zb(i2,i6)*zb(i3,i6))/zb(i5,i6)+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*
     .(zb(i1,i3)*zb(i2,i6)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))-
     .za(i1,i4)*zb(i2,i3)*
     .(zb(i1,i6)**2+((-(za(i2,i3)*zb(i1,i3))-
     .za(i2,i4)*zb(i1,i4))*zb(i2,i6)**2)/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     .zb(i1,i2)*zb(i3,i6)*
     .(-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i5,i6)
     .)
     .+(6._dp*s(i1,i2)*s(i3,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3
     .))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(t(i1,i3,i4)-t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2-
     .(s(i1,i2)*(-((za(i4,i5)**2*zb(i3,i4))/za(i5,i6))-
     .(za(i3,i4)*zb(i3,i6)**2)/zb(i5,i6))*t(i2,i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*s(i3,i4)*zb(i2,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4
     .)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm-
     .(Lnrat(-s(i1,i2),-s(i5,i6))*za(i1,i4)*
     .((za(i1,i5)*za(i4,i5)*zb(i1,i3))/za(i5,i6)+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i2)*
     .(za(i1,i5)*za(i2,i4)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))-
     .za(i1,i4)*zb(i2,i3)*
     .(za(i2,i5)**2+(za(i1,i5)**2*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))-
     .za(i1,i2)*za(i4,i5)*
     .(-(za(i1,i5)*zb(i1,i3))+za(i4,i5)*zb(i3,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)
     .)
     .+(2._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*
     .(za(i1,i5)*za(i2,i4)+za(i1,i4)*za(i2,i5))*zb(i1,i2)*zb(i3,i6)
     .)/(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5
     .,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     .(za(i1,i4)*zb(i2,i3)*(-((za(i1,i5)*za(i2,i5)*zb(i1,i2))/
     .za(i5,i6))-(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*
     .(-(zb(i1,i6)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .zb(i2,i6))+(za(i1,i4)*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .zb(i2,i3)*zb(i2,i6)**2)/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))-
     .za(i1,i4)*zb(i1,i6)*
     .(-(zb(i1,i3)*zb(i2,i6))-zb(i1,i2)*zb(i3,i6))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i5,i6)
     .)
     .-(s(i1,i2)*(-((za(i4,i5)**2*zb(i3,i4))/za(i5,i6))-
     .(za(i3,i4)*zb(i3,i6)**2)/zb(i5,i6))*t(i1,i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     .(6._dp*s(i1,i2)*s(i3,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)
     .)*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(-t(i1,i3,i4)+t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2))/
     .(2._dp*s(i3,i4)*za(i1,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4
     .)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .-(Lnrat(-s(i3,i4),-s(i5,i6))*za(i1,i4)*
     .(-2._dp*(za(i1,i4)*zb(i1,i3)-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))-
     .(8._dp*s(i1,i2)*za(i1,i4)*za(i5,i6)*zb(i2,i6)*zb(i3,i6))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .4._dp*(-(za(i1,i2)*za(i4,i5)*zb(i1,i6)*zb(i2,i3))+
     .za(i1,i5)*za(i2,i4)*zb(i1,i2)*zb(i3,i6))+
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)+
     .(4._dp*s(i5,i6)*za(i1,i4)*zb(i2,i3))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))*
     .((za(i1,i5)*za(i2,i5)*zb(i1,i2))/za(i5,i6)-
     .(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6))-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     .((za(i1,i5)*za(i4,i5)*zb(i1,i3))/za(i5,i6)-
     .(za(i2,i4)*zb(i2,i6)*zb(i3,i6))/zb(i5,i6)-
     .(2._dp*za(i1,i4)*zb(i2,i6)*
     .(-(za(i1,i2)*zb(i1,i3)*zb(i2,i6))-
     .za(i2,i4)*zb(i2,i6)*zb(i3,i4)+
     .za(i1,i5)*zb(i1,i3)*zb(i5,i6)))/
     .((-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i5,i6)))+
     .(3._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(-t(i1,i3,i4)+t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i3)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+(Lnrat(-s(i1,i4),-s(i5,i6))*za(i2,i4)*za(i3,i5)*
     .(zb(i2,i4)*zb(i3,i6)-((-(za(i2,i5)*zb(i2,i4))-
     .za(i3,i5)*zb(i3,i4))*t(i1,i2,i4)-
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))*t(i2,i3,i4))/
     .(za(i2,i3)*za(i5,i6))))/
     .(2._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      tmp=
     .(za(i4,i5)*zb(i1,i6)*(-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*
     .za(i1,i2))+2._dp*za(i1,i4)*za(i2,i3)*zb(i3,i4)))/
     .((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i3)
     .*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      tmp=tmp
     .+(((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i2)*zb(i1,i4)+
     .(-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*zb(i3,i4))*
     .(-((za(i4,i5)**2*zb(i1,i4))/za(i5,i6))-
     .(za(i1,i4)*zb(i1,i6)**2)/zb(i5,i6)))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i3)
     .*
     .zb(i1,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      tmp=tmp
     .+(L0(-t(i1,i2,i3),-s(i2,i3))*za(i2,i3)*zb(i1,i3)*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     .(s(i2,i3)*za(i1,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**2*
     .zb(i5,i6))


      tmp=tmp
     .+(Lsm1_2mh(s(i1,i4),t(i1,i2,i3),s(i2,i3),s(i5,i6))*
     .(-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     .(za(i1,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**3*
     .zb(i5,i6))


      tmp=tmp
     .-(Lnrat(-s(i2,i3),-s(i5,i6))*za(i2,i3)*
     .((6._dp*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i2,i3)*
     .(-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i2))+
     .2._dp*za(i1,i4)*za(i2,i3)*zb(i3,i4))*
     .(-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2-
     .((za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))*
     .(zb(i1,i4)*zb(i3,i6)+zb(i1,i3)*zb(i4,i6)))/
     .(zb(i1,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     .zb(i5,i6))+(-((za(i4,i5)*
     .(za(i4,i5)*zb(i1,i4)*
     .(za(i1,i2)*zb(i2,i3)-za(i1,i4)*zb(i3,i4))+
     .2._dp*(-((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i5)*
     .zb(i1,i3))+
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*
     .zb(i3,i6))))/za(i5,i6))-
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i2)*zb(i2,i3)*
     .(-(za(i1,i5)*zb(i1,i6))+za(i4,i5)*zb(i4,i6)))/
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))+
     .(za(i2,i4)*zb(i1,i6)*zb(i2,i3)*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6)))/
     .zb(i5,i6)+(zb(i2,i3)*zb(i4,i6)*
     .(-((2._dp*(-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i4)+
     .za(i1,i4)*za(i2,i3)*zb(i1,i3))*zb(i1,i6))-
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i2)*
     .(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     .zb(i4,i6))/
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     .(zb(i1,i4)*zb(i5,i6))+
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i1,i6)*zb(i3,i4)*
     .(-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6)-
     .za(i4,i5)*zb(i5,i6)))/(zb(i1,i4)*zb(i5,i6))+
     .za(i4,i5)*(-(za(i1,i2)*zb(i1,i6)*zb(i2,i3))+
     .za(i4,i5)*zb(i3,i4)*zb(i5,i6))-
     .4._dp*za(i1,i5)*zb(i2,i3)*
     .(za(i2,i4)*zb(i1,i6)+
     .(za(i1,i2)*za(i4,i5)*zb(i1,i4)*zb(i5,i6))/
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*za(i1,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      tmp=tmp
     .-(L1(-s(i5,i6),-t(i1,i2,i3))*za(i2,i4)*za(i4,i5)*zb(i4,i6))/
     .(za(i1,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*t(i1,i2,i3))


      tmp=tmp
     .-(L0(-t(i1,i2,i3),-s(i5,i6))*za(i2,i4)*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))*zb(i4,i6)*
     .t(i1,i2,i3))/
     .(s(i5,i6)*za(i1,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**2*
     .zb(i5,i6))


      tmp=tmp
     .+(Lnrat(-s(i1,i4),-s(i5,i6))*
     .((6._dp*za(i1,i4)*za(i2,i3)*
     .(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     .(-2._dp*za(i1,i2)*zb(i1,i4)*zb(i2,i3)+
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i3,i4))*
     .(-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)+
     .(za(i4,i5)*(-3._dp*(s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i2)*
     .za(i4,i5)*zb(i1,i4)+
     .za(i1,i2)*za(i1,i4)*
     .(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     .zb(i1,i4)-
     .za(i1,i4)*za(i2,i3)*
     .(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     .zb(i3,i4)-
     .2._dp*za(i2,i3)*zb(i1,i3)*
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i5)+
     .2._dp*za(i1,i4)*za(i5,i6)*zb(i4,i6))))/za(i5,i6)+
     .(za(i1,i4)*zb(i1,i6)*
     .((-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))*
     .(-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6))+
     .2._dp*((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i4)*
     .zb(i4,i6)-
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i5)*zb(i5,i6))
     .))/zb(i5,i6)-(za(i1,i2)*((za(i1,i5)*za(i4,i5)*zb(i1,i4))/za(i5,i6)
     .+
     .(za(i1,i4)*zb(i1,i6)*zb(i4,i6))/zb(i5,i6))*
     .(-s(i1,i4)**2+2._dp*s(i1,i4)*s(i2,i3)-s(i2,i3)**2+
     .2._dp*s(i1,i4)*s(i5,i6)+2._dp*s(i2,i3)*s(i5,i6)-s(i5,i6)**2+
     .(s(i1,i4)-s(i2,i3)-s(i5,i6))*(t(i1,i2,i3)-t(i2,i3,i4))
     .))/(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     .(2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i3)
     .*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


      tmp=tmp
     .-(I3m(s(i2,i3),s(i1,i4),s(i5,i6))*za(i2,i3)*
     .((3._dp*za(i1,i4)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     .zb(i2,i3)*((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i2)*
     .zb(i1,i4)+
     .(-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*zb(i3,i4))*
     .(-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2-
     .(za(i4,i5)*zb(i1,i3)*
     .(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6)))/t(i1,i2,i3)-
     .(za(i2,i4)*zb(i1,i6)*zb(i2,i3)*
     .((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i5)-
     .za(i4,i5)*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))+
     .2._dp*za(i1,i4)*za(i5,i6)*zb(i4,i6))+
     .(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i4)*zb(i3,i4)*
     .(-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))-
     .3._dp*za(i1,i2)*za(i4,i5)*zb(i2,i3)*zb(i4,i6))-
     .za(i4,i5)*zb(i1,i3)*
     .((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i4,i6)-
     .(-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*zb(i5,i6))+
     .(za(i1,i2)*za(i1,i5)*zb(i2,i3)*
     .(za(i2,i4)*zb(i1,i2)*zb(i4,i6)+
     .za(i3,i4)*zb(i1,i3)*zb(i4,i6))*
     .(-t(i1,i2,i3)+t(i2,i3,i4)))/
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))/
     .(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     .2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     .(za(i1,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))


c --- DSW. Introduce tmp to split up the above expression
      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+(za(i1,i3)*(-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))*
     .(tmp))/
     .(za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))


      zaj_virt_a6_slc_fsc_pm=zaj_virt_a6_slc_fsc_pm
     .+I3m(s(i1,i2),s(i3,i4),s(i5,i6))
     .*(-((za(i1,i2)*za(i4,i5)*zb(i1,i3)*zb(i4,i6))/
     .((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))
     .+
     .(3._dp*za(i1,i5)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i2)*zb(i1,i3)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*zb(i1,i2)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6)))+
     .(za(i1,i2)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .zb(i2,i3)*((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i5)*
     .zb(i1,i2)-(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*
     .zb(i2,i6)))/(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))
     .)*zb(i5,i6))/((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
     .+
     .(3._dp*s(i1,i2)*za(i5,i6)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3
     .))*
     .zb(i4,i6)*(za(i1,i2)*(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*
     .zb(i1,i2))-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))*zb(i5,i6)
     .+((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))*t(i1,i2,i4)
     .))/((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(
     .i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
     .+
     .(za(i1,i2)*zb(i4,i6)*(-(((-(za(i1,i4)*za(i3,i5)*zb(i1,i3))-
     .za(i2,i4)*za(i3,i5)*zb(i2,i3))*
     .(zb(i1,i3)*(-(za(i1,i3)*zb(i1,i4))-
     .za(i2,i3)*zb(i2,i4))+
     .zb(i1,i4)*(-(za(i1,i4)*zb(i1,i4))-
     .za(i2,i4)*zb(i2,i4))-za(i2,i3)*zb(i1,i2)*zb(i3,i4)
     .)*t(i1,i2,i3))/(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     .za(i4,i5)*(zb(i1,i3)*
     .(-(za(i1,i3)*zb(i1,i3))-za(i2,i3)*zb(i2,i3))+
     .zb(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))+
     .za(i2,i4)*zb(i1,i2)*zb(i3,i4))*t(i1,i2,i3)-
     .za(i5,i6)*zb(i1,i2)*
     .(3._dp*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6))-
     .za(i2,i4)*zb(i3,i6)*(t(i1,i2,i3)-t(i1,i2,i4)))+
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*zb(i1,i2)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6)+
     .(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*t(i1,i2,i3)-
     .(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*t(i1,i2,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
     .+
     .(za(i1,i5)*(-(zb(i1,i3)*(-(za(i1,i2)*za(i4,i5)*zb(i1,i3)*zb(i5,i6)
     .)-
     .za(i2,i4)*za(i2,i5)*zb(i2,i3)*zb(i5,i6)-
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i2)*zb(i1,i6))-za(i2,i5)*zb(i5,i6))))+
     .((-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i2,i3)*
     .(-(za(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .zb(i2,i6))-3._dp*za(i1,i5)*za(i2,i4)*zb(i2,i3)*zb(i5,i6)+
     .(za(i1,i4)*za(i1,i5)*zb(i2,i3)*zb(i5,i6)*
     .(-t(i1,i3,i4)+t(i2,i3,i4)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))

      end function

      function zaj_virt_a6_slc_fsc_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_slc_fsc_mp

        complex(dp) :: L0,L1,Lsm1_2mh,I3m,Lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      zaj_virt_a6_slc_fsc_mp= 
     .-(za(i3,i5)*(-(za(i2,i5)*za(i3,i4))+za(i2,i3)*za(i4,i5)))/
     .(2._dp*za(i1,i4)*za(i3,i4)*za(i5,i6)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))+
     .(Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*
     .(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))**2)/
     .(za(i5,i6)*zb(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .*3)-
     .(L0(-t(i1,i2,i3),-s(i1,i2))*za(i2,i3)*zb(i1,i2)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))**2)/
     .(s(i1,i2)*za(i5,i6)*zb(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2)+
     .(Lnrat(-s(i1,i2),-s(i5,i6))*za(i4,i5)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     .(2._dp*s(i3,i4)*za(i5,i6)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3
     .))*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))-
     .(za(i3,i5)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i4)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i3))*zb(i4,i6))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i4)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     .(za(i3,i5)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i3)+
     .2._dp*za(i2,i4)*zb(i1,i2)*zb(i3,i4))*zb(i4,i6))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(za(i1,i3)*(-((za(i3,i5)*(-((za(i2,i5)*
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i3)*
     .zb(i1,i2))+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*zb(i2,i4)
     .))/za(i5,i6))+zb(i1,i6)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i1,i3)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i2,i4))))/za(i3,i4))+
     .(zb(i4,i6)*(-(za(i2,i5)*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i2,i4)+
     .2._dp*za(i1,i3)*zb(i1,i2)*zb(i3,i4)))-
     .(zb(i1,i6)*((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*
     .zb(i2,i4)-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i1,i3)*zb(i3,i4)))/
     .zb(i5,i6)))/zb(i3,i4)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i4)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     .(zb(i2,i4)*(-((za(i3,i5)*((za(i2,i5)*
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i3)*
     .zb(i1,i2))+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*zb(i2,i4)))/
     .za(i5,i6)-zb(i1,i6)*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i1,i3)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i2,i4))))/za(i3,i4))+
     .(zb(i4,i6)*(za(i2,i5)*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i2,i4)+
     .2._dp*za(i1,i3)*zb(i1,i2)*zb(i3,i4))+
     .(zb(i1,i6)*((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*
     .zb(i2,i4)-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i1,i3)*zb(i3,i4)))/
     .zb(i5,i6)))/zb(i3,i4)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i3)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp-
     .(((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i1,i2)-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*zb(i1,i3))*
     .((za(i3,i5)**2*zb(i3,i4))/za(i5,i6)+
     .(za(i3,i4)*zb(i4,i6)**2)/zb(i5,i6)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i3,i4)*zb(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))+
     .((-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i3))+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i3,i4))*
     .((za(i3,i5)**2*zb(i3,i4))/za(i5,i6)+
     .(za(i3,i4)*zb(i4,i6)**2)/zb(i5,i6)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i4)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i3,i4))-
     .(L0(-t(i1,i2,i4),-s(i1,i2))*za(i1,i2)*zb(i1,i4)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))**2)/
     .(s(i1,i2)*za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**
     .2*
     .zb(i5,i6))-(Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),s(i5,i6))*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))**2*
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))/
     .(za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**3*zb(i5,i
     .6))-
     .((zb(i1,i6)*zb(i3,i4)+zb(i1,i4)*zb(i3,i6))*zb(i4,i6))/
     .(2._dp*zb(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i3,
     .i4)*
     .zb(i5,i6))-(Lnrat(-s(i1,i2),-s(i5,i6))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*zb(i3,i6)*
     .(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     .(2._dp*s(i3,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i5,i6))

      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp-
     .(Lnrat(-s(i1,i2),-s(i5,i6))*za(i1,i2)*
     .((-6._dp*zb(i1,i2)*((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i4)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2+
     .((-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*
     .(-(zb(i1,i6)*zb(i3,i4))-zb(i1,i4)*zb(i3,i6)))/
     .((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i3,i4)*
     .zb(i5,i6))+(-((za(i3,i5)*
     .(2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i4,i5)*
     .zb(i1,i4)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))
     .-za(i3,i5)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i3,i4)))/
     .za(i5,i6))-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)*
     .zb(i1,i2)*(za(i3,i5)*zb(i3,i6)-za(i4,i5)*zb(i4,i6)))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     .(za(i2,i3)*zb(i1,i2)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*zb(i4,i6))/
     .zb(i5,i6)+(zb(i1,i2)*zb(i3,i6)*
     .(((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .zb(i3,i6))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     .(2._dp*(s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i3)-
     .za(i1,i2)*za(i3,i4)*zb(i1,i4))*zb(i4,i6)))/
     .(zb(i3,i4)*zb(i5,i6))-
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i3)*zb(i4,i6)*
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)-
     .za(i3,i5)*zb(i5,i6)))/(zb(i3,i4)*zb(i5,i6))+
     .za(i3,i5)*(-(za(i2,i4)*zb(i1,i2)*zb(i4,i6))+
     .za(i3,i5)*zb(i1,i3)*zb(i5,i6))+
     .4._dp*za(i4,i5)*zb(i1,i2)*
     .(za(i2,i3)*zb(i4,i6)+
     .(za(i2,i4)*za(i3,i5)*zb(i3,i4)*zb(i5,i6))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(Lnrat(-s(i1,i2),-s(i5,i6))*zb(i1,i2)*
     .(((za(i2,i5)*za(i3,i4)-za(i2,i3)*za(i4,i5))*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)))/
     .(za(i3,i4)*za(i5,i6)*(-(za(i1,i4)*zb(i1,i3))-
     .za(i2,i4)*zb(i2,i3)))+
     .(6._dp*za(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i3)+
     .2._dp*za(i2,i4)*zb(i1,i2)*zb(i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2+
     .(-((za(i1,i2)*za(i3,i5)*zb(i1,i4)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)))/za(i5,i6))
     .-(za(i1,i2)*za(i4,i5)*(((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i4,i5)*
     .zb(i1,i3)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     .za(i3,i5)*(2._dp*(s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i4)-
     .za(i2,i3)*zb(i1,i2)*zb(i3,i4))))/(za(i3,i4)*za(i5,i6))
     .+((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i3)*
     .(-(za(i3,i5)*zb(i3,i6))+za(i4,i5)*zb(i4,i6)))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)*za(i3,i5)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)+
     .za(i5,i6)*zb(i4,i6)))/(za(i3,i4)*za(i5,i6))+
     .zb(i4,i6)*(za(i1,i2)*za(i3,i5)*zb(i1,i3)-
     .za(i2,i4)*za(i5,i6)*zb(i4,i6))-
     .4._dp*za(i1,i2)*zb(i3,i6)*
     .(za(i3,i5)*zb(i1,i4)+
     .(za(i3,i4)*za(i5,i6)*zb(i1,i3)*zb(i4,i6))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))+
     .(zb(i4,i6)*(za(i3,i4)*
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i4,i6)
     .+2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i3)*zb(i3,i6)-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .zb(i5,i6))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*zb(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))+
     .(L1(-s(i5,i6),-t(i1,i2,i3))*za(i4,i5)*zb(i1,i4)*zb(i4,i6))/
     .(zb(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*t(i1,i2,i3
     .))

      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp-
     .(L0(-t(i1,i2,i3),-s(i5,i6))*za(i4,i5)*zb(i1,i4)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))*t(i1,i2,i3))/
     .(s(i5,i6)*za(i5,i6)*zb(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2)-
     .(Lnrat(-s(i3,i4),-s(i5,i6))*(-((za(i3,i5)*
     .(2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i4,i5)*zb(i1,i4)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))+
     .(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))*zb(i3,i4))/
     .za(i5,i6))+(6._dp*zb(i1,i2)*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i4)+
     .2._dp*za(i1,i2)*za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*zb(i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     .(zb(i4,i6)*(-(za(i2,i4)*zb(i1,i2)*
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*zb(i3,i4))
     .-za(i3,i4)*zb(i1,i3)*(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*
     .zb(i3,i4)+3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*
     .zb(i1,i3)*zb(i4,i6)+
     .2._dp*za(i2,i3)*zb(i1,i2)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i3,i6)-
     .2._dp*za(i4,i5)*zb(i3,i4)*zb(i5,i6))))/zb(i5,i6)+
     .(zb(i1,i3)*(-((za(i3,i5)*za(i4,i5)*zb(i3,i4))/za(i5,i6))-
     .(za(i3,i4)*zb(i3,i6)*zb(i4,i6))/zb(i5,i6))*
     .(-s(i1,i2)**2+2._dp*s(i1,i2)*s(i3,i4)-s(i3,i4)**2+
     .2._dp*s(i1,i2)*s(i5,i6)+2._dp*s(i3,i4)*s(i5,i6)-s(i5,i6)**2+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*(t(i1,i2,i3)-t(i1,i2,i4))))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp-
     .(I3m(s(i1,i2),s(i3,i4),s(i5,i6))*za(i1,i2)*
     .((3._dp*za(i3,i4)*zb(i1,i2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i3))+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2-
     .((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(za(i3,i4)*zb(i1,i3)*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))-
     .3._dp*za(i2,i4)*za(i3,i5)*zb(i1,i2)*zb(i3,i6))-
     .za(i2,i3)*zb(i1,i2)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i4,i5)-
     .za(i3,i5)*(-(za(i1,i4)*zb(i1,i3))-
     .za(i2,i4)*zb(i2,i3))-2._dp*za(i3,i4)*za(i5,i6)*zb(i3,i6))*
     .zb(i4,i6)+za(i3,i5)*zb(i1,i4)*
     .(-((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*zb(i3,i6))-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i4,i5)*zb(i5,i6))+
     .(za(i2,i4)*za(i4,i5)*zb(i1,i2)*
     .(-(za(i1,i3)*zb(i1,i4)*zb(i3,i6))-
     .za(i2,i3)*zb(i2,i4)*zb(i3,i6))*(t(i1,i2,i3)-t(i1,i2,i4)))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     .(za(i3,i5)*zb(i1,i4)*(-(za(i1,i4)*zb(i1,i6))-
     .za(i2,i4)*zb(i2,i6)))/t(i1,i2,i4)))/
     .(za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))

      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(L1(-s(i5,i6),-t(i1,i2,i4))*za(i2,i3)*za(i3,i5)*zb(i3,i6))/
     .(za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*t(i1,i2,i4
     .))+
     .(L0(-t(i1,i2,i4),-s(i5,i6))*za(i2,i3)*
     .(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*zb(i3,i6)*t(i1,i2,i4)
     .)/
     .(s(i5,i6)*za(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**
     .2*
     .zb(i5,i6))-(Lnrat(-s(i3,i4),-s(i5,i6))*
     .((6._dp*za(i1,i2)*za(i3,i4)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i3)+
     .2._dp*za(i2,i4)*zb(i1,i2)*zb(i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     .(za(i3,i5)*(-(za(i1,i2)*za(i3,i4)*zb(i1,i3)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))-
     .3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i4)*za(i3,i5)*
     .zb(i3,i4)-za(i2,i4)*za(i3,i4)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))*zb(i3,i4)-
     .2._dp*za(i1,i2)*zb(i1,i4)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i4,i5)-
     .2._dp*za(i3,i4)*za(i5,i6)*zb(i3,i6))))/za(i5,i6)-
     .(za(i3,i4)*zb(i4,i6)*((-(za(i1,i3)*zb(i1,i6))-
     .za(i2,i3)*zb(i2,i6))*
     .(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))+
     .2._dp*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i3)*zb(i3,i6)-
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6))))/
     .zb(i5,i6)+(za(i2,i4)*
     .(-((za(i3,i5)*za(i4,i5)*zb(i3,i4))/za(i5,i6))-
     .(za(i3,i4)*zb(i3,i6)*zb(i4,i6))/zb(i5,i6))*
     .(-s(i1,i2)**2+2._dp*s(i1,i2)*s(i3,i4)-s(i3,i4)**2+
     .2._dp*s(i1,i2)*s(i5,i6)+2._dp*s(i3,i4)*s(i5,i6)-s(i5,i6)**2+
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*(-t(i1,i2,i3)+t(i1,i2,i4))))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i4)*
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))

      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(I3m(s(i1,i2),s(i3,i4),s(i5,i6))*zb(i1,i2)*
     .((3._dp*za(i1,i2)*((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i1,i2)-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i4)*zb(i1,i3))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*zb(i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2+
     .(za(i2,i3)*(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))*
     .zb(i4,i6))/t(i1,i2,i3)-
     .(za(i2,i3)*((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i4,i5)*zb(i3,i4)+
     .(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i3,i6))*
     .zb(i4,i6)+(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i2,i4)*(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .zb(i3,i4))+3._dp*za(i1,i2)*za(i4,i5)*zb(i1,i3)*zb(i4,i6))+
     .za(i1,i2)*za(i3,i5)*zb(i1,i4)*
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i3,i6)-
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i4,i6)-
     .2._dp*za(i4,i5)*zb(i3,i4)*zb(i5,i6))-
     .(za(i1,i2)*zb(i1,i3)*
     .(-(za(i1,i3)*za(i4,i5)*zb(i1,i4))-
     .za(i2,i3)*za(i4,i5)*zb(i2,i4))*zb(i3,i6)*
     .(-t(i1,i2,i3)+t(i1,i2,i4)))/
     .(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(zb(i2,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))+
     .(za(i1,i3)*zb(i1,i4)*(za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))*
     .((L0(-t(i1,i3,i4),-s(i3,i4))*
     .((-2._dp*za(i1,i3)*zb(i2,i6))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6))/t(i1,i3,i4)))/
     .s(i3,i4)-(L1(-s(i3,i4),-t(i1,i3,i4))*za(i1,i3)*zb(i1,i6))/t(i1,i3,
     .i4)**2
     .))/(2._dp*za(i1,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(
     .i5,i6))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(za(i1,i3)*za(i2,i3)*zb(i2,i6)**2*
     .((-2._dp*L0(-s(i5,i6),-t(i1,i3,i4))*za(i1,i3))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .(L1(-s(i5,i6),-t(i1,i3,i4))*za(i2,i3))/t(i1,i3,i4)))/
     .(2._dp*za(i1,i4)*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i
     .4))*
     .zb(i5,i6))-(za(i1,i3)*zb(i4,i6)*
     .(-(zb(i4,i6)/zb(i3,i4))+
     .(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6))/t(i1,i3,i4)))/
     .(2._dp*za(i1,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i5,
     .i6))-
     .(Lsm1_2mh(s(i1,i2),t(i1,i3,i4),s(i3,i4),s(i5,i6))*za(i1,i3)**3*zb(
     .i2,i6)**2*
     .t(i1,i3,i4))/
     .(za(i1,i4)*za(i3,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .*3*
     .zb(i5,i6))+(I3m(s(i1,i2),s(i3,i4),s(i5,i6))*zb(i2,i4)*
     .((s(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6)))/2-
     .za(i1,i2)*zb(i1,i4)*((-(za(i1,i3)*zb(i1,i2))-
     .za(i3,i4)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))-
     .za(i1,i5)*za(i3,i4)*zb(i1,i2)*zb(i4,i6))-
     .(za(i1,i3)*za(i1,i5)*(-(za(i2,i3)*zb(i1,i3))-
     .za(i2,i4)*zb(i1,i4))*zb(i2,i4)*
     .(-(zb(i2,i3)*(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)))-
     .zb(i2,i4)*(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))-
     .zb(i1,i2)*(za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6))))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .za(i1,i5)*za(i2,i3)*(zb(i2,i4)*
     .(-(zb(i1,i3)*(-(za(i1,i3)*zb(i1,i6))-
     .za(i2,i3)*zb(i2,i6)))-
     .zb(i1,i4)*(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))+
     .zb(i1,i2)*(za(i2,i3)*zb(i3,i6)+za(i2,i4)*zb(i4,i6)))-
     .zb(i1,i2)*(-(zb(i1,i4)*
     .(za(i1,i3)*zb(i3,i6)+za(i1,i4)*zb(i4,i6)))-
     .zb(i2,i4)*(za(i2,i3)*zb(i3,i6)+za(i2,i4)*zb(i4,i6))))+
     .(3._dp*s(i1,i2)*(s(i1,i2)-s(i3,i4)-s(i5,i6))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(t(i1,i3,i4)-t(i2,i3,i4)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i3)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(Lnrat(-s(i3,i4),-s(i5,i6))*zb(i2,i4)*
     .(-2._dp*(-(za(i1,i3)*zb(i1,i4))+za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))+
     .4._dp*(-(za(i1,i2)*za(i3,i5)*zb(i1,i4)*zb(i2,i6))+
     .za(i1,i3)*za(i2,i5)*zb(i1,i2)*zb(i4,i6))+
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)+
     .(4._dp*s(i5,i6)*za(i1,i3)*zb(i2,i4))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))*
     .(-((za(i1,i5)*za(i2,i5)*zb(i1,i2))/za(i5,i6))+
     .(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6))-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     .((za(i1,i5)*za(i3,i5)*zb(i1,i4))/za(i5,i6)+
     .(2._dp*za(i1,i5)*zb(i2,i4)*
     .(za(i1,i5)*za(i2,i3)*zb(i1,i2)-
     .za(i1,i5)*za(i3,i4)*zb(i1,i4)-
     .za(i2,i3)*za(i5,i6)*zb(i2,i6)))/
     .(za(i5,i6)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))-
     .(za(i2,i3)*zb(i2,i6)*zb(i4,i6))/zb(i5,i6))+
     .(8._dp*s(i1,i2)*za(i1,i5)*za(i3,i5)*zb(i2,i4)*zb(i5,i6))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .(3._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(t(i1,i3,i4)-t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i2,i3)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))-
     .(za(i2,i3)*zb(i2,i4)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))*
     .((L0(-t(i2,i3,i4),-s(i3,i4))*
     .((-2._dp*za(i1,i5)*zb(i2,i4))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))/t(i2,i3,i4)))/
     .s(i3,i4)-(L1(-s(i3,i4),-t(i2,i3,i4))*za(i2,i5)*zb(i2,i4))/t(i2,i3,
     .i4)**2
     .))/(2._dp*za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i
     .2,i4)))
      
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(za(i1,i5)**2*zb(i1,i4)*zb(i2,i4)*
     .((-2._dp*L0(-s(i5,i6),-t(i2,i3,i4))*zb(i2,i4))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .(L1(-s(i5,i6),-t(i2,i3,i4))*zb(i1,i4))/t(i2,i3,i4)))/
     .(2._dp*za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i
     .4))*
     .zb(i3,i4))+(za(i3,i5)*zb(i2,i4)*
     .(za(i3,i5)/za(i3,i4)+(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))/
     .t(i2,i3,i4)))/
     .(2._dp*za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i
     .4)))-
     .(Lsm1_2mh(s(i1,i2),t(i2,i3,i4),s(i3,i4),s(i5,i6))*za(i1,i5)**2*zb(
     .i2,i4)**3*
     .t(i2,i3,i4))/
     .(za(i5,i6)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     .*3*
     .zb(i3,i4))+(Lnrat(-s(i1,i2),-s(i5,i6))*zb(i2,i4)*
     .(((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i2)*
     .(-(za(i2,i5)*(-(za(i1,i5)*za(i2,i3))+za(i1,i2)*za(i3,i5))*
     .zb(i2,i4))+(za(i1,i3)*za(i1,i5)**2*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i2,i4)
     .)/(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))-
     .za(i1,i5)*za(i2,i5)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)
     .)
     .-(2._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i1,i2)*za(i3,i5)*
     .(zb(i1,i6)*zb(i2,i4)+zb(i1,i4)*zb(i2,i6)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     .(za(i1,i3)*zb(i2,i4)*(-((za(i1,i5)*za(i2,i5)*zb(i1,i2))/
     .za(i5,i6))-(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))-
     .(za(i2,i3)*zb(i2,i6)*zb(i4,i6))/zb(i5,i6)+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*
     .(zb(i1,i4)*zb(i2,i6)*
     .(-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))-
     .za(i1,i3)*zb(i2,i4)*
     .(zb(i1,i6)**2+((-(za(i2,i3)*zb(i1,i3))-
     .za(i2,i4)*zb(i1,i4))*zb(i2,i6)**2)/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     .zb(i1,i2)*zb(i4,i6)*
     .(-(za(i2,i3)*zb(i2,i6))+za(i3,i4)*zb(i4,i6))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i5,i6)
     .)
     .+(6._dp*s(i1,i2)*s(i3,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4
     .))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(t(i1,i3,i4)-t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2-
     .(s(i1,i2)*((za(i3,i5)**2*zb(i3,i4))/za(i5,i6)+
     .(za(i3,i4)*zb(i4,i6)**2)/zb(i5,i6))*t(i2,i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*s(i3,i4)*zb(i2,i3)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4
     .)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(Lnrat(-s(i1,i2),-s(i5,i6))*za(i1,i3)*
     .((za(i1,i5)*za(i3,i5)*zb(i1,i4))/za(i5,i6)+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i2)*
     .(za(i1,i5)*za(i2,i3)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))-
     .za(i1,i3)*zb(i2,i4)*
     .(za(i2,i5)**2+(za(i1,i5)**2*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))-
     .za(i1,i2)*za(i3,i5)*
     .(-(za(i1,i5)*zb(i1,i4))-za(i3,i5)*zb(i3,i4))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)
     .)
     .+(2._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*
     .(za(i1,i5)*za(i2,i3)+za(i1,i3)*za(i2,i5))*zb(i1,i2)*zb(i4,i6)
     .)/(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5
     .,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     .(za(i1,i3)*zb(i2,i4)*(-((za(i1,i5)*za(i2,i5)*zb(i1,i2))/
     .za(i5,i6))-(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6)))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*
     .(-(zb(i1,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .zb(i2,i6))+(za(i1,i3)*
     .(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .zb(i2,i4)*zb(i2,i6)**2)/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))-
     .za(i1,i3)*zb(i1,i6)*
     .(-(zb(i1,i4)*zb(i2,i6))-zb(i1,i2)*zb(i4,i6))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i5,i6)
     .)
     .-(s(i1,i2)*((za(i3,i5)**2*zb(i3,i4))/za(i5,i6)+
     .(za(i3,i4)*zb(i4,i6)**2)/zb(i5,i6))*t(i1,i3,i4))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     .(6._dp*s(i1,i2)*s(i3,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)
     .)*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(-t(i1,i3,i4)+t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2))/
     .(2._dp*s(i3,i4)*za(i1,i4)*(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4
     .)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(I3m(s(i1,i2),s(i3,i4),s(i5,i6))*za(i1,i3)*
     .(-((za(i1,i3)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     .zb(i2,i4)*(-(za(i1,i3)*
     .(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)))-
     .za(i1,i4)*(-(za(i1,i5)*zb(i1,i4))-
     .za(i2,i5)*zb(i2,i4))+
     .za(i1,i2)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))*
     .zb(i2,i6))/(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))+
     .zb(i1,i4)*(za(i1,i3)*(-(za(i1,i2)*
     .(za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4)))-
     .za(i2,i3)*(-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))-
     .za(i2,i4)*(-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))+
     .za(i1,i2)*(-(za(i1,i3)*
     .(za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4)))-
     .za(i2,i3)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))))*
     .zb(i2,i6)+(s(i1,i2)*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6)))/2+
     .za(i2,i3)*zb(i1,i2)*(-(za(i1,i2)*za(i3,i5)*zb(i2,i6)*zb(i3,i4))+
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     .(3._dp*s(i1,i2)*(s(i1,i2)-s(i3,i4)-s(i5,i6))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(-t(i1,i3,i4)+t(i2,i3,i4)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     .2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2))))/
     .((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     .i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i4)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))
     
      zaj_virt_a6_slc_fsc_mp=zaj_virt_a6_slc_fsc_mp+
     .(Lnrat(-s(i3,i4),-s(i5,i6))*za(i1,i3)*
     .(-2._dp*(za(i1,i3)*zb(i1,i4)-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))-
     .(8._dp*s(i1,i2)*za(i1,i3)*za(i5,i6)*zb(i2,i6)*zb(i4,i6))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))+
     .4._dp*(-(za(i1,i2)*za(i3,i5)*zb(i1,i6)*zb(i2,i4))+
     .za(i1,i5)*za(i2,i3)*zb(i1,i2)*zb(i4,i6))+
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)+
     .(4._dp*s(i5,i6)*za(i1,i3)*zb(i2,i4))/
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))*
     .((za(i1,i5)*za(i2,i5)*zb(i1,i2))/za(i5,i6)-
     .(za(i1,i2)*zb(i1,i6)*zb(i2,i6))/zb(i5,i6))-
     .(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     .((za(i1,i5)*za(i3,i5)*zb(i1,i4))/za(i5,i6)-
     .(za(i2,i3)*zb(i2,i6)*zb(i4,i6))/zb(i5,i6)-
     .(2._dp*za(i1,i3)*zb(i2,i6)*
     .(-(za(i1,i2)*zb(i1,i4)*zb(i2,i6))+
     .za(i2,i3)*zb(i2,i6)*zb(i3,i4)+
     .za(i1,i5)*zb(i1,i4)*zb(i5,i6)))/
     .((-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*zb(i5,i6)))+
     .(3._dp*(-s(i1,i2)-s(i3,i4)+s(i5,i6))*
     .(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     .(-(za(i1,i5)*zb(i1,i6))-za(i2,i5)*zb(i2,i6))*
     .(-t(i1,i3,i4)+t(i2,i3,i4)))/
     .(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     .6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)))/
     .(2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     .(i5,i6)-
     .2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i4)*
     .(-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4)))


      end function

c     primitive amplitudes for q;1 leading color piece
c     !!  note that the argument ordering is changed to match
c     !!  q,qb,g,a,lb,l for the leading color pieces
      function zaj_virt_a6_lc_pole_pm(i1,i3,i2,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_tree_pm
        complex(dp) :: zaj_virt_a6_lc_pole_pm
        complex(dp) :: Vcc, Vsc

        complex(dp) :: xl12, xl23, xl56
        complex(dp) :: lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)


        xl12=Lnrat(musq,real(-s(i1,i2),dp))
        xl23=Lnrat(musq,real(-s(i2,i3),dp))
        xl56=Lnrat(musq,real(-s(i5,i6),dp))

        zaj_virt_a6_lc_tree_pm = 
     &    (((za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4))*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)))/
     &    (s(i5,i6)*za(i1,i2)*za(i2,i3)*zb(i1,i4)*zb(i3,i4))+
     &    (za(i3,i5)*zb(i1,i2)*(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6)))
     &    /
     &    (s(i5,i6)*za(i1,i2)*zb(i1,i4)*t(i1,i2,i4))-
     &    (za(i3,i4)*zb(i1,i6)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/
     &    (s(i5,i6)*za(i2,i3)*zb(i3,i4)*t(i2,i3,i4)))

        Vcc=-four
     &   -(+(epinv*epinv2+xl12*epinv+half*xl12**2)
     &     +(epinv*epinv2+xl23*epinv+half*xl23**2))-two*(epinv+xl56)
        Vsc=half*(one+epinv+xl56)

        zaj_virt_a6_lc_pole_pm = zaj_virt_a6_lc_tree_pm*(Vcc+Vsc)

      end function

      ! this helper routine makes it simpler to use the same
      ! crossing code from the LO amplitudes
      function zaj_virt_a6_lc_pole_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_pole_mp

        zaj_virt_a6_lc_pole_mp =
     &    zaj_virt_a6_lc_pole_pm(i2,i1,i3,i4,i6,i5,zb,za)

      end function

      function zaj_virt_a6_lc_pole_pp(i1,i3,i2,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_tree_pp
        complex(dp) :: zaj_virt_a6_lc_pole_pp
        complex(dp) :: Vcc, Vsc

        complex(dp) :: xl12, xl23, xl56
        complex(dp) :: lnrat

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)


        xl12=Lnrat(musq,real(-s(i1,i2),dp))
        xl23=Lnrat(musq,real(-s(i2,i3),dp))
        xl56=Lnrat(musq,real(-s(i5,i6),dp))

        zaj_virt_a6_lc_tree_pp = 
     &    (-za(i1,i3)*za(i3,i5)**2)/
     &    (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))

        Vcc=-four
     &   -(epinv*epinv2+xl12*epinv+half*xl12**2)
     &   -(epinv*epinv2+xl23*epinv+half*xl23**2)-two*(epinv+xl56)
        Vsc=half*(one+epinv+xl56)

        zaj_virt_a6_lc_pole_pp = zaj_virt_a6_lc_tree_pp*(Vcc+Vsc)

      end function

      function zaj_virt_a6_lc_fcc_pp(i1,i3,i2,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_fcc_pp

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        complex(dp) :: L0, Lsm1, Lsm1_2me


        zaj_virt_a6_lc_fcc_pp =
     &    (Lsm1(-s(i1,i4),-t(i1,i2,i4),-s(i1,i2),-t(i1,i2,i4))*za(i3,i5)**2)
     &    /
     &    (za(i1,i4)*za(i2,i3)*za(i2,i4)*za(i5,i6))+
     &    (Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))*za(i1,i3)*za(
     &    i3,i5)**2)/
     &    (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1_2me(t(i1,i2,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*za(i1,i3)*
     &    za(i3,i5)**2)/(za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1_2me(t(i1,i2,i3),t(i2,i3,i4),s(i2,i3),s(i5,i6))*za(i1,i3)*za(
     &    i3,i5)*
     &    (-(za(i1,i5)*za(i3,i4))+za(i1,i3)*za(i4,i5)))/
     &    (za(i1,i2)*za(i1,i4)**2*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (Lsm1(-s(i2,i3),-t(i2,i3,i4),-s(i3,i4),-t(i2,i3,i4))*za(i3,i5)*
     &    (-(za(i2,i5)*za(i3,i4))+za(i2,i3)*za(i4,i5)))/
     &    (za(i1,i2)*za(i2,i4)**2*za(i3,i4)*za(i5,i6))-
     &    (2._dp*L0(-t(i2,i3,i4),-s(i5,i6))*za(i1,i3)*za(i1,i5)*za(i3,i5)*
     &    (za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4)))/
     &    (s(i5,i6)*za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &    (2._dp*L0(-t(i2,i3,i4),-s(i3,i4))*za(i2,i5)*za(i3,i5)*zb(i2,i4))/
     &    (s(i3,i4)*za(i1,i2)*za(i2,i4)*za(i5,i6))+
     &    (2._dp*L0(-t(i2,i3,i4),-s(i2,i3))*za(i3,i5)*za(i4,i5)*zb(i2,i4))/
     &    (s(i2,i3)*za(i1,i4)*za(i2,i4)*za(i5,i6))

      end function

      function zaj_virt_a6_lc_fcc_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_lc_fcc_mp

        zaj_virt_a6_lc_fcc_mp = zaj_virt_a6_lc_fcc_pm(i2,i1,i3,i4,i6,i5,zb,za)

      end function

      function zaj_virt_a6_lc_fcc_pm(i1,i3,i2,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_fcc_pm

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        complex(dp) :: L0,Lsm1,Lsm1_2mh,Lnrat,I3m


        zaj_virt_a6_lc_fcc_pm =
     &    (2._dp*L0(-t(i2,i3,i4),-s(i2,i3))*za(i4,i5)*zb(i2,i4)*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/
     &    (s(i2,i3)*za(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3
     &    ,i4)))
     &    +(Lsm1(-s(i1,i2),-t(i1,i2,i3),-s(i2,i3),-t(i1,i2,i3))*za(i1,i3)*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))**2)/
     &    (za(i1,i2)*za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))+
     &    (Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*za(i1,i3)*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))**2)/
     &    (za(i1,i2)*za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))-
     &    (2._dp*L0(-t(i1,i2,i4),-s(i5,i6))*za(i3,i4)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*zb(i3,i6))/
     &    (s(i5,i6)*za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    zb(i5,i6))+Lsm1_2mh(s(i1,i4),t(i1,i2,i3),s(i2,i3),s(i5,i6))*
     &    (-((za(i1,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     &    (za(i1,i2)*za(i2,i3)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**3*
     &    zb(i5,i6)))+(za(i1,i3)**3*zb(i4,i6)**2*t(i1,i2,i3)**2)/
     &    (za(i1,i2)*za(i2,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**3*zb(i5,i6)))-
     &    (2._dp*Lnrat(-s(i2,i3),-s(i5,i6))*(2._dp*(s(i1,i3)-s(i2,i4))*za(i3,i4)
     &    *zb(i2,i3)*
     &    (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6))+
     &    (za(i3,i5)*zb(i1,i2)*(((-s(i1,i4)+s(i2,i3)-s(i5,i6))*
     &    za(i3,i5)*zb(i2,i3)+
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*zb(i2,i6))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))+
     &    (s(i1,i3)-s(i2,i4))*
     &    (-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i5)*zb(i2,i3))+
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*zb(i3,i6))))/
     &    (za(i5,i6)*zb(i1,i4))-
     &    (za(i3,i4)*(zb(i1,i2)*
     &    (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))-
     &    zb(i2,i4)*(-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6))+
     &    zb(i2,i3)*(-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))*
     &    (-(za(i2,i5)*zb(i2,i3)*zb(i5,i6))+zb(i3,i6)*t(i1,i2,i3)))/
     &    zb(i5,i6)))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))

          zaj_virt_a6_lc_fcc_pm=zaj_virt_a6_lc_fcc_pm-
     &    (Lsm1(-s(i1,i4),-t(i1,i2,i4),-s(i1,i2),-t(i1,i2,i4))*za(i3,i5)**2*
     &    zb(i1,i2)**2)/
     &    (za(i5,i6)*zb(i1,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    t(i1,i2,i4))-(2._dp*L0(-t(i1,i2,i4),-s(i1,i4))*za(i1,i4)*zb(i1,i2)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    (s(i1,i4)*za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    zb(i5,i6)*t(i1,i2,i4))-Lsm1_2mh(s(i2,i3),t(i1,i2,i4),s(i1,i4),
     &    s(i5,i6))*((za(i3,i5)**2*zb(i1,i2)**2)/
     &    (za(i5,i6)*zb(i1,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    t(i1,i2,i4))+((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))**2)/
     &    (za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))**3*
     &    zb(i5,i6)*t(i1,i2,i4))-
     &    (za(i2,i4)**2*zb(i3,i6)**2*t(i1,i2,i4))/
     &    (za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))**3*zb(i5,i
     &    6))
     &    )-(I3m(s(i1,i4),s(i2,i3),s(i5,i6))*zb(i1,i2)*
     &    (-(s(i5,i6)*za(i3,i4)*((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*
     &    zb(i3,i6)+2._dp*za(i1,i4)*za(i2,i5)*zb(i2,i3)*zb(i4,i6)))+
     &    2._dp*za(i1,i4)*za(i1,i5)*za(i3,i5)*
     &    ((-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i1,i3)+
     &    2._dp*za(i2,i4)*zb(i1,i4)*zb(i2,i3))*zb(i5,i6)-
     &    2._dp*za(i1,i4)*za(i4,i5)*zb(i4,i6)*
     &    (s(i5,i6)*(-s(i1,i4)-s(i2,i3)+s(i5,i6))+
     &    (s(i1,i4)-s(i2,i3)-s(i5,i6))*t(i1,i2,i4))))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))
     
          zaj_virt_a6_lc_fcc_pm=zaj_virt_a6_lc_fcc_pm+
     &    I3m(s(i1,i4),s(i2,i3),s(i5,i6))*
     &    ((-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i3,i5)*zb(i1,i2)*
     &    (za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4)))/
     &    (2._dp*za(i5,i6)*zb(i1,i4))-
     &    za(i3,i4)*zb(i1,i6)*(-(za(i1,i5)*zb(i1,i2))+
     &    za(i3,i5)*zb(i2,i3))+
     &    (za(i2,i4)*za(i3,i5)*zb(i1,i2)-za(i3,i4)*za(i4,i5)*zb(i1,i4))*
     &    zb(i2,i6)+(za(i1,i4)*za(i2,i3)*za(i3,i5)*zb(i1,i2)*zb(i3,i6))/
     &    za(i1,i2)+(za(i1,i3)*za(i4,i5)*zb(i1,i2)*
     &    ((s(i1,i3)-s(i2,i4))*za(i1,i2)+
     &    za(i1,i3)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))*
     &    zb(i4,i6))/
     &    (za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))-
     &    ((-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*
     &    (-2._dp*s(i5,i6)*za(i1,i3)*zb(i1,i6)+
     &    (-s(i1,i4)+s(i2,i3)-s(i5,i6))*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))))/
     &    (2._dp*za(i1,i2)*zb(i5,i6)))/
     &    ((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))+
     &    (za(i1,i3)*za(i4,i5)*zb(i1,i2)*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)))/
     &    (za(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*t(i1,i2,i3
     &    ))
     &    -(2._dp*za(i1,i4)*za(i2,i4)*za(i3,i5)*zb(i1,i2)*zb(i3,i6))/
     &    (za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*t(i1,i2,i4
     &    ))
     &    +(((za(i3,i5)**2*zb(i1,i2)**2)/
     &    (za(i5,i6)*zb(i1,i4)*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     &    (-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))**2/
     &    (za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    zb(i5,i6)))*(2._dp*s(i1,i4)*s(i5,i6)+
     &    (-s(i1,i4)+s(i2,i3)-s(i5,i6))*t(i1,i2,i4)))/(2._dp*t(i1,i2,i4)**2))-
     &    (2._dp*Lnrat(-s(i1,i4),-s(i5,i6))*(-2._dp*(s(i1,i3)-s(i2,i4))*za(i1,i4
     &    )*zb(i1,i2)*
     &    (-(za(i1,i5)*zb(i1,i6))-za(i4,i5)*zb(i4,i6))+
     &    (za(i3,i4)*zb(i1,i6)*((s(i1,i3)-s(i2,i4))*
     &    ((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i4,i6)-
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*zb(i5,i6))+
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    (-((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i1,i6))-
     &    (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i4,i5)*zb(i5,i6))))/
     &    (za(i2,i3)*zb(i5,i6))+
     &    (zb(i1,i2)*(-(za(i1,i4)*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3)))+
     &    za(i2,i4)*(-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))+
     &    za(i3,i4)*(-(za(i1,i5)*zb(i1,i3))+za(i4,i5)*zb(i3,i4)))*
     &    (-(za(i1,i4)*za(i5,i6)*zb(i4,i6))+za(i1,i5)*t(i1,i3,i4)))/
     &    za(i5,i6)))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))
     
          zaj_virt_a6_lc_fcc_pm=zaj_virt_a6_lc_fcc_pm+
     &    Lsm1_2mh(s(i1,i2),t(i2,i3,i4),s(i3,i4),s(i5,i6))*
     &    ((za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2/
     &    (za(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    t(i2,i3,i4))-(za(i3,i4)**2*zb(i1,i6)**2)/
     &    (za(i2,i3)*(-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i5,i6)*
     &    t(i2,i3,i4)))+I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    ((2._dp*za(i3,i4)**2*zb(i1,i6)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,
     &    i4)))/
     &    (za(i2,i3)*t(i2,i3,i4)**2)+
     &    (za(i1,i3)*(-(za(i3,i4)*za(i4,i5)*zb(i1,i2)*zb(i1,i6))+
     &    ((za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)))/
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     &    (za(i2,i3)*t(i1,i2,i3)*t(i2,i3,i4)))+
     &    (2._dp*L0(-t(i2,i3,i4),-s(i5,i6))*zb(i1,i6)*
     &    (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/
     &    (s(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    t(i2,i3,i4))-(2._dp*Lnrat(-s(i2,i3),-s(i5,i6))*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2)/
     &    (za(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    t(i2,i3,i4))+(Lsm1(-s(i2,i3),-t(i2,i3,i4),-s(i3,i4),-t(i2,i3,i4))*
     &    (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4))**2)/
     &    (za(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    t(i2,i3,i4))-(I3m(s(i2,i3),s(i1,i4),s(i5,i6))*za(i3,i4)*
     &    (2._dp*za(i5,i6)*zb(i1,i6)*zb(i2,i3)*
     &    (-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i3))-
     &    2._dp*za(i1,i4)*za(i2,i3)*zb(i2,i4))*zb(i3,i6)-
     &    s(i5,i6)*zb(i1,i2)*((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*
     &    zb(i3,i6)+2._dp*za(i1,i4)*za(i2,i5)*zb(i2,i3)*zb(i4,i6))+
     &    2._dp*za(i2,i5)*zb(i2,i3)*zb(i2,i6)*
     &    (s(i5,i6)*(-s(i1,i4)-s(i2,i3)+s(i5,i6))+
     &    (-s(i1,i4)+s(i2,i3)-s(i5,i6))*t(i2,i3,i4))))/
     &    ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &    i6)-
     &    2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &    (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))

      end function

      function zaj_virt_a6_lc_fsc_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a6_lc_fsc_mp

        zaj_virt_a6_lc_fsc_mp = zaj_virt_a6_lc_fsc_pm(i2,i1,i3,i4,i6,i5,zb,za)

      end function

      function zaj_virt_a6_lc_fsc_pm(i1,i3,i2,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_fsc_pm
        complex(dp) :: L0,L1,Lsm1_2mh,I3m,Lnrat
        complex(dp) :: Fsc4, fac

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        Fsc4= 
     &  -(Lnrat(-s(i2,i3),-s(i5,i6))*zb(i1,i2)*
     &  ((za(i1,i5)*za(i3,i5)*zb(i1,i3)+za(i1,i5)*za(i4,i5)*zb(i1,i4))/
     &  za(i5,i6)-(-(za(i1,i4)*zb(i1,i6)*zb(i4,i6))-
     &  za(i2,i4)*zb(i2,i6)*zb(i4,i6))/zb(i5,i6)))/
     &  (2._dp*zb(i1,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))

        fac=
     &  -(za(i2,i4)*((za(i4,i5)*
     &  (-((za(i3,i5)*
     &  ((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*
     &  zb(i1,i3)-
     &  (-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i4)*
     &  zb(i2,i3)))/za(i5,i6))+
     &  ((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i2,i4)+
     &  2._dp*za(i1,i4)*za(i2,i3)*zb(i1,i3))*zb(i2,i6)))/
     &  za(i1,i4)-(zb(i1,i6)*
     &  (-(za(i3,i5)*
     &  (-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i1,i3))-
     &  2._dp*za(i2,i4)*zb(i1,i4)*zb(i2,i3)))-
     &  ((-((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*
     &  zb(i1,i3))+
     &  (s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i2,i4)*
     &  zb(i1,i4))*zb(i2,i6))/zb(i5,i6)))/zb(i1,i4)))/
     &  (2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))+
     &  (za(i2,i4)*zb(i1,i2)*(-(za(i1,i2)*zb(i1,i6))+
     &  za(i2,i4)*zb(i4,i6))*
     &  ((L0(-t(i1,i2,i4),-s(i1,i4))*
     &  ((-2._dp*za(i2,i4)*zb(i3,i6))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))+
     &  (-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))/
     &  t(i1,i2,i4)))/s(i1,i4)-
     &  (L1(-s(i1,i4),-t(i1,i2,i4))*za(i2,i4)*zb(i2,i6))/t(i1,i2,i4)**2))/
     &  (2._dp*za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &  zb(i5,i6))
        fac=fac+(za(i2,i4)*za(i3,i4)*zb(i3,i6)**2*
     &  ((-2._dp*L0(-s(i5,i6),-t(i1,i2,i4))*za(i2,i4))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))+
     &  (L1(-s(i5,i6),-t(i1,i2,i4))*za(i3,i4))/t(i1,i2,i4)))/
     &  (2._dp*za(i1,i2)*za(i1,i4)*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i5,i6))-
     &  (Lsm1_2mh(s(i2,i3),t(i1,i2,i4),s(i1,i4),s(i5,i6))*za(i2,i4)**3*
     &  zb(i3,i6)**2*t(i1,i2,i4))/
     &  (za(i1,i2)*za(i1,i4)*(-(za(i1,i2)*zb(i1,i3))-
     &  za(i2,i4)*zb(i3,i4))**3*zb(i5,i6))-
     &  (Lnrat(-s(i2,i3),-s(i5,i6))*za(i2,i4)*
     &  (-((za(i2,i5)*za(i4,i5)*zb(i1,i2))/za(i5,i6))+
     &  (2._dp*(-s(i1,i4)-s(i2,i3)+s(i5,i6))*
     &  (za(i2,i5)*za(i3,i4)+za(i2,i4)*za(i3,i5))*zb(i1,i6)*
     &  zb(i2,i3))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)+
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i2,i3)*
     &  (za(i2,i5)*za(i3,i4)*
     &  (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))-
     &  za(i2,i3)*za(i4,i5)*
     &  (za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))+
     &  za(i2,i4)*zb(i1,i3)*
     &  (za(i3,i5)**2+
     &  (za(i2,i5)**2*
     &  (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))))/
     &  ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &  za(i5,i6))+(za(i2,i4)*zb(i1,i3)*
     &  (-((za(i2,i5)*za(i3,i5)*zb(i2,i3))/za(i5,i6))-
     &  (za(i2,i3)*zb(i2,i6)*zb(i3,i6))/zb(i5,i6)))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))+
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*
     &  (-((za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  zb(i2,i6)*zb(i3,i6))-
     &  (za(i2,i4)*zb(i1,i3)*
     &  (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &  zb(i3,i6)**2)/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))-
     &  za(i2,i4)*zb(i2,i6)*
     &  (-(zb(i1,i6)*zb(i2,i3))+zb(i1,i2)*zb(i3,i6))))/
     &  ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*
     &  zb(i5,i6))-(s(i2,i3)*
     &  (-((za(i4,i5)**2*zb(i1,i4))/za(i5,i6))-
     &  (za(i1,i4)*zb(i1,i6)**2)/zb(i5,i6))*t(i1,i2,i4))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)+
     &  (6._dp*s(i1,i4)*s(i2,i3)*
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6))*
     &  (-t(i1,i2,i4)+t(i1,i3,i4)))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2))/
     &  (2._dp*s(i1,i4)*za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-
     &  za(i2,i4)*zb(i3,i4)))
        fac=fac-
     &  (I3m(s(i2,i3),s(i1,i4),s(i5,i6))*za(i2,i4)*
     &  ((za(i2,i4)*zb(i1,i3)*
     &  (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4))*
     &  (za(i1,i2)*(za(i2,i5)*zb(i1,i2)+
     &  za(i3,i5)*zb(i1,i3))-
     &  za(i2,i4)*(-(za(i2,i5)*zb(i2,i4))-
     &  za(i3,i5)*zb(i3,i4))+
     &  za(i2,i3)*(-(za(i1,i5)*zb(i1,i3))+
     &  za(i4,i5)*zb(i3,i4)))*zb(i3,i6))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))-
     &  zb(i1,i2)*(za(i2,i4)*
     &  (za(i1,i3)*(za(i2,i5)*zb(i1,i2)+
     &  za(i3,i5)*zb(i1,i3))-
     &  za(i2,i3)*
     &  (-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4))-
     &  za(i3,i4)*
     &  (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))+
     &  za(i2,i3)*(-(za(i2,i4)*
     &  (-(za(i1,i5)*zb(i1,i2))+za(i4,i5)*zb(i2,i4)))-
     &  za(i3,i4)*
     &  (-(za(i1,i5)*zb(i1,i3))+za(i4,i5)*zb(i3,i4))))*
     &  zb(i3,i6)+(s(i2,i3)*
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/2+
     &  za(i3,i4)*zb(i2,i3)*
     &  (za(i2,i3)*za(i4,i5)*zb(i1,i4)*zb(i3,i6)+
     &  (-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))+
     &  (3._dp*s(i2,i3)*(-s(i1,i4)+s(i2,i3)-s(i5,i6))*
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6))*
     &  (-t(i1,i2,i4)+t(i1,i3,i4)))/
     &  (2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2))))/
     &  ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))
        fac=fac-
     &  (Lnrat(-s(i1,i4),-s(i5,i6))*za(i2,i4)*
     &  (4._dp*(za(i2,i5)*za(i3,i4)*zb(i1,i6)*zb(i2,i3)+
     &  za(i2,i3)*za(i4,i5)*zb(i1,i3)*zb(i2,i6))-
     &  (8._dp*s(i2,i3)*za(i2,i4)*za(i5,i6)*zb(i1,i6)*zb(i3,i6))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))-
     &  2._dp*(-(za(i2,i4)*zb(i1,i2))+za(i3,i4)*zb(i1,i3))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6))+
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3)-
     &  (4._dp*s(i5,i6)*za(i2,i4)*zb(i1,i3))/
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))*
     &  ((za(i2,i5)*za(i3,i5)*zb(i2,i3))/za(i5,i6)-
     &  (za(i2,i3)*zb(i2,i6)*zb(i3,i6))/zb(i5,i6))-
     &  (s(i1,i4)-s(i2,i3)-s(i5,i6))*
     &  (-((za(i2,i5)*za(i4,i5)*zb(i1,i2))/za(i5,i6))-
     &  (za(i3,i4)*zb(i1,i6)*zb(i3,i6))/zb(i5,i6)-
     &  (2._dp*za(i2,i4)*zb(i3,i6)*
     &  (za(i2,i3)*zb(i1,i2)*zb(i3,i6)-
     &  za(i3,i4)*zb(i1,i4)*zb(i3,i6)-
     &  za(i2,i5)*zb(i1,i2)*zb(i5,i6)))/
     &  ((-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &  zb(i5,i6)))+
     &  (3._dp*(-s(i1,i4)-s(i2,i3)+s(i5,i6))*
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6))*
     &  (-t(i1,i2,i4)+t(i1,i3,i4)))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     &  (2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))

        Fsc4=Fsc4+za(i1,i4)/za(i2,i4)*fac

        Fsc4=Fsc4-
     &  (za(i1,i5)*zb(i1,i2)*((L0(-t(i2,i3,i4),-s(i5,i6))*
     &  (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/(s(i5,i6)*za(i5,i6))
     &  +(L1(-s(i5,i6),-t(i2,i3,i4))*(-(za(i1,i3)*zb(i1,i6)*zb(i2,i3))-
     &  za(i1,i4)*zb(i1,i6)*zb(i2,i4)))/t(i2,i3,i4)**2))/
     &  (2._dp*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     &  (za(i3,i4)*zb(i2,i3)*(-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))*
     &  ((L1(-t(i2,i3,i4),-s(i2,i3))*za(i4,i5)*zb(i2,i4))/s(i2,i3)**2+
     &  (L0(-t(i2,i3,i4),-s(i2,i3))*
     &  (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/(s(i2,i3)*t(i2,i3,i4))
     &  ))/(2._dp*za(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i
     &  4)))+
     &  (-((za(i1,i3)*za(i3,i5)*za(i4,i5))/
     &  (za(i1,i2)*za(i2,i3)*za(i5,i6)*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))-
     &  (za(i1,i5)*za(i3,i5)*zb(i1,i2))/
     &  (za(i2,i3)*za(i5,i6)*zb(i3,i4)*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))-
     &  (za(i3,i5)**2*zb(i1,i3)**2)/
     &  (za(i2,i3)*za(i5,i6)*zb(i1,i4)*zb(i3,i4)*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))-
     &  (s(i3,i4)*za(i1,i5)*za(i4,i5))/
     &  (za(i1,i2)*za(i5,i6)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))+
     &  (za(i1,i5)*za(i3,i4)*za(i3,i5)*zb(i1,i3))/
     &  (za(i2,i3)*za(i5,i6)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4)))-
     &  (za(i1,i3)*zb(i1,i6)**2)/
     &  (za(i1,i2)*za(i2,i3)*zb(i1,i4)*zb(i3,i4)*zb(i5,i6))-
     &  (za(i1,i3)*zb(i1,i2)*zb(i1,i6)*zb(i4,i6))/
     &  (za(i2,i3)*zb(i1,i4)*zb(i3,i4)*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))+
     &  (za(i1,i4)*za(i3,i4)*zb(i3,i6)*zb(i4,i6))/
     &  (za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i5,i6))-
     &  (s(i1,i3)*za(i3,i4)*zb(i3,i6)*zb(i4,i6))/
     &  (za(i2,i3)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i5,i6))+
     &  (za(i3,i4)*(za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))*zb(i4,i6))/
     &  (za(i1,i2)*za(i2,i3)*zb(i3,i4)*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*zb(i5,i6))-
     &  (za(i2,i4)*zb(i2,i6)*(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6)))
     &  /
     &  (za(i1,i2)*(-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*zb(i5,i6)*
     &  t(i1,i2,i4))-(za(i1,i5)*zb(i1,i2)*
     &  (za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/
     &  (za(i5,i6)*zb(i3,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  t(i2,i3,i4)))/2._dp
     
        fac=
     &  (za(i4,i5)*zb(i1,i6)*(-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i3))-
     &  2._dp*za(i1,i4)*za(i2,i3)*zb(i2,i4)))/
     &  ((s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-2._dp*s(i1,i4)*s(i5,
     &  i6)-
     &  2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     &  (((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i3)*zb(i1,i4)-
     &  (-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*zb(i2,i4))*
     &  (-((za(i4,i5)**2*zb(i1,i4))/za(i5,i6))-
     &  (za(i1,i4)*zb(i1,i6)**2)/zb(i5,i6)))/
     &  (2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  zb(i1,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))+
     &  (Lsm1_2mh(s(i1,i4),t(i1,i2,i3),s(i2,i3),s(i5,i6))*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &  (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     &  (za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**3*zb(i5,i6))
     &  -(L0(-t(i1,i2,i3),-s(i2,i3))*za(i2,i3)*zb(i1,i2)*
     &  (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))**2)/
     &  (s(i2,i3)*za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**2*
     &  zb(i5,i6))
        fac=fac+(Lnrat(-s(i2,i3),-s(i5,i6))*za(i2,i3)*
     &  ((-6._dp*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*zb(i2,i3)*
     &  (-((-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i3))-
     &  2._dp*za(i1,i4)*za(i2,i3)*zb(i2,i4))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2-
     &  ((za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))*
     &  (zb(i1,i4)*zb(i2,i6)+zb(i1,i2)*zb(i4,i6)))/
     &  (zb(i1,i4)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  zb(i5,i6))+(-((za(i4,i5)*
     &  (za(i4,i5)*zb(i1,i4)*
     &  (-(za(i1,i3)*zb(i2,i3))-za(i1,i4)*zb(i2,i4))
     &  +2._dp*(-((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i5)*zb(i1,i2))+
     &  (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i5,i6)*
     &  zb(i2,i6))))/za(i5,i6))+
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i3)*zb(i2,i3)*
     &  (-(za(i1,i5)*zb(i1,i6))+za(i4,i5)*zb(i4,i6)))/
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))-
     &  (za(i3,i4)*zb(i1,i6)*zb(i2,i3)*
     &  (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6)))/zb(i5,i6)
     &  -(zb(i2,i3)*zb(i4,i6)*(-((2._dp*(-s(i1,i4)+s(i2,i3)-s(i5,i6))*
     &  za(i3,i4)-za(i1,i4)*za(i2,i3)*zb(i1,i2))*
     &  zb(i1,i6))-
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i3)*
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  zb(i4,i6))/
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     &  (zb(i1,i4)*zb(i5,i6))+
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*zb(i1,i6)*zb(i2,i4)*
     &  (-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6)-
     &  za(i4,i5)*zb(i5,i6)))/(zb(i1,i4)*zb(i5,i6))+
     &  za(i4,i5)*(za(i1,i3)*zb(i1,i6)*zb(i2,i3)+
     &  za(i4,i5)*zb(i2,i4)*zb(i5,i6))+
     &  4._dp*za(i1,i5)*zb(i2,i3)*
     &  (za(i3,i4)*zb(i1,i6)+
     &  (za(i1,i3)*za(i4,i5)*zb(i1,i4)*zb(i5,i6))/
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     &  (2._dp*za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
        fac=fac-
     &  (L1(-s(i5,i6),-t(i1,i2,i3))*za(i3,i4)*za(i4,i5)*zb(i4,i6))/
     &  (za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*t(i1,i2,i3))
     &  -(L0(-t(i1,i2,i3),-s(i5,i6))*za(i3,i4)*
     &  (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6))*zb(i4,i6)*t(i1,i2,i3)
     &  )/(s(i5,i6)*za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))**2
     &  *
     &  zb(i5,i6))+(Lnrat(-s(i1,i4),-s(i5,i6))*
     &  ((-6._dp*za(i1,i4)*za(i2,i3)*
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  (2._dp*za(i1,i3)*zb(i1,i4)*zb(i2,i3)+
     &  (-s(i1,i4)-s(i2,i3)+s(i5,i6))*zb(i2,i4))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)+
     &  (za(i4,i5)*(-3._dp*(s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i3)*
     &  za(i4,i5)*zb(i1,i4)+
     &  za(i1,i3)*za(i1,i4)*
     &  (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*zb(i1,i4)
     &  +za(i1,i4)*za(i2,i3)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &  zb(i2,i4)+
     &  2._dp*za(i2,i3)*zb(i1,i2)*
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i5)+
     &  2._dp*za(i1,i4)*za(i5,i6)*zb(i4,i6))))/za(i5,i6)+
     &  (za(i1,i4)*zb(i1,i6)*
     &  ((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &  (-(za(i2,i4)*zb(i2,i6))-za(i3,i4)*zb(i3,i6))+
     &  2._dp*((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i3,i4)*zb(i4,i6)-
     &  (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i3,i5)*zb(i5,i6)))
     &  )/zb(i5,i6)-(za(i1,i3)*((za(i1,i5)*za(i4,i5)*zb(i1,i4))/za(i5,i6)+
     &  (za(i1,i4)*zb(i1,i6)*zb(i4,i6))/zb(i5,i6))*
     &  (-s(i1,i4)**2+2._dp*s(i1,i4)*s(i2,i3)-s(i2,i3)**2+
     &  2._dp*s(i1,i4)*s(i5,i6)+2._dp*s(i2,i3)*s(i5,i6)-s(i5,i6)**2+
     &  (s(i1,i4)-s(i2,i3)-s(i5,i6))*(t(i1,i2,i3)-t(i2,i3,i4)))
     &  )/(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))))/
     &  (2._dp*(s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))
        fac=fac+
     &  (I3m(s(i2,i3),s(i1,i4),s(i5,i6))*za(i2,i3)*
     &  ((-3._dp*za(i1,i4)*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  zb(i2,i3)*((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i3)*
     &  zb(i1,i4)-
     &  (-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i2,i3)*zb(i2,i4))*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)**2-
     &  (za(i4,i5)*zb(i1,i2)*
     &  (za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6)))/t(i1,i2,i3)-
     &  (-(za(i3,i4)*zb(i1,i6)*zb(i2,i3)*
     &  ((-s(i1,i4)+s(i2,i3)-s(i5,i6))*za(i1,i5)-
     &  za(i4,i5)*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))+
     &  2._dp*za(i1,i4)*za(i5,i6)*zb(i4,i6)))+
     &  (za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))*
     &  (-(za(i1,i4)*zb(i2,i4)*
     &  (-(za(i2,i5)*zb(i2,i6))-za(i3,i5)*zb(i3,i6)))+
     &  3._dp*za(i1,i3)*za(i4,i5)*zb(i2,i3)*zb(i4,i6))-
     &  za(i4,i5)*zb(i1,i2)*
     &  ((s(i1,i4)-s(i2,i3)-s(i5,i6))*za(i1,i4)*zb(i4,i6)-
     &  (-s(i1,i4)-s(i2,i3)+s(i5,i6))*za(i1,i5)*zb(i5,i6))-
     &  (za(i1,i3)*za(i1,i5)*zb(i2,i3)*
     &  (za(i2,i4)*zb(i1,i2)*zb(i4,i6)+
     &  za(i3,i4)*zb(i1,i3)*zb(i4,i6))*
     &  (-t(i1,i2,i3)+t(i2,i3,i4)))/
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))/
     &  (s(i1,i4)**2-2._dp*s(i1,i4)*s(i2,i3)+s(i2,i3)**2-
     &  2._dp*s(i1,i4)*s(i5,i6)-2._dp*s(i2,i3)*s(i5,i6)+s(i5,i6)**2)))/
     &  (za(i1,i2)*(za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)))

      zaj_virt_a6_lc_fsc_pm = Fsc4+za(i1,i3)/za(i2,i3)*fac

      end function

      function zaj_virt_a6_lc_fsc_pp(i1,i3,i2,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a6_lc_fsc_pp

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        complex(dp)  :: L0,L1,Lsm1_2me,Lnrat,Ls1

        zaj_virt_a6_lc_fsc_pp = 
     &  (Lsm1_2me(t(i1,i2,i3),t(i2,i3,i4),s(i2,i3),s(i5,i6))*za(i1,i3)*za(
     &  i1,i5)**2*
     &  za(i3,i4))/(za(i1,i2)*za(i1,i4)**3*za(i2,i3)*za(i5,i6))-
     &  (za(i1,i3)*za(i3,i5)**2)/
     &  (2._dp*za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))+
     &  (Lnrat(-t(i2,i3,i4),-s(i5,i6))*za(i1,i3)*za(i3,i5)**2)/
     &  (2._dp*za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))-
     &  (L0(-t(i1,i2,i3),-s(i2,i3))*za(i1,i3)*za(i1,i5)**2*zb(i1,i2))/
     &  (s(i2,i3)*za(i1,i2)*za(i1,i4)**2*za(i5,i6))+
     &  (L0(-t(i2,i3,i4),-s(i5,i6))*za(i1,i3)*za(i1,i5)**2*
     &  (za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4)))/
     &  (s(i5,i6)*za(i1,i2)*za(i1,i4)**2*za(i2,i3)*za(i5,i6))+
     &  (L0(-t(i2,i3,i4),-s(i5,i6))*za(i1,i3)**2*za(i3,i5)*zb(i1,i6))/
     &  (s(i5,i6)*za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4))+
     &  (L1(-t(i2,i3,i4),-s(i5,i6))*za(i1,i3)**3*za(i5,i6)*zb(i1,i6)**2)/
     &  (2._dp*s(i5,i6)**2*za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4))+
     &  (L0(-t(i2,i3,i4),-s(i2,i3))*za(i1,i5)**2*za(i3,i4)*zb(i2,i4))/
     &  (s(i2,i3)*za(i1,i2)*za(i1,i4)**2*za(i5,i6))-
     &  (L0(-t(i1,i2,i3),-s(i5,i6))*za(i1,i3)*za(i1,i5)*za(i4,i5)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))/
     &  (s(i5,i6)*za(i1,i2)*za(i1,i4)**2*za(i2,i3)*za(i5,i6))-
     &  (L0(-t(i1,i2,i3),-s(i5,i6))*za(i1,i3)*za(i3,i5)*zb(i4,i6))/
     &  (s(i5,i6)*za(i1,i2)*za(i1,i4)*za(i2,i3))+
     &  (L1(-t(i1,i2,i3),-s(i5,i6))*za(i1,i3)*za(i3,i4)*za(i5,i6)*zb(i4,i6
     &  )**2)/
     &  (s(i5,i6)**2*za(i1,i2)*za(i1,i4)*za(i2,i3))+
     &  (L1(-s(i3,i4),-t(i2,i3,i4))*za(i2,i5)**2*za(i3,i4)*zb(i2,i4)**2)/
     &  (2._dp*za(i1,i2)*za(i2,i4)*za(i5,i6)*t(i2,i3,i4)**2)+
     &  (Ls1(-s(i2,i3),-t(i2,i3,i4),-s(i3,i4),-t(i2,i3,i4))*za(i2,i5)**2*z
     &  a(i3,i4)*
     &  zb(i2,i4)**2)/(za(i1,i2)*za(i2,i4)*za(i5,i6)*t(i2,i3,i4)**2)-
     &  (L1(-s(i2,i3),-t(i2,i3,i4))*za(i2,i3)*za(i4,i5)**2*zb(i2,i4)**2)/
     &  (2._dp*za(i1,i4)*za(i2,i4)*za(i5,i6)*t(i2,i3,i4)**2)-
     &  (za(i1,i5)*za(i2,i4)*za(i3,i5)*zb(i2,i4)**2)/
     &  (2._dp*za(i1,i2)*za(i1,i4)*za(i5,i6)*t(i1,i2,i4)*t(i2,i3,i4))-
     &  (zb(i2,i4)*(za(i1,i2)*zb(i2,i6)+za(i1,i4)*zb(i4,i6))*
     &  (-(za(i2,i3)*zb(i2,i6))+za(i3,i4)*zb(i4,i6)))/
     &  (2._dp*za(i1,i2)*za(i1,i4)*zb(i5,i6)*t(i1,i2,i4)*t(i2,i3,i4))
 
      end function

      function zaj_virt_a64v_fvf_pp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a64v_fvf_pp
        complex(dp) :: s,t,Lsm1_2me

        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_virt_a64v_fvf_pp =
     &    -za(i2,i5)**2/(za(i1,i2)*za(i5,i6)*za(i3,i4)**2)
     &    *Lsm1_2me(t(i1,i2,i3),t(i1,i2,i4),s(i1,i2),s(i5,i6))

      end function

      function zaj_virt_a64v_fvf_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a64v_fvf_mp
       
        zaj_virt_a64v_fvf_mp = zaj_virt_a64v_fvf_pm(i1,i2,i4,i3,i5,i6,za,zb)
    
      end function

      function zaj_virt_a64v_fvf_pm(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a64v_fvf_pm
        complex(dp) :: s,t
        complex(dp):: Lsm1_2mh,I3m,zab2,I3m123456,I3m563412

        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        I3m123456=I3m(s(i1,i2),s(i3,i4),s(i5,i6))
        I3m563412=I3m123456


        zaj_virt_a64v_fvf_pm = 
     &   -zab2(i5,i2,i4,i1)**2/(za(i5,i6)*zb(i1,i2)*zab2(i3,i1,i2,i4)**2)
     &   *Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),s(i5,i6))
     &   -zab2(i2,i1,i3,i6)**2/(za(i1,i2)*zb(i5,i6)*zab2(i3,i1,i2,i4)**2)
     &   *Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))

     &   -I3m123456*za(i4,i5)*zb(i1,i3)
     &   *zab2(i2,i1,i3,i6)/(2._dp*zab2(i3,i1,i2,i4)*t(i1,i2,i3))

     &   -I3m123456*za(i2,i4)*zb(i3,i6)
     &   *zab2(i5,i2,i4,i1)/(2._dp*zab2(i3,i1,i2,i4)*t(i1,i2,i4))

     &   +I3m563412*za(i2,i4)*zb(i3,i6)
     &   *zab2(i5,i3,i6,i1)/(2._dp*zab2(i3,i1,i2,i4)*t(i3,i5,i6))

     &   +I3m563412*za(i4,i5)*zb(i1,i3)
     &   *zab2(i2,i4,i5,i6)/(2._dp*zab2(i3,i1,i2,i4)*t(i4,i5,i6))


      end function

      function zaj_virt_a64v_fvs_pp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a64v_fvs_pp
        complex(dp) :: L0,L1,Lsm1_2me,Lnrat

        complex(dp) :: s,t, Brackpp, Brackppa
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

C---This is the square bracket in Eq.(11.1) subject to exchange 16,25
        Brackppa(i1,i2,i3,i4,i5,i6)=
     &   za(i1,i2)*(za(i5,i3)*zb(i3,i1))**2/za(i5,i6)/za(i3,i4)**2
     &   *L1(-t(i1,i2,i3),-s(i1,i2))/s(i1,i2)**2
     &   +za(i5,i2)*za(i5,i3)*zb(i3,i1)/za(i5,i6)/za(i3,i4)**2
     &   *L0(-t(i1,i2,i3),-s(i1,i2))/s(i1,i2)

C---This is the whole expression in Eq.(11.1) subject to exchange 34
        Brackpp(i1,i2,i3,i4,i5,i6)=
     &   -za(i2,i3)*za(i2,i4)*za(i3,i5)*za(i4,i5)
     &   /za(i1,i2)/za(i5,i6)/za(i3,i4)**4
     &   *Lsm1_2me(t(i1,i2,i3),t(i1,i2,i4),s(i1,i2),s(i5,i6))
     &    -(za(i2,i4)*za(i3,i5)+za(i2,i3)*za(i4,i5))
     &   /za(i1,i2)/za(i5,i6)/za(i3,i4)**3
     &   *(za(i5,i3)*zb(i3,i1)*za(i1,i2)
     &   *L0(-t(i1,i2,i3),-s(i1,i2))/s(i1,i2)
     &   +za(i5,i6)*zb(i6,i4)*za(i4,i2)
     &   *L0(-t(i1,i2,i3),-s(i5,i6))/s(i5,i6)
     &   +0.5_dp*za(i5,i2)*(Lnrat(-t(i1,i2,i3),-t(i1,i2,i4))))
     &   -Brackppa(i1,i2,i3,i4,i5,i6)-Brackppa(i6,i5,i3,i4,i2,i1)
     &   +0.5_dp/za(i3,i4)**2
     &   *(za(i2,i5)**2/za(i1,i2)/za(i5,i6)
     &   -zb(i1,i6)**2/zb(i1,i2)/zb(i5,i6))

c---exch34:(3<-->4)
       zaj_virt_a64v_fvs_pp =  Brackpp(i1,i2,i3,i4,i5,i6)
     &                       + Brackpp(i1,i2,i4,i3,i5,i6)

      end function

      function zaj_virt_a64v_fvs_mp(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a64v_fvs_mp

        zaj_virt_a64v_fvs_mp = zaj_virt_a64v_fvs_pm(i1,i2,i4,i3,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvs_pm(i1,i2,i3,i4,i5,i6,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: zaj_virt_a64v_fvs_pm

        zaj_virt_a64v_fvs_pm =  Brackpm(i1,i2,i3,i4,i5,i6,za,zb)
     &                        + Brackpm(i2,i1,i4,i3,i6,i5,zb,za)
 
      contains
        function Brackpm(i1,i2,i3,i4,i5,i6,za,zb) 
          implicit none
          complex(dp):: Brackpm
C---This is the whole expression in Eq.(11.6) subject to exchange flip_2
c---ie  flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)
      
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex(dp):: Lsm1_2mh,I3m,Lnrat,zab2
          real(dp):: s,t,delta,IDelta

C---statement functions
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
          delta(i1,i2,i3,i4,i5,i6)=s(i1,i2)-s(i3,i4)-s(i5,i6)
          zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

          IDelta=1._dp/(s(i1,i2)**2+s(i3,i4)**2+s(i5,i6)**2
     &     -2._dp*(+s(i1,i2)*s(i3,i4)+s(i1,i2)*s(i5,i6)+s(i5,i6)*s(i3,i4)))

          Brackpm=
     &     -2._dp*za(i2,i3)*zb(i4,i6)*zab2(i3,i1,i2,i6)
     &     *zab2(i2,i5,i6,i4)*t(i1,i2,i3)
     &     /za(i1,i2)/zb(i5,i6)/zab2(i3,i1,i2,i4)**4
     &     *Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))

     &    +(2._dp*za(i2,i3)*zb(i4,i6)*(zb(i1,i2)*za(i5,i6)*za(i2,i3)*zb(i4,i6)

     &     -zab2(i5,i2,i4,i1)*zab2(i3,i1,i2,i4))*zab2(i4,i1,i2,i3)
     &     *(t(i1,i2,i3)-t(i1,i2,i4))/zab2(i3,i1,i2,i4)**3*IDelta

     &     -3._dp*(s(i3,i4)*delta(i3,i4,i1,i2,i5,i6)*zab2(i2,i3,i4,i1)
     &      *zab2(i5,i3,i4,i6)*IDelta
     &     -za(i2,i3)*zb(i3,i1)*za(i5,i4)*zb(i4,i6)
     &     -za(i2,i4)*zb(i4,i1)*za(i5,i3)*zb(i3,i6))

     &     *zab2(i4,i1,i2,i3)/zab2(i3,i1,i2,i4)*IDelta
     &     -zb(i1,i4)*za(i3,i5)*za(i2,i3)*zb(i4,i6)
     &     *zab2(i4,i1,i2,i3)**2/zab2(i3,i1,i2,i4)**2*IDelta
     &     -zb(i1,i3)*za(i4,i5)*za(i2,i4)*zb(i3,i6)*IDelta)
     &     *I3m(s(i1,i2),s(i3,i4),s(i5,i6))

     &     +(2._dp*za(i2,i3)*zb(i4,i6)*t(i1,i2,i3)*zab2(i2,i1,i4,i6)
     &     /za(i1,i2)/zb(i5,i6)/zab2(i3,i1,i2,i4)**3

     &     -(zab2(i2,i1,i4,i6)**2
     &      +2._dp*za(i2,i3)*zb(i3,i6)*za(i2,i4)*zb(i4,i6))
     &     /za(i1,i2)/zb(i5,i6)/zab2(i3,i1,i2,i4)**2)
     &     *Lnrat(-t(i1,i2,i3),-s(i3,i4))

c---ie      flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
     &      +Brackpma(i1,i2,i3,i4,i5,i6,za,zb)
     &      +Brackpma(i5,i6,i4,i3,i1,i2,zb,za)

     &      +zb(i1,i6)*zab2(i4,i1,i2,i3)
     &      *(zb(i1,i6)*delta(i3,i4,i5,i6,i1,i2)
     &     -2._dp*zb(i1,i2)*za(i2,i5)*zb(i5,i6))
     &     /zb(i1,i2)/zb(i5,i6)/zab2(i3,i1,i2,i4)*IDelta

     &      +(zab2(i2,i1,i4,i6)**2+za(i2,i1)*zb(i1,i6)*za(i2,i5)*zb(i5,i6))
     &     /za(i1,i2)/zb(i5,i6)/zab2(i3,i1,i2,i4)**2

        return
        end
       
        function Brackpma(i1,i2,i3,i4,i5,i6,za,zb) 
          complex(dp):: Brackpma
C---This is the curly braces expression in Eq.(11.6) subject to flip_3
c---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      
          complex(dp) :: L0,L1,Lnrat,zab2
          real(dp) :: s,t,delta,IDelta

C---statement functions
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
          delta(i1,i2,i3,i4,i5,i6)=s(i1,i2)-s(i3,i4)-s(i5,i6)
          zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

          IDelta=1._dp/(s(i1,i2)**2+s(i3,i4)**2+s(i5,i6)**2
     &     -2._dp*(s(i1,i2)*s(i3,i4)+s(i1,i2)*s(i5,i6)+s(i5,i6)*s(i3,i4)))

           Brackpma= 
     &     2._dp*zb(i1,i2)*za(i2,i3)*zb(i3,i6)
     &     /zb(i5,i6)/zab2(i3,i1,i2,i4)**2
     &     *(zab2(i3,i1,i2,i6)*zab2(i2,i1,i3,i4)/zab2(i3,i1,i2,i4)
     &     *L0(-t(i1,i2,i3),-s(i1,i2))/s(i1,i2)
     &     -za(i2,i3)*zb(i3,i6)
     &     *(L0(-t(i1,i2,i3),-s(i1,i2))/s(i1,i2)
     &     -1/2._dp*L1(-t(i1,i2,i3),-s(i1,i2))/s(i1,i2)))

     &     -(3._dp*delta(i5,i6,i3,i4,i1,i2)*zab2(i2,i3,i4,i1)
     &     *zab2(i5,i3,i4,i6)*zab2(i4,i1,i2,i3)
     &     /zab2(i3,i1,i2,i4)*IDelta**2

     &     +2._dp*za(i3,i2)*zb(i2,i1)*zb(i4,i6)*(t(i1,i2,i3)-t(i1,i2,i4))
     &     *(za(i2,i5)*t(i1,i2,i3)+za(i2,i1)*zb(i1,i6)*za(i6,i5))
     &     /zab2(i3,i1,i2,i4)**3*IDelta

     &     -(2._dp*za(i2,i3)*zb(i3,i6)*(t(i1,i2,i3)-t(i1,i2,i4))
     &     +delta(i1,i2,i3,i4,i5,i6)
     &     *(za(i2,i5)*t(i1,i2,i3)
     &     +za(i2,i1)*zb(i1,i6)*za(i6,i5))/za(i5,i6))
     &     *za(i5,i2)*zb(i2,i1)/zab2(i3,i1,i2,i4)**2*IDelta
 
     &     +0.5_dp*(za(i1,i2)*zb(i1,i6)**2/zb(i5,i6)
     &      +zb(i1,i2)*za(i2,i5)**2/za(i5,i6)
     &     -2._dp*za(i2,i5)*zb(i1,i6))*zab2(i4,i1,i2,i3)
     &     /zab2(i3,i1,i2,i4)*IDelta)
     &     *(Lnrat(-s(i1,i2),-s(i3,i4)))
     
        return
        end

      end function

      function zaj_virt_a64ax_pp(i1,i2,i3,i4,i5,i6,za,zb)
        include 'masses.f'
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a64ax_pp
        complex(dp) :: L0,L1,Lsm1_2me,Lnrat
        real(dp) :: s,t,mtsq  

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
        
        mtsq=mt**2

        zaj_virt_a64ax_pp =
     &    (Lnrat(-t(i1,i2,i3),-s(i5,i6))*za(i2,i5)**2)/(za(i1,i2)*za(i3,i4)*
     &    *2*za(i5,i6))-
     &    (Lnrat(-t(i1,i2,i4),-s(i5,i6))*za(i2,i5)**2)/(za(i1,i2)*za(i3,i4)*
     &    *2*za(i5,i6))-
     &    (Lsm1_2me(t(i1,i2,i3),t(i1,i2,i4),s(i1,i2),s(i5,i6))*za(i2,i5)*
     &    (za(i2,i4)*za(i3,i5)+za(i2,i3)*za(i4,i5)))/
     &    (2._dp*za(i1,i2)*za(i3,i4)**3*za(i5,i6))-
     &    (L0(-t(i1,i2,i3),-s(i1,i2))*za(i2,i5)*za(i3,i5)*zb(i1,i3))/
     &    (s(i1,i2)*za(i3,i4)**2*za(i5,i6))+
     &    (L0(-t(i1,i2,i4),-s(i1,i2))*za(i2,i5)*za(i4,i5)*zb(i1,i4))/
     &    (s(i1,i2)*za(i3,i4)**2*za(i5,i6))+
     &    (((L1(-t(i1,i2,i4),-s(i5,i6))*s(i3,i4))/s(i5,i6)**2+
     &    L0(-t(i1,i2,i4),-s(i5,i6))/s(i5,i6))*za(i2,i3)*za(i2,i5)*zb(i3,i6)
     &    )/
     &    (za(i1,i2)*za(i3,i4)**2)-(L1(-t(i1,i2,i4),-s(i5,i6))*za(i2,i3)*za(
     &    i2,i5)*
     &    zb(i1,i3)*zb(i3,i6))/(s(i5,i6)**2*za(i2,i4)*za(i3,i4))-
     &    (((L1(-t(i1,i2,i3),-s(i5,i6))*s(i3,i4))/s(i5,i6)**2+
     &    L0(-t(i1,i2,i3),-s(i5,i6))/s(i5,i6))*za(i2,i4)*za(i2,i5)*zb(i4,i6)
     &    )/
     &    (za(i1,i2)*za(i3,i4)**2)+(L1(-t(i1,i2,i3),-s(i5,i6))*(s(i1,i4)+s(i
     &    3,i4))*
     &    za(i2,i5)*zb(i4,i6))/(s(i5,i6)**2*za(i1,i3)*za(i3,i4))-
     &    (za(i2,i5)*(((-(za(i2,i3)*zb(i1,i3))-za(i2,i4)*zb(i1,i4))*zb(i3,i6
     &    ))/
     &    za(i2,i4)+(zb(i4,i6)*t(i1,i3,i4))/za(i1,i3)))/
     &    (12._dp*mtsq*s(i5,i6)*za(i3,i4))

      end function

      function zaj_virt_a64ax_pm(i1,i2,i3,i4,i5,i6,za,zb)
        include 'masses.f'
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a64ax_pm
        complex(dp) :: L0,L1,Lnrat,Lsm1_2mh,i3m
        real(dp) :: s,t,mtsq  

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

        mtsq=mt**2

        zaj_virt_a64ax_pm= 
     &  (Lnrat(-s(i5,i6),-s(i3,i4))*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i
     &  4))**2)/
     &  (za(i5,i6)*zb(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &  *2)+
     &  (L0(-t(i1,i2,i4),-s(i1,i2))*za(i2,i4)*
     &  (za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &  (-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4)))/
     &  (s(i1,i2)*za(i5,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**
     &  2)-
     &  (L0(-t(i1,i2,i4),-s(i5,i6))*zb(i1,i3)*
     &  (za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &  (-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))*
     &  (-(za(i1,i3)*zb(i1,i2))-za(i3,i4)*zb(i2,i4)))/
     &  (s(i5,i6)*za(i5,i6)*zb(i1,i2)*zb(i2,i4)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2)+
     &  (zb(i1,i3)*(-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4))*zb(i3,i6))
     &  /
     &  (12._dp*mtsq*s(i5,i6)*zb(i2,i4)*zb(i3,i4))-
     &  (za(i2,i4)*za(i4,i5)*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6)))
     &  /
     &  (12._dp*mtsq*s(i5,i6)*za(i1,i3)*za(i3,i4))+
     &  (za(i2,i4)*za(i3,i5)*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6)))
     &  /
     &  (s(i5,i6)*za(i1,i3)*za(i3,i4)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     &  (zb(i1,i3)*(-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4))*zb(i4,i6))
     &  /
     &  (s(i5,i6)*zb(i2,i4)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &  zb(i3,i4))+(za(i4,i5)*zb(i1,i3)*
     &  (-(((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i4)*za(i3,i5))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &  2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &  za(i3,i4)*za(i5,i6)))+
     &  ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i3,i5)*zb(i1,i3))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &  2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)
     &  *
     &  zb(i1,i2))-(2._dp*za(i2,i4)*zb(i4,i6))/
     &  (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &  6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     &  (zb(i1,i3)*zb(i4,i6))/(s(i5,i6)*zb(i1,i2)*zb(i3,i4))-
     &  ((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i3)*zb(i4,i6))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i1,i2)*zb(i3,i4))))/
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))+
     &  (za(i2,i4)*zb(i3,i6)*((za(i2,i4)*za(i3,i5))/
     &  (s(i5,i6)*za(i1,i2)*za(i3,i4))-
     &  ((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i4)*za(i3,i5))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &  2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  za(i3,i4))-(2._dp*za(i3,i5)*zb(i1,i3))/
     &  (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &  6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     &  ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i4)*zb(i4,i6))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &  2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &  *
     &  zb(i5,i6))-((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i3)*zb(i4,i6))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i3,i4)*zb(i5,i6))))/
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))-
     &  (L0(-t(i1,i2,i3),-s(i1,i2))*zb(i1,i3)*
     &  (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))/
     &  (s(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2*zb(i5,i6
     &  ))-
     &  (L0(-t(i1,i2,i3),-s(i5,i6))*za(i2,i4)*
     &  (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))*
     &  (za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4))*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))/
     &  (s(i5,i6)*za(i1,i2)*za(i1,i3)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2*zb(i5,i6))+
     &  (Lnrat(-s(i5,i6),-s(i3,i4))*(-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i
     &  3,i6))**2)/
     &  (za(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2*zb(i5,i
     &  6))-
     &  (L1(-s(i5,i6),-t(i1,i2,i3))*za(i1,i4)*za(i2,i4)*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6))*zb(i4,i6))/
     &  (za(i1,i2)*za(i1,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &  zb(i5,i6)*t(i1,i2,i3))-(Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),
     &  s(i5,i6))*((-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6))**2*
     &  (-(za(i1,i2)*zb(i1,i4))+za(i2,i3)*zb(i3,i4))**2-
     &  za(i2,i3)**2*zb(i4,i6)**2*t(i1,i2,i3)**2))/
     &  (2._dp*za(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**4*zb(
     &  i5,i6))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  +Lnrat(-s(i1,i2),-s(i3,i4))*((-6._dp*zb(i1,i2)*
     &  ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i5)-
     &  2._dp*za(i1,i2)*za(i5,i6)*zb(i1,i6))*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     &  (za(i2,i4)*zb(i1,i3)*zb(i3,i6)**2)/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i3,i4)*zb(i5,i6))-
     &  (zb(i1,i3)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6))*zb(i4,i6))
     &  /
     &  ((-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2*zb(i3,i4)*
     &  zb(i5,i6))+(zb(i1,i4)*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6))*
     &  (3._dp*(-(za(i1,i3)*zb(i1,i4)*zb(i3,i6))-
     &  za(i2,i3)*zb(i2,i4)*zb(i3,i6))-
     &  zb(i4,i6)*(t(i1,i2,i3)-t(i1,i2,i4))))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2*zb(i3,i4)*zb(i5,i6
     &  )
     &  ))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  +I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &  (-((za(i2,i4)*za(i4,i5)*zb(i1,i3)*zb(i3,i6))/
     &  (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &  6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2))-
     &  (3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     &  ((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*zb(i1,i2)+
     &  (-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     &  (3._dp*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-(za(i1,i2)*za(i2,i5)*zb(i1,i2)*zb(i1,i6))-
     &  za(i2,i3)*za(i4,i5)*zb(i1,i4)*zb(i3,i6)-
     &  za(i2,i4)*za(i3,i5)*zb(i1,i3)*zb(i4,i6)-
     &  za(i2,i5)*za(i5,i6)*zb(i1,i6)*zb(i5,i6)))/
     &  (2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     &  (i5,i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     &  (za(i4,i5)*zb(i1,i3)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6)))
     &  /
     &  (2._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*t(i1,i2,i3))+
     &  (za(i3,i5)*zb(i1,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-(za(i1,i2)*zb(i1,i6))+za(i2,i3)*zb(i3,i6))*
     &  (t(i1,i2,i3)-t(i1,i2,i4)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  -(L1(-s(i5,i6),-t(i1,i2,i4))*za(i3,i5)*zb(i1,i3)*
     &  (za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*zb(i2,i3))/
     &  (za(i5,i6)*zb(i1,i2)*zb(i2,i4)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*t(i1,i2,i4))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  -(Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),s(i5,i6))*
     &  ((za(i2,i3)*zb(i1,i2)-za(i3,i4)*zb(i1,i4))**2*
     &  (-(za(i1,i5)*zb(i1,i4))-za(i2,i5)*zb(i2,i4))**2-
     &  za(i3,i5)**2*zb(i1,i4)**2*t(i1,i2,i4)**2))/
     &  (2._dp*za(i5,i6)*zb(i1,i2)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i
     &  4))**4)


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  +I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &  (-((za(i2,i4)*za(i4,i5)*zb(i1,i3)*zb(i3,i6))/
     &  (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &  6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2))-
     &  (3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     &  (za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6))-
     &  (-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     &  (3._dp*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (-(za(i1,i2)*za(i2,i5)*zb(i1,i2)*zb(i1,i6))-
     &  za(i2,i3)*za(i4,i5)*zb(i1,i4)*zb(i3,i6)-
     &  za(i2,i4)*za(i3,i5)*zb(i1,i3)*zb(i4,i6)-
     &  za(i2,i5)*za(i5,i6)*zb(i1,i6)*zb(i5,i6)))/
     &  (2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     &  (i5,i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))-
     &  (za(i2,i4)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*zb(i3,i6))/
     &  (2._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*t(i1,i2,i4))+
     &  (za(i2,i3)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*zb(i4,i6)*
     &  (-t(i1,i2,i3)+t(i1,i2,i4)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  +Lnrat(-s(i1,i2),-s(i3,i4))*(-((za(i2,i4)*za(i4,i5)**2*zb(i1,i3))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i3,i4)*za(i5,i6)))-
     &  (za(i2,i4)*za(i3,i5)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4)))/
     &  (za(i3,i4)*za(i5,i6)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &  *
     &  2)+(6._dp*za(i1,i2)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*zb(i1,i6)-
     &  2._dp*za(i2,i5)*zb(i1,i2)*zb(i5,i6)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4)))+
     &  (za(i2,i3)*(za(i2,i5)*zb(i1,i2)+za(i4,i5)*zb(i1,i4))*
     &  (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &  (3._dp*(-(za(i1,i3)*za(i4,i5)*zb(i1,i4))-
     &  za(i2,i3)*za(i4,i5)*zb(i2,i4))-
     &  za(i3,i5)*(-t(i1,i2,i3)+t(i1,i2,i4))))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i3,i4)*za(i5,i6)*
     &  (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))**2))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  +Lnrat(-s(i5,i6),-s(i3,i4))*((za(i4,i5)*zb(i1,i3)**2*zb(i3,i6))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i1,i2)*zb(i3,i4))+
     &  (zb(i1,i4)*(za(i3,i5)*zb(i1,i3)-za(i5,i6)*zb(i1,i6))*zb(i3,i6))/
     &  (zb(i1,i2)*zb(i3,i4)*(-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))*
     &  *
     &  2)+(6._dp*(za(i3,i5)*zb(i1,i3)-za(i5,i6)*zb(i1,i6))*
     &  (-((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i5))+
     &  2._dp*za(i1,i2)*za(i5,i6)*zb(i1,i6))*
     &  (-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))*zb(i5,i6))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &  (-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6)))+
     &  ((za(i3,i5)*zb(i1,i3)-za(i5,i6)*zb(i1,i6))*
     &  (-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))*zb(i4,i6)*
     &  (3._dp*(za(i3,i5)*zb(i1,i3)*zb(i4,i5)+
     &  za(i3,i6)*zb(i1,i3)*zb(i4,i6))+
     &  zb(i1,i4)*(t(i3,i5,i6)-t(i4,i5,i6))))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i1,i2)*zb(i3,i4)*
     &  (-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))**2))


        zaj_virt_a64ax_pm=zaj_virt_a64ax_pm
     &  +Lnrat(-s(i5,i6),-s(i3,i4))*((za(i2,i4)**2*za(i4,i5)*zb(i3,i6))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)*za(i3,i4))+
     &  (za(i2,i3)*za(i4,i5)*(za(i2,i4)*zb(i4,i6)+za(i2,i5)*zb(i5,i6)))/
     &  (za(i1,i2)*za(i3,i4)*(-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))*
     &  *2)
     &  -(6._dp*za(i5,i6)*(-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))*
     &  (za(i2,i4)*zb(i4,i6)+za(i2,i5)*zb(i5,i6))*
     &  (-((-s(i1,i2)+s(i3,i4)-s(i5,i6))*zb(i1,i6))+
     &  2._dp*za(i2,i5)*zb(i1,i2)*zb(i5,i6)))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &  (-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6)))+
     &  (za(i3,i5)*(-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))*
     &  (za(i2,i4)*zb(i4,i6)+za(i2,i5)*zb(i5,i6))*
     &  (3._dp*(za(i2,i4)*za(i3,i5)*zb(i4,i5)+
     &  za(i2,i4)*za(i3,i6)*zb(i4,i6))+
     &  za(i2,i3)*(-t(i3,i5,i6)+t(i4,i5,i6))))/
     &  ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &  i6)-
     &  2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)*za(i3,i4)*
     &  (-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))**2))

      end function

      function zaj_virt_a64ax_mp(i1,i2,i3,i4,i5,i6,za,zb)
        include 'masses.f'
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_a64ax_mp
        complex(dp) :: L0,L1,Lnrat,Lsm1_2mh,i3m
        real(dp) :: s,t,mtsq  

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

        mtsq=mt**2

        zaj_virt_a64ax_mp=
     &    (za(i2,i3)**2*za(i3,i5)*zb(i1,i6))/(12._dp*mtsq*s(i5,i6)*za(i2,i4)*z
     &    a(i3,i4))-
     &    (Lnrat(-s(i5,i6),-s(i3,i4))*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i
     &    3))**2)/
     &    (za(i5,i6)*zb(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    *2)-
     &    (za(i2,i3)**2*za(i4,i5)*zb(i1,i6))/
     &    (s(i5,i6)*za(i2,i4)*za(i3,i4)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (L0(-t(i1,i2,i3),-s(i1,i2))*za(i2,i3)*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &    (-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)))/
     &    (s(i1,i2)*za(i5,i6)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**
     &    2)+
     &    (L0(-t(i1,i2,i3),-s(i5,i6))*(za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i
     &    3))*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*zb(i1,i4)*
     &    (-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3)))/
     &    (s(i5,i6)*za(i5,i6)*zb(i1,i2)*zb(i1,i3)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2)+
     &    (za(i2,i5)*zb(i1,i4)**2*zb(i3,i6))/
     &    (s(i5,i6)*zb(i1,i3)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    zb(i3,i4))-(za(i3,i5)*zb(i1,i4)*
     &    (((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i3)*za(i4,i5))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &    2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i3,i4)
     &    *
     &    za(i5,i6))+((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i4,i5)*
     &    zb(i1,i4))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &    2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i5,i6)
     &    *
     &    zb(i1,i2))-(2._dp*za(i2,i3)*zb(i3,i6))/
     &    (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &    6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)-
     &    (zb(i1,i4)*zb(i3,i6))/(s(i5,i6)*zb(i1,i2)*zb(i3,i4))+
     &    ((-s(i1,i2)-s(i3,i4)+s(i5,i6))*zb(i1,i4)*zb(i3,i6))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i1,i2)*zb(i3,i4))))/
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))-
     &    (za(i2,i5)*zb(i1,i4)**2*zb(i4,i6))/(12._dp*mtsq*s(i5,i6)*zb(i1,i3)*z
     &    b(i3,i4))-
     &    (za(i2,i3)*zb(i4,i6)*(-((za(i2,i3)*za(i4,i5))/
     &    (s(i5,i6)*za(i1,i2)*za(i3,i4)))+
     &    ((-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i3)*za(i4,i5))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &    2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &    *
     &    za(i3,i4))-(2._dp*za(i4,i5)*zb(i1,i4))/
     &    (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &    6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)+
     &    ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i3)*zb(i3,i6))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-
     &    2._dp*s(i1,i2)*s(i5,i6)-2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)
     &    *
     &    zb(i5,i6))+((s(i1,i2)-s(i3,i4)-s(i5,i6))*zb(i1,i4)*zb(i3,i6))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i3,i4)*zb(i5,i6))))/
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))+
     &    (L0(-t(i1,i2,i4),-s(i1,i2))*zb(i1,i4)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    (s(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2*zb(i5,i6
     &    ))+
     &    (L0(-t(i1,i2,i4),-s(i5,i6))*za(i2,i3)*
     &    (-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    (s(i5,i6)*za(i1,i2)*za(i2,i4)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2*zb(i5,i6))-
     &    (Lnrat(-s(i5,i6),-s(i3,i4))*(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i
     &    4,i6))**2)/
     &    (za(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2*zb(i5,i
     &    6))+
     &    (L1(-s(i5,i6),-t(i1,i2,i3))*za(i4,i5)*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*zb(i1,i4)**2)/
     &    (za(i5,i6)*zb(i1,i2)*zb(i1,i3)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*t(i1,i2,i3))+
     &    (Lsm1_2mh(s(i3,i4),t(i1,i2,i3),s(i1,i2),s(i5,i6))*
     &    ((za(i2,i4)*zb(i1,i2)+za(i3,i4)*zb(i1,i3))**2*
     &    (-(za(i1,i5)*zb(i1,i3))-za(i2,i5)*zb(i2,i3))**2-
     &    za(i4,i5)**2*zb(i1,i3)**2*t(i1,i2,i3)**2))/
     &    (2._dp*za(i5,i6)*zb(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i
     &    3))**4)


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    -Lnrat(-s(i1,i2),-s(i3,i4))*((za(i2,i3)*za(i3,i5)**2*zb(i1,i4))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i3,i4)*za(i5,i6))+
     &    (za(i2,i3)*za(i4,i5)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3)))/
     &    (za(i3,i4)*za(i5,i6)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    *
     &    2)+(6._dp*za(i1,i2)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*zb(i1,i6)-
     &    2._dp*za(i2,i5)*zb(i1,i2)*zb(i5,i6)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (za(i2,i4)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (3._dp*(-(za(i1,i4)*za(i3,i5)*zb(i1,i3))-
     &    za(i2,i4)*za(i3,i5)*zb(i2,i3))-
     &    za(i4,i5)*(t(i1,i2,i3)-t(i1,i2,i4))))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i3,i4)*za(i5,i6)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2))


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    -I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    (-((za(i2,i3)*za(i3,i5)*zb(i1,i4)*zb(i4,i6))/
     &    (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &    6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2))-
     &    (3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     &    (za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i1,i2)*zb(i1,i6))-
     &    (-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i2,i5)*zb(i5,i6)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (3._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*za(i2,i5)*zb(i1,i2)*zb(i1,i6))-
     &    za(i2,i3)*za(i4,i5)*zb(i1,i4)*zb(i3,i6)-
     &    za(i2,i4)*za(i3,i5)*zb(i1,i3)*zb(i4,i6)-
     &    za(i2,i5)*za(i5,i6)*zb(i1,i6)*zb(i5,i6)))/
     &    (2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     &    (i5,i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (za(i2,i3)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*zb(i4,i6))/
     &    (2._dp*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*t(i1,i2,i3))+
     &    (za(i2,i4)*(za(i2,i5)*zb(i1,i2)+za(i3,i5)*zb(i1,i3))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*zb(i3,i6)*
     &    (t(i1,i2,i3)-t(i1,i2,i4)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2))


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    +(L1(-s(i5,i6),-t(i1,i2,i4))*za(i2,i3)**2*zb(i3,i6)*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    (za(i1,i2)*za(i2,i4)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*
     &    zb(i5,i6)*t(i1,i2,i4))


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    +(Lsm1_2mh(s(i3,i4),t(i1,i2,i4),s(i1,i2),
     &    s(i5,i6))*((-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6))**2*
     &    (-(za(i1,i2)*zb(i1,i3))-za(i2,i4)*zb(i3,i4))**2-
     &    za(i2,i4)**2*zb(i3,i6)**2*t(i1,i2,i4)**2))/
     &    (2._dp*za(i1,i2)*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**4*zb(
     &    i5,i6))


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    -I3m(s(i1,i2),s(i3,i4),s(i5,i6))*
     &    (-((za(i2,i3)*za(i3,i5)*zb(i1,i4)*zb(i4,i6))/
     &    (s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,i
     &    6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2))-
     &    (3._dp*(-s(i1,i2)+s(i3,i4)-s(i5,i6))*
     &    ((s(i1,i2)-s(i3,i4)-s(i5,i6))*za(i2,i5)*zb(i1,i2)+
     &    (-s(i1,i2)-s(i3,i4)+s(i5,i6))*za(i5,i6)*zb(i1,i6))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (3._dp*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*za(i2,i5)*zb(i1,i2)*zb(i1,i6))-
     &    za(i2,i3)*za(i4,i5)*zb(i1,i4)*zb(i3,i6)-
     &    za(i2,i4)*za(i3,i5)*zb(i1,i3)*zb(i4,i6)-
     &    za(i2,i5)*za(i5,i6)*zb(i1,i6)*zb(i5,i6)))/
     &    (2._dp*(s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s
     &    (i5,i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))-
     &    (za(i3,i5)*zb(i1,i4)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))
     &    /
     &    (2._dp*(-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))*t(i1,i2,i4))+
     &    (za(i4,i5)*zb(i1,i3)*(-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))*
     &    (-t(i1,i2,i3)+t(i1,i2,i4)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2))

      
        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    -Lnrat(-s(i1,i2),-s(i3,i4))*((-6._dp*zb(i1,i2)*
     &    ((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i5)-
     &    2._dp*za(i1,i2)*za(i5,i6)*zb(i1,i6))*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3)))+
     &    (za(i2,i3)*zb(i1,i4)*zb(i4,i6)**2)/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i3,i4)*zb(i5,i6))+
     &    (zb(i1,i4)*zb(i3,i6)*(-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6)))
     &    /
     &    ((-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2*zb(i3,i4)*
     &    zb(i5,i6))-(zb(i1,i3)*
     &    (-(za(i1,i3)*zb(i1,i4))-za(i2,i3)*zb(i2,i4))*
     &    (-(za(i1,i2)*zb(i1,i6))+za(i2,i4)*zb(i4,i6))*
     &    (3._dp*(-(za(i1,i4)*zb(i1,i3)*zb(i4,i6))-
     &    za(i2,i4)*zb(i2,i3)*zb(i4,i6))-
     &    zb(i3,i6)*(-t(i1,i2,i3)+t(i1,i2,i4))))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*
     &    (-(za(i1,i4)*zb(i1,i3))-za(i2,i4)*zb(i2,i3))**2*zb(i3,i4)*zb(i5,i6
     &    )
     &    ))


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    -Lnrat(-s(i5,i6),-s(i3,i4))*(-((za(i2,i3)**2*za(i3,i5)*zb(i4,i6)
     &    )/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)*za(i3,i4)))-
     &    (za(i2,i4)*za(i3,i5)*(za(i2,i3)*zb(i3,i6)+za(i2,i5)*zb(i5,i6)))/
     &    (za(i1,i2)*za(i3,i4)*(-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))*
     &    *
     &    2)-(6._dp*za(i5,i6)*(-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))*
     &    (za(i2,i3)*zb(i3,i6)+za(i2,i5)*zb(i5,i6))*
     &    (-((-s(i1,i2)+s(i3,i4)-s(i5,i6))*zb(i1,i6))+
     &    2._dp*za(i2,i5)*zb(i1,i2)*zb(i5,i6)))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &    (-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6)))-
     &    (za(i4,i5)*(-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))*
     &    (za(i2,i3)*zb(i3,i6)+za(i2,i5)*zb(i5,i6))*
     &    (3._dp*(za(i2,i3)*za(i4,i5)*zb(i3,i5)+
     &    za(i2,i3)*za(i4,i6)*zb(i3,i6))+
     &    za(i2,i4)*(t(i3,i5,i6)-t(i4,i5,i6))))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*za(i1,i2)*za(i3,i4)*
     &    (-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))**2))


        zaj_virt_a64ax_mp=zaj_virt_a64ax_mp
     &    -Lnrat(-s(i5,i6),-s(i3,i4))*(-((za(i3,i5)*zb(i1,i4)**2*zb(i4,i6))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i1,i2)*zb(i3,i4)))-
     &    (zb(i1,i3)*(za(i4,i5)*zb(i1,i4)-za(i5,i6)*zb(i1,i6))*zb(i4,i6))/
     &    (zb(i1,i2)*zb(i3,i4)*(-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))*
     &    *2)
     &    +(6._dp*(za(i4,i5)*zb(i1,i4)-za(i5,i6)*zb(i1,i6))*
     &    (-((-s(i1,i2)+s(i3,i4)-s(i5,i6))*za(i2,i5))+
     &    2._dp*za(i1,i2)*za(i5,i6)*zb(i1,i6))*
     &    (-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))*zb(i5,i6))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)**2*
     &    (-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6)))-
     &    ((za(i4,i5)*zb(i1,i4)-za(i5,i6)*zb(i1,i6))*zb(i3,i6)*
     &    (-(za(i3,i5)*zb(i4,i5))-za(i3,i6)*zb(i4,i6))*
     &    (3._dp*(za(i4,i5)*zb(i1,i4)*zb(i3,i5)+
     &    za(i4,i6)*zb(i1,i4)*zb(i3,i6))+
     &    zb(i1,i3)*(-t(i3,i5,i6)+t(i4,i5,i6))))/
     &    ((s(i1,i2)**2-2._dp*s(i1,i2)*s(i3,i4)+s(i3,i4)**2-2._dp*s(i1,i2)*s(i5,
     &    i6)-
     &    2._dp*s(i3,i4)*s(i5,i6)+s(i5,i6)**2)*zb(i1,i2)*zb(i3,i4)*
     &    (-(za(i4,i5)*zb(i3,i5))-za(i4,i6)*zb(i3,i6))**2))

      end function

      end module
