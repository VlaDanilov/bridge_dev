      module zajj_treeamps_m
        implicit none
        include 'types.f'
        include 'mxpart.f'

      contains

      ! amplitudes for 
      ! 0 -> q(p1) + qb(p2) + l(p3) + lb(p4) + γ(p5) + g(p6) + g(p7)

      ! suffix specifies helicities of (g,g,a)
      ! +(p) helicity corresponds to index 2 in the old Zgamjet code, -(m) to 1

      ! additionally to the anomalous coupling amplitudes I give the
      ! final state radiation amplitudes here which can serve as a check

      ! amplitudes obtained by contracting the qqgg V(mu) and qqQQ V(mu)
      ! production currents with V(mu) -> (lla) and anomalous coupling decay
      ! currents

      ! V production current from Berends, Giele, Kuijf, Nucl.Phys. B321 (1989)
      ! I think that the S(+;-;+;-) current on p.66 is broken and had to
      ! recalculate it myself
      ! (see also Bern, Forde, Kosower, Mastrolia, Phys.Rev. D72 (2005))


c =========== 
c lla decay amplitudes
c =========== 


      pure function zajj_tree_qqQQ_fsr_plus(i1,i2,i3,i4,i5,i6,i7,za,zb) ! + photon
        complex(dp) :: zajj_tree_qqQQ_fsr_plus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_fsr_plus =
     -  (zb(i5,i4)*(t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -       (-(za(i2,i7)*(za(i3,i4)*zb(i4,i2) + 
     -              za(i3,i5)*zb(i5,i2))) + 
     -         za(i6,i7)*(za(i3,i4)*zb(i6,i4) + za(i3,i5)*zb(i6,i5))) + 
     -      t(i2,i6,i7)*za(i1,i7)*
     -       (za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*
     -       (za(i1,i3)*zb(i6,i1) + za(i3,i7)*zb(i7,i6))))/
     -  (s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i5))

      end function

      pure function zajj_tree_qqQQ_fsr_minus(i1,i2,i3,i4,i5,i6,i7,za,zb) ! - photon
        complex(dp) :: zajj_tree_qqQQ_fsr_minus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_fsr_minus =
     -  (s(i3,i5)*za(i4,i5)*zb(i4,i3)*
     -     (t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i4,i2)) + za(i6,i7)*zb(i6,i4)) + 
     -       t(i2,i6,i7)*za(i1,i7)*zb(i4,i2)*
     -        (za(i1,i3)*zb(i6,i1) + za(i3,i7)*zb(i7,i6))) - 
     -    s(i4,i5)*za(i3,i5)*zb(i5,i3)*
     -     (t(i1,i6,i7)*za(i1,i5)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i4,i2)) + za(i6,i7)*zb(i6,i4)) + 
     -       t(i2,i6,i7)*za(i1,i7)*zb(i4,i2)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))))/
     -  (s(i3,i5)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3))

      end function


      pure function zajj_tree_qqgg_fsr_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 22222
        complex(dp) :: zajj_tree_qqgg_fsr_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_ppp =
     &  -((za(i1,i3)*(-(za(i1,i2)*
     &           (za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))) + 
     &        za(i1,i6)*(za(i3,i4)*zb(i6,i4) + za(i3,i5)*zb(i6,i5)) + 
     &        za(i1,i7)*(za(i3,i4)*zb(i7,i4) + za(i3,i5)*zb(i7,i5))))/
     &    (za(i1,i7)*za(i2,i6)*za(i3,i5)*za(i4,i5)*za(i6,i7)))

      end function

      pure function zajj_tree_qqgg_fsr_mpp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 21222
        complex(dp) :: zajj_tree_qqgg_fsr_mpp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_mpp =
     &  ((-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*
     &    (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4) + 
     &      za(i1,i7)*zb(i7,i4)))/
     &  (za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3)*zb(i5,i4))

      end function

      pure function zajj_tree_qqgg_fsr_pmm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 22112
        complex(dp) :: zajj_tree_qqgg_fsr_pmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_pmm =
     &  ((za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*
     &    (za(i1,i3)*zb(i2,i1) + za(i3,i6)*zb(i6,i2) + 
     &      za(i3,i7)*zb(i7,i2)))/
     &  (za(i3,i5)*za(i4,i5)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))

      end function

      pure function zajj_tree_qqgg_fsr_mmm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 21112
        complex(dp) :: zajj_tree_qqgg_fsr_mmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_mmm =
     &  (zb(i4,i2)*(za(i1,i3)*zb(i2,i1)*zb(i4,i3) - 
     &      za(i1,i5)*zb(i2,i1)*zb(i5,i4) + 
     &      za(i3,i6)*zb(i4,i3)*zb(i6,i2) - 
     &      za(i5,i6)*zb(i5,i4)*zb(i6,i2) + 
     &      za(i3,i7)*zb(i4,i3)*zb(i7,i2) - 
     &      za(i5,i7)*zb(i5,i4)*zb(i7,i2)))/
     &  (zb(i5,i3)*zb(i5,i4)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))

      end function

      pure function zajj_tree_qqgg_fsr_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 22212
        complex(dp) :: zajj_tree_qqgg_fsr_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_ppm =
     &  (-((za(i1,i7)*(za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*
     &         zb(i6,i1)*(za(i1,i3)*zb(i6,i1) + za(i3,i7)*zb(i7,i6)))/
     &       (za(i1,i6)*zb(i6,i1) + za(i1,i7)*zb(i7,i1) + 
     &         za(i6,i7)*zb(i7,i6))) + 
     &    ((za(i2,i7)*(za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2)) - 
     &         za(i6,i7)*(za(i3,i4)*zb(i6,i4) + za(i3,i5)*zb(i6,i5)))*
     &       (za(i3,i7)*zb(i7,i6)*
     &          (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &            za(i6,i7)*zb(i7,i6)) + 
     &         za(i1,i3)*(za(i2,i6)*zb(i6,i1)*zb(i6,i2) + 
     &            za(i2,i7)*(-(zb(i6,i2)*zb(i7,i1)) + 
     &               zb(i6,i1)*zb(i7,i2)) + 
     &            za(i6,i7)*zb(i6,i1)*zb(i7,i6))))/
     &     (za(i2,i6)*(za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &         za(i6,i7)*zb(i7,i6))))/
     &  (za(i3,i5)*za(i4,i5)*za(i6,i7)*zb(i7,i1)*zb(i7,i6))

      end function

      pure function zajj_tree_qqgg_fsr_mpm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 21212
        complex(dp) :: zajj_tree_qqgg_fsr_mpm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_mpm =
     &  ((za(i1,i7)*zb(i4,i2)*zb(i6,i1)*
     &       (-(za(i1,i3)*zb(i4,i3)*zb(i6,i1)) + 
     &         za(i1,i5)*zb(i5,i4)*zb(i6,i1) + 
     &         (-(za(i3,i7)*zb(i4,i3)) + za(i5,i7)*zb(i5,i4))*zb(i7,i6))
     &       )/
     &     (za(i1,i6)*zb(i6,i1) + za(i1,i7)*zb(i7,i1) + 
     &       za(i6,i7)*zb(i7,i6)) + 
     &    ((za(i2,i7)*zb(i4,i2) - za(i6,i7)*zb(i6,i4))*
     &       ((za(i3,i7)*zb(i4,i3) - za(i5,i7)*zb(i5,i4))*zb(i7,i6)*
     &          (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &            za(i6,i7)*zb(i7,i6)) + 
     &         za(i1,i3)*zb(i4,i3)*
     &          (za(i2,i6)*zb(i6,i1)*zb(i6,i2) + 
     &            za(i2,i7)*(-(zb(i6,i2)*zb(i7,i1)) + 
     &               zb(i6,i1)*zb(i7,i2)) + 
     &            za(i6,i7)*zb(i6,i1)*zb(i7,i6)) - 
     &         za(i1,i5)*zb(i5,i4)*
     &          (za(i2,i6)*zb(i6,i1)*zb(i6,i2) + 
     &            za(i2,i7)*(-(zb(i6,i2)*zb(i7,i1)) + 
     &               zb(i6,i1)*zb(i7,i2)) + 
     &            za(i6,i7)*zb(i6,i1)*zb(i7,i6))))/
     &     (za(i2,i6)*(za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &         za(i6,i7)*zb(i7,i6))))/
     &  (za(i6,i7)*zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*zb(i7,i6))

      end function


      ! these two are probably broken in the 1989 paper, so I had to
      ! recompute them
      
      pure function zajj_tree_qqgg_fsr_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 22122
        complex(dp) :: zajj_tree_qqgg_fsr_pmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_fsr_pmp =
     -  -((zb(i5,i4)*(-(t(i2,i6,i7)*za(i1,i6)*
     -           (za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*zb(i7,i1)*
     -           (s(i2,i6)*za(i1,i6)*
     -              (s(i6,i7)*za(i3,i6) - za(i1,i3)*za(i6,i7)*zb(i7,i1))
     -               + t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*za(i6,i7)*
     -              zb(i7,i2))*zb(i7,i6)) + 
     -        s(i1,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*za(i6,i7)*
     -         zb(i7,i2)**2*(s(i6,i7)*
     -            (za(i3,i4)*zb(i7,i4) + za(i3,i5)*zb(i7,i5)) + 
     -           za(i2,i6)*(za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*
     -            zb(i7,i6))))/
     -    (s(i1,i7)*s(i2,i6)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -      za(i3,i5)*za(i6,i7)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_fsr_mmp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! 21122
        complex(dp) :: zajj_tree_qqgg_fsr_mmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zajj_tree_qqgg_fsr_mmp =
     &  (za(i1,i6)**2*zb(i4,i2)*
     &     (za(i3,i6)*zb(i4,i3) - za(i5,i6)*zb(i5,i4))*zb(i6,i2)*
     &     zb(i7,i6)*(za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &       za(i6,i7)*zb(i7,i6)) + 
     &    za(i1,i5)*zb(i5,i4)*
     &     (za(i1,i7)*zb(i7,i2)**2*
     &        (za(i2,i6)*zb(i4,i2) + za(i6,i7)*zb(i7,i4))*
     &        (za(i1,i7)*zb(i7,i1) + za(i6,i7)*zb(i7,i6)) + 
     &       za(i1,i6)**2*zb(i4,i2)*
     &        (zb(i6,i2)*zb(i7,i1) - zb(i6,i1)*zb(i7,i2))*
     &        (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &          za(i6,i7)*zb(i7,i6)) - 
     &       za(i1,i6)*zb(i7,i2)*
     &        (za(i6,i7)*zb(i4,i2)*zb(i7,i6)*
     &           (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &             za(i6,i7)*zb(i7,i6)) + 
     &          za(i1,i7)*(za(i2,i7)*zb(i4,i2)*zb(i7,i1)*zb(i7,i2) + 
     &             za(i2,i6)*zb(i4,i2)*
     &              (zb(i6,i2)*zb(i7,i1) - zb(i6,i1)*zb(i7,i2)) + 
     &             za(i6,i7)*
     &              (-(zb(i6,i1)*zb(i7,i2)*zb(i7,i4)) + 
     &                zb(i4,i2)*zb(i7,i1)*zb(i7,i6))))) + 
     &    za(i1,i3)*zb(i4,i3)*
     &     (-(za(i1,i7)*zb(i7,i2)**2*
     &          (za(i2,i6)*zb(i4,i2) + za(i6,i7)*zb(i7,i4))*
     &          (za(i1,i7)*zb(i7,i1) + za(i6,i7)*zb(i7,i6))) - 
     &       za(i1,i6)**2*zb(i4,i2)*
     &        (zb(i6,i2)*zb(i7,i1) - zb(i6,i1)*zb(i7,i2))*
     &        (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &          za(i6,i7)*zb(i7,i6)) + 
     &       za(i1,i6)*zb(i7,i2)*
     &        (za(i6,i7)*zb(i4,i2)*zb(i7,i6)*
     &           (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &             za(i6,i7)*zb(i7,i6)) + 
     &          za(i1,i7)*(za(i2,i7)*zb(i4,i2)*zb(i7,i1)*zb(i7,i2) + 
     &             za(i2,i6)*zb(i4,i2)*
     &              (zb(i6,i2)*zb(i7,i1) - zb(i6,i1)*zb(i7,i2)) + 
     &             za(i6,i7)*
     &              (-(zb(i6,i1)*zb(i7,i2)*zb(i7,i4)) + 
     &                zb(i4,i2)*zb(i7,i1)*zb(i7,i6))))))/
     &  (za(i1,i7)*za(i6,i7)*zb(i5,i3)*zb(i5,i4)*zb(i6,i2)*zb(i7,i6)*
     &    (za(i1,i6)*zb(i6,i1) + za(i1,i7)*zb(i7,i1) + 
     &      za(i6,i7)*zb(i7,i6))*
     &    (za(i2,i6)*zb(i6,i2) + za(i2,i7)*zb(i7,i2) + 
     &      za(i6,i7)*zb(i7,i6)))

      end function

c ================= 
c the arguments ha,hb for the anomalous coupling amplitudes are
c ha = I*h1tZ + h3tZ, hb = I*h2tZ + h4tZ for the ZZ amps
c and ha = I*h1tgam + h3tgam, hb = I*h2tgam + h4tgam for the aZ and Za amps
c ================= 

c =========== 
c ZZ anomalous coupling decay amplitudes
c =========== 

      pure function zajj_tree_qqQQ_anomZZ_plus(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! + photon
        complex(dp) :: zajj_tree_qqQQ_anomZZ_plus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_anomZZ_plus = -(
     -  (ha*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (-2*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*zb(i5,i2)*zb(i6,i1) + 
     -       2*t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) - 
     -       2*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*zb(i5,i2)*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)) + 
     -  (hb*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (s(i3,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*zb(i5,i2)*
     -        zb(i6,i1) + s(i4,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*
     -        zb(i5,i2)*zb(i6,i1) + 
     -       t(i2,i6,i7)*za(i1,i5)*za(i1,i7)*za(i3,i4)*zb(i5,i2)*
     -        zb(i5,i4)*zb(i6,i1) - 
     -       s(i3,i5)*t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) - 
     -       s(i4,i5)*t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) - 
     -       t(i1,i6,i7)*za(i1,i5)*za(i3,i4)*zb(i5,i4)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) + 
     -       s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*zb(i5,i2)*
     -        zb(i7,i6) + s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*
     -        zb(i5,i2)*zb(i7,i6) + 
     -       t(i2,i6,i7)*za(i1,i7)*za(i3,i4)*za(i5,i7)*zb(i5,i2)*
     -        zb(i5,i4)*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)))

      end function

      pure function zajj_tree_qqQQ_anomZZ_minus(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! - photon
        complex(dp) :: zajj_tree_qqQQ_anomZZ_minus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_anomZZ_minus =
     -  -(ha*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -      (-2*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i4,i3)*zb(i5,i2)*
     -         zb(i6,i2) - 2*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i3,i2)*
     -         zb(i5,i4)*zb(i6,i2) + 
     -        2*t(i1,i6,i7)*za(i1,i5)*za(i6,i7)*zb(i6,i2)*
     -         (zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)) + 
     -        2*t(i2,i6,i7)*za(i1,i7)*zb(i4,i3)*zb(i5,i2)*
     -         (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -        2*t(i2,i6,i7)*za(i1,i7)*zb(i3,i2)*zb(i5,i4)*
     -         (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3)) - 
     -  (hb*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i3,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i2) + 
     -       s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i2) + 
     -       t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*za(i4,i5)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i5,i4)*zb(i6,i2) - 
     -       s(i3,i5)*t(i1,i6,i7)*za(i1,i5)*za(i6,i7)*zb(i5,i4)*
     -        zb(i6,i2)*zb(i6,i3) - 
     -       s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i6,i7)*zb(i6,i2)*
     -        (zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)) - 
     -       s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -       t(i2,i6,i7)*za(i1,i7)*za(i4,i5)*zb(i4,i3)*zb(i5,i2)*
     -        zb(i5,i4)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3))

      end function

      pure function zajj_tree_qqgg_anomZZ_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22222
        complex(dp) :: zajj_tree_qqgg_anomZZ_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_ppp = -(
     -  (ha*(s(i3,i5) + s(i4,i5))*za(i1,i3)*zb(i5,i4)*
     -     (-(za(i1,i2)*zb(i5,i2)) + za(i1,i6)*zb(i6,i5) + 
     -       za(i1,i7)*zb(i7,i5)))/(2.*za(i1,i7)*za(i2,i6)*za(i6,i7)) - 
     -  (hb*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (s(i3,i5)*za(i1,i3) + s(i4,i5)*za(i1,i3) + 
     -       za(i1,i5)*za(i3,i4)*zb(i5,i4))*
     -     (-(za(i1,i2)*zb(i5,i2)) + za(i1,i6)*zb(i6,i5) + 
     -       za(i1,i7)*zb(i7,i5)))/(4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)))

      end function

      pure function zajj_tree_qqgg_anomZZ_mpp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21222
        complex(dp) :: zajj_tree_qqgg_anomZZ_mpp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_mpp =
     -  -(ha*(s(i3,i5) + s(i4,i5))*za(i1,i5)*za(i3,i5)*
     -      (2*za(i1,i2)*zb(i4,i3)*zb(i5,i2) + 
     -        2*za(i1,i2)*zb(i3,i2)*zb(i5,i4) - 
     -        2*za(i1,i6)*zb(i5,i4)*zb(i6,i3) - 
     -        2*za(i1,i6)*zb(i4,i3)*zb(i6,i5) - 
     -        2*za(i1,i7)*zb(i5,i4)*zb(i7,i3) - 
     -        2*za(i1,i7)*zb(i4,i3)*zb(i7,i5)))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3)) - 
     -  (hb*(s(i3,i5) + s(i4,i5))*za(i1,i5)*za(i3,i5)*
     -     (-(s(i4,i5)*za(i1,i2)*zb(i4,i3)*zb(i5,i2)) - 
     -       s(i3,i5)*za(i1,i2)*zb(i3,i2)*zb(i5,i4) - 
     -       s(i4,i5)*za(i1,i2)*zb(i3,i2)*zb(i5,i4) + 
     -       s(i3,i5)*za(i1,i6)*zb(i5,i4)*zb(i6,i3) + 
     -       s(i4,i5)*za(i1,i6)*zb(i5,i4)*zb(i6,i3) + 
     -       za(i1,i6)*za(i4,i5)*zb(i4,i3)*zb(i5,i4)*zb(i6,i5) + 
     -       s(i3,i5)*za(i1,i7)*zb(i5,i4)*zb(i7,i3) + 
     -       s(i4,i5)*za(i1,i7)*zb(i5,i4)*zb(i7,i3) + 
     -       za(i1,i7)*za(i4,i5)*zb(i4,i3)*zb(i5,i4)*zb(i7,i5)))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3))

      end function

      pure function zajj_tree_qqgg_anomZZ_pmm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22112
        complex(dp) :: zajj_tree_qqgg_anomZZ_pmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_pmm = -(
     -  (ha*(s(i3,i5) + s(i4,i5))*zb(i5,i2)*zb(i5,i4)*
     -     (-2*za(i1,i3)*zb(i2,i1) - 2*za(i3,i6)*zb(i6,i2) - 
     -       2*za(i3,i7)*zb(i7,i2)))/(4._dp*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))
     -   + (hb*(s(i3,i5) + s(i4,i5))*zb(i5,i2)*zb(i5,i4)*
     -     (s(i3,i5)*za(i1,i3)*zb(i2,i1) + 
     -       s(i4,i5)*za(i1,i3)*zb(i2,i1) + 
     -       za(i1,i5)*za(i3,i4)*zb(i2,i1)*zb(i5,i4) + 
     -       s(i3,i5)*za(i3,i6)*zb(i6,i2) + 
     -       s(i4,i5)*za(i3,i6)*zb(i6,i2) + 
     -       za(i3,i4)*za(i5,i6)*zb(i5,i4)*zb(i6,i2) + 
     -       s(i3,i5)*za(i3,i7)*zb(i7,i2) + 
     -       s(i4,i5)*za(i3,i7)*zb(i7,i2) + 
     -       za(i3,i4)*za(i5,i7)*zb(i5,i4)*zb(i7,i2)))/
     -   (4._dp*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_anomZZ_mmm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21112
        complex(dp) :: zajj_tree_qqgg_anomZZ_mmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_mmm =
     -  (ha*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (-2*zb(i4,i3)*zb(i5,i2) - 2*zb(i3,i2)*zb(i5,i4))*
     -     (za(i1,i5)*zb(i2,i1) + za(i5,i6)*zb(i6,i2) + 
     -       za(i5,i7)*zb(i7,i2)))/
     -   (4._dp*zb(i5,i3)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)) + 
     -  (hb*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i4,i5)*zb(i4,i3)*zb(i5,i2) + 
     -       s(i3,i5)*zb(i3,i2)*zb(i5,i4) + s(i4,i5)*zb(i3,i2)*zb(i5,i4)
     -       )*(za(i1,i5)*zb(i2,i1) + za(i5,i6)*zb(i6,i2) + 
     -       za(i5,i7)*zb(i7,i2)))/
     -   (4._dp*zb(i5,i3)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))

      end function

      pure function zajj_tree_qqgg_anomZZ_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22212
        complex(dp) :: zajj_tree_qqgg_anomZZ_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_ppm = -(
     -  -(ha*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -      (-2*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*zb(i5,i2)*
     -         zb(i6,i1)**2 - 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*zb(i6,i1)*
     -         (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5)) + 
     -        2*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i6,i2)*
     -         (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i1)
     -         - 2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i7)*zb(i5,i2)*
     -         zb(i6,i1)*zb(i7,i6) - 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*
     -         (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6))
     -      )/(4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i7,i1))
     -   - (hb*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (s(i3,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*
     -        zb(i5,i2)*zb(i6,i1)**2 + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*
     -        zb(i5,i2)*zb(i6,i1)**2 + 
     -       t(i2,i6,i7)*za(i1,i5)*za(i1,i7)*za(i2,i6)*za(i3,i4)*
     -        zb(i5,i2)*zb(i5,i4)*zb(i6,i1)**2 + 
     -       s(i3,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*zb(i6,i1)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5)) + 
     -       s(i4,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*zb(i6,i1)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5)) + 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i5)*za(i3,i4)*zb(i5,i4)*
     -        zb(i6,i1)*(-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))
     -        - s(i3,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i1)
     -        - s(i4,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i1)
     -        - t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*za(i3,i4)*zb(i5,i4)*
     -        zb(i6,i2)*(-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*
     -        zb(i7,i1) + s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*
     -        za(i3,i7)*zb(i5,i2)*zb(i6,i1)*zb(i7,i6) + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i7)*
     -        zb(i5,i2)*zb(i6,i1)*zb(i7,i6) + 
     -       t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i4)*za(i5,i7)*
     -        zb(i5,i2)*zb(i5,i4)*zb(i6,i1)*zb(i7,i6) + 
     -       s(i3,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6)
     -        + s(i4,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6)
     -        + t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i4)*za(i5,i7)*zb(i5,i4)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6)))
     -    /(4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i7,i1)))

      end function

      pure function zajj_tree_qqgg_anomZZ_mpm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21212
        complex(dp) :: zajj_tree_qqgg_anomZZ_mpm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_mpm =
     -  -(ha*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -      (-2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i4,i3)*zb(i5,i2)*
     -         zb(i6,i1)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -        2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i3,i2)*zb(i5,i4)*
     -         zb(i6,i1)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -        2*t(i1,i6,i7)*za(i2,i7)*zb(i4,i3)*zb(i5,i2)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           ) + 2*t(i1,i6,i7)*za(i2,i7)*zb(i3,i2)*zb(i5,i4)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           ) - 2*t(i1,i6,i7)*za(i6,i7)*zb(i5,i4)*zb(i6,i3)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           ) - 2*t(i1,i6,i7)*za(i6,i7)*zb(i4,i3)*zb(i6,i5)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           )))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i5,i3)*
     -     zb(i7,i1)) - (hb*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -       s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -       s(i4,i5)*t(i1,i6,i7)*za(i2,i7)*zb(i4,i3)*zb(i5,i2)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         - s(i3,i5)*t(i1,i6,i7)*za(i2,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         - s(i4,i5)*t(i1,i6,i7)*za(i2,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         + s(i3,i5)*t(i1,i6,i7)*za(i6,i7)*zb(i5,i4)*zb(i6,i3)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         + s(i4,i5)*t(i1,i6,i7)*za(i6,i7)*zb(i5,i4)*zb(i6,i3)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         + t(i1,i6,i7)*za(i4,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i4)*
     -        zb(i6,i5)*(-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -       ))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i5,i3)*
     -     zb(i7,i1))

      end function


      pure function zajj_tree_qqgg_anomZZ_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22122
        complex(dp) :: zajj_tree_qqgg_anomZZ_pmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_pmp = -(
     -  -(ha*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -      (2*s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i6)*
     -         zb(i5,i2)*zb(i7,i1)*zb(i7,i6) - 
     -        2*s(i2,i6)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*za(i6,i7)*
     -         zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) + 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)*za(i2,i6)*
     -         za(i6,i7)*zb(i5,i2)*zb(i7,i1)*zb(i7,i2)*zb(i7,i6) - 
     -        2*s(i1,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*za(i6,i7)*
     -         zb(i7,i2)**2*(s(i6,i7)*zb(i7,i5) + 
     -           za(i2,i6)*zb(i5,i2)*zb(i7,i6))))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i7,i6)) - 
     -  (hb*(s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (-(s(i2,i6)*s(i3,i5)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*
     -          za(i3,i6)*zb(i5,i2)*zb(i7,i1)*zb(i7,i6)) - 
     -       s(i2,i6)*s(i4,i5)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*
     -        za(i3,i6)*zb(i5,i2)*zb(i7,i1)*zb(i7,i6) - 
     -       s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i4)*
     -        za(i5,i6)*zb(i5,i2)*zb(i5,i4)*zb(i7,i1)*zb(i7,i6) + 
     -       s(i2,i6)*s(i3,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*
     -        za(i6,i7)*zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) + 
     -       s(i2,i6)*s(i4,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*
     -        za(i6,i7)*zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) + 
     -       s(i2,i6)*t(i2,i6,i7)*za(i1,i5)*za(i1,i6)**2*za(i3,i4)*
     -        za(i6,i7)*zb(i5,i2)*zb(i5,i4)*zb(i7,i1)**2*zb(i7,i6) - 
     -       s(i3,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)*
     -        za(i2,i6)*za(i6,i7)*zb(i5,i2)*zb(i7,i1)*zb(i7,i2)*
     -        zb(i7,i6) - s(i4,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*
     -        za(i1,i6)*za(i2,i6)*za(i6,i7)*zb(i5,i2)*zb(i7,i1)*
     -        zb(i7,i2)*zb(i7,i6) - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i5)*za(i1,i6)*za(i2,i6)*
     -        za(i3,i4)*za(i6,i7)*zb(i5,i2)*zb(i5,i4)*zb(i7,i1)*
     -        zb(i7,i2)*zb(i7,i6) + 
     -       s(i1,i7)*s(i3,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*
     -        za(i6,i7)*zb(i7,i2)**2*
     -        (s(i6,i7)*zb(i7,i5) + za(i2,i6)*zb(i5,i2)*zb(i7,i6)) + 
     -       s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*
     -        za(i6,i7)*zb(i7,i2)**2*
     -        (s(i6,i7)*zb(i7,i5) + za(i2,i6)*zb(i5,i2)*zb(i7,i6)) + 
     -       s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i3,i4)*
     -        za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*
     -        (s(i6,i7)*zb(i7,i5) + za(i2,i6)*zb(i5,i2)*zb(i7,i6))))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_anomZZ_mmp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21122
        complex(dp) :: zajj_tree_qqgg_anomZZ_mmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZZ_mmp =
     -  (ha*za(i3,i5)*(2*s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*
     -        za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i3) + 2*s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i3,i5)*za(i6,i7)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i3) + 
     -       2*s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i6,i7)*zb(i4,i3)*zb(i7,i2)**2*zb(i7,i5) + 
     -       2*s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i2)**2*
     -        zb(i7,i5) + 2*s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)**2*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*zb(i7,i2)**2*
     -        zb(i7,i6) + 2*s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*zb(i5,i3)*
     -        zb(i7,i2)**2*zb(i7,i6) + 
     -       2*s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i6) + 
     -       2*s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*za(i3,i5)*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i6) - 2*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i5,i3)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*zb(i3,i2)*
     -        zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6)))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i5,i3)*zb(i7,i6)) + 
     -  (hb*za(i3,i5)*(-(s(i1,i7)*s(i3,i5)*s(i6,i7)*t(i1,i6,i7)*
     -          za(i1,i5)*za(i2,i6)*za(i3,i5)*za(i6,i7)*zb(i5,i3)*
     -          zb(i5,i4)*zb(i7,i2)**2*zb(i7,i3)) - 
     -       2*s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i3,i5)*za(i6,i7)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i3) - 
     -       s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i2)**2*
     -        zb(i7,i5) - s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)**2*za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*
     -        zb(i5,i3)*zb(i7,i2)**2*zb(i7,i6) - 
     -       s(i1,i7)*s(i3,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i6) - 
     -       2*s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i6) + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i5,i3)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i3,i5)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*
     -        zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + 2*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*
     -        zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i4,i5)**2*t(i2,i6,i7)*za(i1,i6)*
     -        (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - s(i1,i7)*s(i4,i5)**2*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i6,i7)*zb(i7,i2)**2*
     -        (s(i6,i7)*(zb(i5,i4)*zb(i7,i3) + zb(i4,i3)*zb(i7,i5)) + 
     -          za(i2,i6)*(zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -           zb(i7,i6))))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i5,i3)*zb(i7,i6))

      end function

c =========== 
c Za anomalous coupling decay amplitudes
c =========== 

      pure function zajj_tree_qqQQ_anomZa_plus(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! + photon
        complex(dp) :: zajj_tree_qqQQ_anomZa_plus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_anomZa_plus = -(
     -  (ha*s(i3,i4)*zb(i5,i4)*
     -     (2*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*zb(i5,i2)*zb(i6,i1) - 
     -       2*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i5,i2)*zb(i6,i2) + 
     -       2*t(i1,i6,i7)*za(i1,i3)*za(i6,i7)*zb(i6,i2)*zb(i6,i5) + 
     -       2*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*zb(i5,i2)*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)) + 
     -  (hb*s(i3,i4)*zb(i5,i4)*
     -     (t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i3,i4)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i1) + 
     -       t(i2,i6,i7)*za(i1,i4)*za(i1,i7)*za(i3,i4)*zb(i4,i2)*
     -        zb(i5,i4)*zb(i6,i1) - 
     -       t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*za(i3,i4)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i2) + 
     -       t(i1,i6,i7)*za(i1,i3)*za(i3,i4)*za(i6,i7)*zb(i5,i4)*
     -        zb(i6,i2)*zb(i6,i3) + 
     -       t(i1,i6,i7)*za(i1,i4)*za(i3,i4)*zb(i5,i4)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i4,i2)) + za(i6,i7)*zb(i6,i4)) + 
     -       t(i2,i6,i7)*za(i1,i7)*za(i3,i4)*za(i3,i7)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i7,i6) + 
     -       t(i2,i6,i7)*za(i1,i7)*za(i3,i4)*za(i4,i7)*zb(i4,i2)*
     -        zb(i5,i4)*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)))

      end function

      pure function zajj_tree_qqQQ_anomZa_minus(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! - photon
        complex(dp) :: zajj_tree_qqQQ_anomZa_minus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_anomZa_minus =
     -  (ha*s(i3,i4)*za(i3,i5)*
     -     (t(i1,i6,i7)*za(i1,i5)*zb(i6,i2)*
     -        (-(za(i2,i7)*(zb(i4,i3)*zb(i5,i2) + 
     -               zb(i3,i2)*zb(i5,i4))) + 
     -          za(i6,i7)*(zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)))
     -        + t(i2,i6,i7)*za(i1,i7)*
     -        (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))))/
     -   (2._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3)) - 
     -  (hb*s(i3,i4)*s(i3,i5)*za(i3,i5)*zb(i4,i3)*
     -     (t(i1,i6,i7)*zb(i6,i2)*
     -        (za(i1,i3)*(-(za(i2,i7)*zb(i3,i2)) + 
     -             za(i6,i7)*zb(i6,i3)) + 
     -          za(i1,i4)*(-(za(i2,i7)*zb(i4,i2)) + za(i6,i7)*zb(i6,i4))
     -          ) + t(i2,i6,i7)*za(i1,i7)*
     -        (za(i1,i3)*zb(i3,i2)*zb(i6,i1) + 
     -          za(i1,i4)*zb(i4,i2)*zb(i6,i1) + 
     -          (za(i3,i7)*zb(i3,i2) + za(i4,i7)*zb(i4,i2))*zb(i7,i6))))
     -    /(4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3))

      end function


      pure function zajj_tree_qqgg_anomZa_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22222
        complex(dp) :: zajj_tree_qqgg_anomZa_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_ppp = -(
     -  -(hb*s(i3,i4)*zb(i5,i4)*
     -      (-(za(i1,i2)*za(i1,i3)*za(i3,i4)*zb(i3,i2)*zb(i5,i4)) - 
     -        za(i1,i2)*za(i1,i4)*za(i3,i4)*zb(i4,i2)*zb(i5,i4) + 
     -        za(i1,i3)*za(i1,i6)*za(i3,i4)*zb(i5,i4)*zb(i6,i3) + 
     -        za(i1,i3)*za(i1,i7)*za(i3,i4)*zb(i5,i4)*zb(i7,i3) + 
     -        za(i1,i4)*za(i3,i4)*zb(i5,i4)*
     -         (za(i1,i6)*zb(i6,i4) + za(i1,i7)*zb(i7,i4))))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)) - 
     -  (ha*s(i3,i4)*zb(i5,i4)*
     -     (-2*za(i1,i2)*za(i1,i3)*zb(i5,i2) + 
     -       2*za(i1,i3)*za(i1,i6)*zb(i6,i5) + 
     -       2*za(i1,i3)*za(i1,i7)*zb(i7,i5)))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)))

      end function

      pure function zajj_tree_qqgg_anomZa_mpp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21222
        complex(dp) :: zajj_tree_qqgg_anomZa_mpp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_mpp =
     -  (hb*s(i3,i4)*s(i3,i5)*za(i3,i5)*zb(i4,i3)*
     -     (-(za(i1,i2)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))) + 
     -       za(i1,i3)*(za(i1,i6)*zb(i6,i3) + za(i1,i7)*zb(i7,i3)) + 
     -       za(i1,i4)*(za(i1,i6)*zb(i6,i4) + za(i1,i7)*zb(i7,i4))))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3)) + 
     -  (ha*s(i3,i4)*za(i1,i5)*za(i3,i5)*
     -     (za(i1,i2)*(zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4)) - 
     -       za(i1,i6)*(zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)) - 
     -       za(i1,i7)*(zb(i5,i4)*zb(i7,i3) + zb(i4,i3)*zb(i7,i5))))/
     -   (2.*za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3))

      end function

      pure function zajj_tree_qqgg_anomZa_pmm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22112
        complex(dp) :: zajj_tree_qqgg_anomZa_pmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_pmm = -(
     -  (ha*s(i3,i4)*zb(i5,i4)*
     -     (2*za(i1,i3)*zb(i2,i1)*zb(i5,i2) + 
     -       2*za(i3,i6)*zb(i5,i2)*zb(i6,i2) + 
     -       2*za(i3,i7)*zb(i5,i2)*zb(i7,i2)))/
     -   (4._dp*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)) + 
     -  (hb*s(i3,i4)*zb(i5,i4)*
     -     (za(i1,i3)*za(i3,i4)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4) + 
     -       za(i1,i4)*za(i3,i4)*zb(i2,i1)*zb(i4,i2)*zb(i5,i4) + 
     -       za(i3,i4)*za(i3,i6)*zb(i3,i2)*zb(i5,i4)*zb(i6,i2) + 
     -       za(i3,i4)*za(i4,i6)*zb(i4,i2)*zb(i5,i4)*zb(i6,i2) + 
     -       za(i3,i4)*za(i3,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i2) + 
     -       za(i3,i4)*za(i4,i7)*zb(i4,i2)*zb(i5,i4)*zb(i7,i2)))/
     -   (4._dp*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_anomZa_mmm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21112
        complex(dp) :: zajj_tree_qqgg_anomZa_mmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_mmm =
     -  (ha*s(i3,i4)*za(i3,i5)*
     -     (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -     (za(i1,i5)*zb(i2,i1) + za(i5,i6)*zb(i6,i2) + 
     -       za(i5,i7)*zb(i7,i2)))/
     -   (2.*zb(i5,i3)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)) - 
     -  (hb*s(i3,i4)*s(i3,i5)*za(i3,i5)*zb(i4,i3)*
     -     (za(i1,i3)*zb(i2,i1)*zb(i3,i2) + 
     -       za(i1,i4)*zb(i2,i1)*zb(i4,i2) + 
     -       za(i3,i6)*zb(i3,i2)*zb(i6,i2) + 
     -       za(i4,i6)*zb(i4,i2)*zb(i6,i2) + 
     -       za(i3,i7)*zb(i3,i2)*zb(i7,i2) + 
     -       za(i4,i7)*zb(i4,i2)*zb(i7,i2)))/
     -   (4._dp*zb(i5,i3)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))

      end function

      pure function zajj_tree_qqgg_anomZa_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22212
        complex(dp) :: zajj_tree_qqgg_anomZa_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_ppm = -(
     -  (hb*s(i3,i4)*zb(i5,i4)*
     -     (t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i2,i7)*za(i3,i4)*
     -        zb(i3,i2)*zb(i5,i4)*zb(i6,i1) - 
     -       t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*za(i3,i4)*
     -        zb(i3,i2)*zb(i5,i4)*zb(i6,i1)**2 - 
     -       t(i2,i6,i7)*za(i1,i4)*za(i1,i7)*za(i2,i6)*za(i3,i4)*
     -        zb(i4,i2)*zb(i5,i4)*zb(i6,i1)**2 - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i3,i4)*za(i6,i7)*
     -        zb(i5,i4)*zb(i6,i1)*zb(i6,i3) + 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i4)*za(i3,i4)*zb(i5,i4)*
     -        zb(i6,i1)*(za(i2,i7)*zb(i4,i2) - za(i6,i7)*zb(i6,i4)) - 
     -       t(i1,i6,i7)*za(i1,i3)*za(i2,i7)**2*za(i3,i4)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i2)*zb(i7,i1) + 
     -       t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*za(i3,i4)*za(i6,i7)*
     -        zb(i5,i4)*zb(i6,i2)*zb(i6,i3)*zb(i7,i1) + 
     -       t(i1,i6,i7)*za(i1,i4)*za(i2,i7)*za(i3,i4)*zb(i5,i4)*
     -        zb(i6,i2)*(-(za(i2,i7)*zb(i4,i2)) + za(i6,i7)*zb(i6,i4))*
     -        zb(i7,i1) + t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i7)*za(i3,i4)*
     -        za(i3,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i6) + 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i7)*za(i3,i4)*za(i4,i7)*
     -        zb(i4,i2)*zb(i5,i4)*zb(i7,i6) - 
     -       t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i4)*za(i3,i7)*
     -        zb(i3,i2)*zb(i5,i4)*zb(i6,i1)*zb(i7,i6) - 
     -       t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i4)*za(i4,i7)*
     -        zb(i4,i2)*zb(i5,i4)*zb(i6,i1)*zb(i7,i6) - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i4)*za(i6,i7)*zb(i5,i4)*
     -        (za(i3,i7)*zb(i6,i3) + za(i4,i7)*zb(i6,i4))*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i7,i1)) + 
     -  (ha*s(i3,i4)*zb(i5,i4)*
     -     (2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i5,i2)*
     -        zb(i6,i1) - 2*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*
     -        zb(i5,i2)*zb(i6,i1)**2 - 
     -       2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i6,i7)*zb(i6,i1)*
     -        zb(i6,i5) - 2*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)**2*
     -        zb(i5,i2)*zb(i6,i2)*zb(i7,i1) + 
     -       2*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*za(i6,i7)*zb(i6,i2)*
     -        zb(i6,i5)*zb(i7,i1) + 
     -       2*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i7)*za(i3,i7)*zb(i5,i2)*
     -        zb(i7,i6) - 2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i7)*
     -        zb(i5,i2)*zb(i6,i1)*zb(i7,i6) - 
     -       2*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*za(i6,i7)*zb(i6,i5)*
     -        zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i7,i1)))

      end function

      pure function zajj_tree_qqgg_anomZa_mpm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21212
        complex(dp) :: zajj_tree_qqgg_anomZa_mpm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_mpm =
     -  (hb*s(i3,i4)*s(i3,i5)*za(i3,i5)*zb(i4,i3)*
     -     (t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i6,i1)*
     -        (za(i1,i3)*zb(i3,i2)*zb(i6,i1) + 
     -          za(i1,i4)*zb(i4,i2)*zb(i6,i1) + 
     -          (za(i3,i7)*zb(i3,i2) + za(i4,i7)*zb(i4,i2))*zb(i7,i6))
     -        + t(i1,i6,i7)*(za(i2,i7)*zb(i6,i2)*
     -           (za(i1,i3)*(za(i2,i7)*zb(i3,i2) - 
     -                za(i6,i7)*zb(i6,i3)) + 
     -             za(i1,i4)*(za(i2,i7)*zb(i4,i2) - za(i6,i7)*zb(i6,i4))
     -             )*zb(i7,i1) + 
     -          t(i2,i6,i7)*(s(i6,i7)*za(i3,i7)*zb(i6,i3) + 
     -             za(i1,i3)*zb(i6,i1)*
     -              (-(za(i2,i7)*zb(i3,i2)) + za(i6,i7)*zb(i6,i3)) + 
     -             s(i6,i7)*za(i4,i7)*zb(i6,i4) + 
     -             za(i1,i4)*zb(i6,i1)*
     -              (-(za(i2,i7)*zb(i4,i2)) + za(i6,i7)*zb(i6,i4)) - 
     -             za(i2,i7)*za(i3,i7)*zb(i3,i2)*zb(i7,i6) - 
     -             za(i2,i7)*za(i4,i7)*zb(i4,i2)*zb(i7,i6)))))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i5,i3)*
     -     zb(i7,i1)) - (ha*s(i3,i4)*za(i3,i5)*
     -     (t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*
     -        (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -       t(i1,i6,i7)*(za(i1,i5)*za(i2,i7)*zb(i6,i2)*
     -           (za(i2,i7)*(zb(i4,i3)*zb(i5,i2) + 
     -                zb(i3,i2)*zb(i5,i4)) - 
     -             za(i6,i7)*(zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5))
     -             )*zb(i7,i1) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1)*
     -              (-(za(i2,i7)*
     -                   (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))) + 
     -                za(i6,i7)*
     -                 (zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5))) + 
     -             za(i5,i7)*
     -              (s(i6,i7)*
     -                 (zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)) - 
     -                za(i2,i7)*
     -                 (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -                 zb(i7,i6))))))/
     -   (2.*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i5,i3)*
     -     zb(i7,i1))

      end function

      pure function zajj_tree_qqgg_anomZa_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22122
        complex(dp) :: zajj_tree_qqgg_anomZa_pmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_pmp = -(
     -  -(ha*s(i3,i4)*zb(i5,i4)*
     -      (2*s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*
     -         za(i6,i7)*zb(i7,i2)**2*zb(i7,i5) - 
     -        2*s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i6)*
     -         zb(i5,i2)*zb(i7,i1)*zb(i7,i6) + 
     -        2*s(i2,i6)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*za(i6,i7)*
     -         zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) - 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)*za(i2,i6)*
     -         za(i6,i7)*zb(i5,i2)*zb(i7,i1)*zb(i7,i2)*zb(i7,i6) + 
     -        2*s(i1,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)**2*za(i6,i7)*
     -         zb(i5,i2)*zb(i7,i2)**2*zb(i7,i6)))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i7,i6)) - 
     -  (hb*s(i3,i4)*zb(i5,i4)*
     -     (s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*za(i3,i4)*
     -        za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i3) + 
     -       s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i4)*za(i2,i6)*
     -        za(i3,i4)*za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i4) - 
     -       s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i4)*
     -        za(i3,i6)*zb(i3,i2)*zb(i5,i4)*zb(i7,i1)*zb(i7,i6) - 
     -       s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i4)*
     -        za(i4,i6)*zb(i4,i2)*zb(i5,i4)*zb(i7,i1)*zb(i7,i6) + 
     -       s(i2,i6)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*za(i3,i4)*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i1)**2*zb(i7,i6) + 
     -       s(i2,i6)*t(i2,i6,i7)*za(i1,i4)*za(i1,i6)**2*za(i3,i4)*
     -        za(i6,i7)*zb(i4,i2)*zb(i5,i4)*zb(i7,i1)**2*zb(i7,i6) - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)*za(i2,i6)*
     -        za(i3,i4)*za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i1)*
     -        zb(i7,i2)*zb(i7,i6) - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i4)*za(i1,i6)*za(i2,i6)*
     -        za(i3,i4)*za(i6,i7)*zb(i4,i2)*zb(i5,i4)*zb(i7,i1)*
     -        zb(i7,i2)*zb(i7,i6) + 
     -       s(i1,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)**2*za(i3,i4)*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i6) + 
     -       s(i1,i7)*t(i1,i6,i7)*za(i1,i4)*za(i2,i6)**2*za(i3,i4)*
     -        za(i6,i7)*zb(i4,i2)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i6)))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_anomZa_mmp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21122
        complex(dp) :: zajj_tree_qqgg_anomZa_mmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomZa_mmp =
     -  (hb*s(i3,i4)*za(i3,i5)*
     -     (s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*za(i3,i5)*
     -        za(i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i2)**2*zb(i7,i3) + 
     -       s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i4)*za(i2,i6)*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i2)**2*
     -        zb(i7,i4) - s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*
     -        za(i3,i5)*(za(i3,i6)*zb(i3,i2) + za(i4,i6)*zb(i4,i2))*
     -        zb(i4,i3)*zb(i5,i3)*zb(i7,i1)*zb(i7,i6) + 
     -       s(i2,i6)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i5)*za(i6,i7)*
     -        (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3)*
     -        zb(i5,i3)*zb(i7,i1)**2*zb(i7,i6) - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i6)*za(i2,i6)*za(i3,i5)*
     -        za(i6,i7)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*
     -        zb(i4,i3)*zb(i5,i3)*zb(i7,i1)*zb(i7,i2)*zb(i7,i6) + 
     -       s(i1,i7)*t(i1,i6,i7)*za(i2,i6)**2*za(i3,i5)*za(i6,i7)*
     -        (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3)*
     -        zb(i5,i3)*zb(i7,i2)**2*zb(i7,i6)))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i5,i3)*zb(i7,i6)) + 
     -  (ha*s(i3,i4)*za(i3,i5)*
     -     (-2*s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*
     -        za(i6,i7)*zb(i7,i2)**2*
     -        (zb(i5,i4)*zb(i7,i3) + zb(i4,i3)*zb(i7,i5)) + 
     -       2*s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i5,i6)*
     -        (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*zb(i7,i1)*
     -        zb(i7,i6) - 2*s(i2,i6)*t(i2,i6,i7)*za(i1,i5)*za(i1,i6)**2*
     -        za(i6,i7)*(zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -        zb(i7,i1)**2*zb(i7,i6) + 
     -       2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i5)*za(i1,i6)*za(i2,i6)*
     -        za(i6,i7)*(zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -        zb(i7,i1)*zb(i7,i2)*zb(i7,i6) - 
     -       2*s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*za(i6,i7)*
     -        (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*zb(i7,i2)**2*
     -        zb(i7,i6)))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i5,i3)*zb(i7,i6))

      end function

c =========== 
c aZ anomalous coupling decay amplitudes
c =========== 
      pure function zajj_tree_qqQQ_anomaZ_plus(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! + photon
        complex(dp) :: zajj_tree_qqQQ_anomaZ_plus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_anomaZ_plus = -(
     -  (ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (-2*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*zb(i5,i2)*zb(i6,i1) + 
     -       2*t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) - 
     -       2*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*zb(i5,i2)*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)) + 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (s(i3,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*zb(i5,i2)*
     -        zb(i6,i1) + s(i4,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*
     -        zb(i5,i2)*zb(i6,i1) + 
     -       t(i2,i6,i7)*za(i1,i5)*za(i1,i7)*za(i3,i4)*zb(i5,i2)*
     -        zb(i5,i4)*zb(i6,i1) - 
     -       s(i3,i5)*t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) - 
     -       s(i4,i5)*t(i1,i6,i7)*za(i1,i3)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) - 
     -       t(i1,i6,i7)*za(i1,i5)*za(i3,i4)*zb(i5,i4)*zb(i6,i2)*
     -        (za(i2,i7)*zb(i5,i2) - za(i6,i7)*zb(i6,i5)) + 
     -       s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*zb(i5,i2)*
     -        zb(i7,i6) + s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i3,i7)*
     -        zb(i5,i2)*zb(i7,i6) + 
     -       t(i2,i6,i7)*za(i1,i7)*za(i3,i4)*za(i5,i7)*zb(i5,i2)*
     -        zb(i5,i4)*zb(i7,i6)))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)))

      end function

      pure function zajj_tree_qqQQ_anomaZ_minus(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! - photon
        complex(dp) :: zajj_tree_qqQQ_anomaZ_minus
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqQQ_anomaZ_minus =
     -  -(ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -      (-2*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i4,i3)*zb(i5,i2)*
     -         zb(i6,i2) - 2*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i3,i2)*
     -         zb(i5,i4)*zb(i6,i2) + 
     -        2*t(i1,i6,i7)*za(i1,i5)*za(i6,i7)*zb(i6,i2)*
     -         (zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)) + 
     -        2*t(i2,i6,i7)*za(i1,i7)*zb(i4,i3)*zb(i5,i2)*
     -         (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -        2*t(i2,i6,i7)*za(i1,i7)*zb(i3,i2)*zb(i5,i4)*
     -         (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3)) - 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i3,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i2) + 
     -       s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i2) + 
     -       t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*za(i4,i5)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i5,i4)*zb(i6,i2) - 
     -       s(i3,i5)*t(i1,i6,i7)*za(i1,i5)*za(i6,i7)*zb(i5,i4)*
     -        zb(i6,i2)*zb(i6,i3) - 
     -       s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i6,i7)*zb(i6,i2)*
     -        (zb(i5,i4)*zb(i6,i3) + zb(i4,i3)*zb(i6,i5)) - 
     -       s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -       t(i2,i6,i7)*za(i1,i7)*za(i4,i5)*zb(i4,i3)*zb(i5,i2)*
     -        zb(i5,i4)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*zb(i5,i3))

      end function

      pure function zajj_tree_qqgg_anomaZ_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22222
        complex(dp) :: zajj_tree_qqgg_anomaZ_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_ppp = -(
     -  (ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i1,i3)*zb(i5,i4)*
     -     (-(za(i1,i2)*zb(i5,i2)) + za(i1,i6)*zb(i6,i5) + 
     -       za(i1,i7)*zb(i7,i5)))/(2.*za(i1,i7)*za(i2,i6)*za(i6,i7)) - 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (s(i3,i5)*za(i1,i3) + s(i4,i5)*za(i1,i3) + 
     -       za(i1,i5)*za(i3,i4)*zb(i5,i4))*
     -     (-(za(i1,i2)*zb(i5,i2)) + za(i1,i6)*zb(i6,i5) + 
     -       za(i1,i7)*zb(i7,i5)))/(4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)))

      end function

      pure function zajj_tree_qqgg_anomaZ_mpp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21222
        complex(dp) :: zajj_tree_qqgg_anomaZ_mpp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_mpp =
     -  -(ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i1,i5)*za(i3,i5)*
     -      (2*za(i1,i2)*zb(i4,i3)*zb(i5,i2) + 
     -        2*za(i1,i2)*zb(i3,i2)*zb(i5,i4) - 
     -        2*za(i1,i6)*zb(i5,i4)*zb(i6,i3) - 
     -        2*za(i1,i6)*zb(i4,i3)*zb(i6,i5) - 
     -        2*za(i1,i7)*zb(i5,i4)*zb(i7,i3) - 
     -        2*za(i1,i7)*zb(i4,i3)*zb(i7,i5)))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3)) - 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i1,i5)*za(i3,i5)*
     -     (-(s(i4,i5)*za(i1,i2)*zb(i4,i3)*zb(i5,i2)) - 
     -       s(i3,i5)*za(i1,i2)*zb(i3,i2)*zb(i5,i4) - 
     -       s(i4,i5)*za(i1,i2)*zb(i3,i2)*zb(i5,i4) + 
     -       s(i3,i5)*za(i1,i6)*zb(i5,i4)*zb(i6,i3) + 
     -       s(i4,i5)*za(i1,i6)*zb(i5,i4)*zb(i6,i3) + 
     -       za(i1,i6)*za(i4,i5)*zb(i4,i3)*zb(i5,i4)*zb(i6,i5) + 
     -       s(i3,i5)*za(i1,i7)*zb(i5,i4)*zb(i7,i3) + 
     -       s(i4,i5)*za(i1,i7)*zb(i5,i4)*zb(i7,i3) + 
     -       za(i1,i7)*za(i4,i5)*zb(i4,i3)*zb(i5,i4)*zb(i7,i5)))/
     -   (4._dp*za(i1,i7)*za(i2,i6)*za(i6,i7)*zb(i5,i3))

      end function

      pure function zajj_tree_qqgg_anomaZ_pmm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22112
        complex(dp) :: zajj_tree_qqgg_anomaZ_pmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_pmm = -(
     -  (ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i2)*zb(i5,i4)*
     -     (-2*za(i1,i3)*zb(i2,i1) - 2*za(i3,i6)*zb(i6,i2) - 
     -       2*za(i3,i7)*zb(i7,i2)))/(4._dp*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))
     -   + (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i2)*zb(i5,i4)*
     -     (s(i3,i5)*za(i1,i3)*zb(i2,i1) + 
     -       s(i4,i5)*za(i1,i3)*zb(i2,i1) + 
     -       za(i1,i5)*za(i3,i4)*zb(i2,i1)*zb(i5,i4) + 
     -       s(i3,i5)*za(i3,i6)*zb(i6,i2) + 
     -       s(i4,i5)*za(i3,i6)*zb(i6,i2) + 
     -       za(i3,i4)*za(i5,i6)*zb(i5,i4)*zb(i6,i2) + 
     -       s(i3,i5)*za(i3,i7)*zb(i7,i2) + 
     -       s(i4,i5)*za(i3,i7)*zb(i7,i2) + 
     -       za(i3,i4)*za(i5,i7)*zb(i5,i4)*zb(i7,i2)))/
     -   (4._dp*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_anomaZ_mmm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21112
        complex(dp) :: zajj_tree_qqgg_anomaZ_mmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_mmm =
     -  (ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (-2*zb(i4,i3)*zb(i5,i2) - 2*zb(i3,i2)*zb(i5,i4))*
     -     (za(i1,i5)*zb(i2,i1) + za(i5,i6)*zb(i6,i2) + 
     -       za(i5,i7)*zb(i7,i2)))/
     -   (4._dp*zb(i5,i3)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6)) + 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i4,i5)*zb(i4,i3)*zb(i5,i2) + 
     -       s(i3,i5)*zb(i3,i2)*zb(i5,i4) + s(i4,i5)*zb(i3,i2)*zb(i5,i4)
     -       )*(za(i1,i5)*zb(i2,i1) + za(i5,i6)*zb(i6,i2) + 
     -       za(i5,i7)*zb(i7,i2)))/
     -   (4._dp*zb(i5,i3)*zb(i6,i2)*zb(i7,i1)*zb(i7,i6))

      end function

      pure function zajj_tree_qqgg_anomaZ_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22212
        complex(dp) :: zajj_tree_qqgg_anomaZ_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_ppm = -(
     -  -(ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -      (-2*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*zb(i5,i2)*
     -         zb(i6,i1)**2 - 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*zb(i6,i1)*
     -         (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5)) + 
     -        2*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i6,i2)*
     -         (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i1)
     -         - 2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i7)*zb(i5,i2)*
     -         zb(i6,i1)*zb(i7,i6) - 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*
     -         (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6))
     -      )/(4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i7,i1))
     -   - (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (s(i3,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*
     -        zb(i5,i2)*zb(i6,i1)**2 + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i7)*za(i2,i6)*
     -        zb(i5,i2)*zb(i6,i1)**2 + 
     -       t(i2,i6,i7)*za(i1,i5)*za(i1,i7)*za(i2,i6)*za(i3,i4)*
     -        zb(i5,i2)*zb(i5,i4)*zb(i6,i1)**2 + 
     -       s(i3,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*zb(i6,i1)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5)) + 
     -       s(i4,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*zb(i6,i1)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5)) + 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i5)*za(i3,i4)*zb(i5,i4)*
     -        zb(i6,i1)*(-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))
     -        - s(i3,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i1)
     -        - s(i4,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i7)*zb(i6,i2)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i1)
     -        - t(i1,i6,i7)*za(i1,i5)*za(i2,i7)*za(i3,i4)*zb(i5,i4)*
     -        zb(i6,i2)*(-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*
     -        zb(i7,i1) + s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*
     -        za(i3,i7)*zb(i5,i2)*zb(i6,i1)*zb(i7,i6) + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i7)*
     -        zb(i5,i2)*zb(i6,i1)*zb(i7,i6) + 
     -       t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*za(i3,i4)*za(i5,i7)*
     -        zb(i5,i2)*zb(i5,i4)*zb(i6,i1)*zb(i7,i6) + 
     -       s(i3,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6)
     -        + s(i4,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i7)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6)
     -        + t(i1,i6,i7)*t(i2,i6,i7)*za(i3,i4)*za(i5,i7)*zb(i5,i4)*
     -        (-(za(i2,i7)*zb(i5,i2)) + za(i6,i7)*zb(i6,i5))*zb(i7,i6)))
     -    /(4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i7,i1)))

      end function

      pure function zajj_tree_qqgg_anomaZ_mpm(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21212
        complex(dp) :: zajj_tree_qqgg_anomaZ_mpm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_mpm =
     -  -(ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -      (-2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i4,i3)*zb(i5,i2)*
     -         zb(i6,i1)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -        2*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i3,i2)*zb(i5,i4)*
     -         zb(i6,i1)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -        2*t(i1,i6,i7)*za(i2,i7)*zb(i4,i3)*zb(i5,i2)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           ) + 2*t(i1,i6,i7)*za(i2,i7)*zb(i3,i2)*zb(i5,i4)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           ) - 2*t(i1,i6,i7)*za(i6,i7)*zb(i5,i4)*zb(i6,i3)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           ) - 2*t(i1,i6,i7)*za(i6,i7)*zb(i4,i3)*zb(i6,i5)*
     -         (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -           t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))
     -           )))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i5,i3)*
     -     zb(i7,i1)) - (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -       s(i3,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) + 
     -       s(i4,i5)*t(i2,i6,i7)*za(i1,i7)*za(i2,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i6,i1)*
     -        (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)) - 
     -       s(i4,i5)*t(i1,i6,i7)*za(i2,i7)*zb(i4,i3)*zb(i5,i2)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         - s(i3,i5)*t(i1,i6,i7)*za(i2,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         - s(i4,i5)*t(i1,i6,i7)*za(i2,i7)*zb(i3,i2)*zb(i5,i4)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         + s(i3,i5)*t(i1,i6,i7)*za(i6,i7)*zb(i5,i4)*zb(i6,i3)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         + s(i4,i5)*t(i1,i6,i7)*za(i6,i7)*zb(i5,i4)*zb(i6,i3)*
     -        (-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -         + t(i1,i6,i7)*za(i4,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i4)*
     -        zb(i6,i5)*(-(za(i1,i5)*za(i2,i7)*zb(i6,i2)*zb(i7,i1)) + 
     -          t(i2,i6,i7)*(za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6)))
     -       ))/
     -   (4._dp*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*za(i2,i6)*zb(i5,i3)*
     -     zb(i7,i1))

      end function

      pure function zajj_tree_qqgg_anomaZ_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 22122
        complex(dp) :: zajj_tree_qqgg_anomaZ_pmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_pmp = -(
     -  -(ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -      (2*s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i6)*
     -         zb(i5,i2)*zb(i7,i1)*zb(i7,i6) - 
     -        2*s(i2,i6)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*za(i6,i7)*
     -         zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) + 
     -        2*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)*za(i2,i6)*
     -         za(i6,i7)*zb(i5,i2)*zb(i7,i1)*zb(i7,i2)*zb(i7,i6) - 
     -        2*s(i1,i7)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*za(i6,i7)*
     -         zb(i7,i2)**2*(s(i6,i7)*zb(i7,i5) + 
     -           za(i2,i6)*zb(i5,i2)*zb(i7,i6))))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i7,i6)) - 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*zb(i5,i4)*
     -     (-(s(i2,i6)*s(i3,i5)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*
     -          za(i3,i6)*zb(i5,i2)*zb(i7,i1)*zb(i7,i6)) - 
     -       s(i2,i6)*s(i4,i5)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*
     -        za(i3,i6)*zb(i5,i2)*zb(i7,i1)*zb(i7,i6) - 
     -       s(i2,i6)*s(i6,i7)*t(i2,i6,i7)*za(i1,i6)**2*za(i3,i4)*
     -        za(i5,i6)*zb(i5,i2)*zb(i5,i4)*zb(i7,i1)*zb(i7,i6) + 
     -       s(i2,i6)*s(i3,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*
     -        za(i6,i7)*zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) + 
     -       s(i2,i6)*s(i4,i5)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)**2*
     -        za(i6,i7)*zb(i5,i2)*zb(i7,i1)**2*zb(i7,i6) + 
     -       s(i2,i6)*t(i2,i6,i7)*za(i1,i5)*za(i1,i6)**2*za(i3,i4)*
     -        za(i6,i7)*zb(i5,i2)*zb(i5,i4)*zb(i7,i1)**2*zb(i7,i6) - 
     -       s(i3,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*za(i1,i6)*
     -        za(i2,i6)*za(i6,i7)*zb(i5,i2)*zb(i7,i1)*zb(i7,i2)*
     -        zb(i7,i6) - s(i4,i5)*t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i3)*
     -        za(i1,i6)*za(i2,i6)*za(i6,i7)*zb(i5,i2)*zb(i7,i1)*
     -        zb(i7,i2)*zb(i7,i6) - 
     -       t(i1,i6,i7)*t(i2,i6,i7)*za(i1,i5)*za(i1,i6)*za(i2,i6)*
     -        za(i3,i4)*za(i6,i7)*zb(i5,i2)*zb(i5,i4)*zb(i7,i1)*
     -        zb(i7,i2)*zb(i7,i6) + 
     -       s(i1,i7)*s(i3,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*
     -        za(i6,i7)*zb(i7,i2)**2*
     -        (s(i6,i7)*zb(i7,i5) + za(i2,i6)*zb(i5,i2)*zb(i7,i6)) + 
     -       s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i3)*za(i2,i6)*
     -        za(i6,i7)*zb(i7,i2)**2*
     -        (s(i6,i7)*zb(i7,i5) + za(i2,i6)*zb(i5,i2)*zb(i7,i6)) + 
     -       s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i3,i4)*
     -        za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*
     -        (s(i6,i7)*zb(i7,i5) + za(i2,i6)*zb(i5,i2)*zb(i7,i6))))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i7,i6)))

      end function

      pure function zajj_tree_qqgg_anomaZ_mmp(i1,i2,i3,i4,i5,i6,i7,za,zb,ha,hb) ! 21122
        complex(dp) :: zajj_tree_qqgg_anomaZ_mmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zajj_tree_qqgg_anomaZ_mmp =
     -  (ha*za(i3,i5)*(2*s(i1,i7)*s(i3,i4)*s(i6,i7)*t(i1,i6,i7)*
     -        za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i3) + 2*s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*
     -        za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i3) + 2*s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i3,i5)*za(i6,i7)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i3) + 
     -       2*s(i1,i7)*s(i3,i4)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i6,i7)*zb(i4,i3)*zb(i7,i2)**2*zb(i7,i5) + 
     -       2*s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i6,i7)*zb(i4,i3)*zb(i7,i2)**2*zb(i7,i5) + 
     -       2*s(i1,i7)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i2)**2*
     -        zb(i7,i5) + 2*s(i1,i7)*s(i3,i4)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)**2*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*zb(i7,i2)**2*
     -        zb(i7,i6) + 2*s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)**2*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*zb(i7,i2)**2*
     -        zb(i7,i6) + 2*s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*zb(i5,i3)*
     -        zb(i7,i2)**2*zb(i7,i6) + 
     -       2*s(i1,i7)*s(i3,i4)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i6) + 
     -       2*s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i2)**2*zb(i7,i6) + 
     -       2*s(i1,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*za(i3,i5)*
     -        za(i6,i7)*zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i6) - 2*s(i3,i4)*t(i2,i6,i7)*za(i1,i6)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i5,i3)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*s(i3,i4)*t(i2,i6,i7)*za(i1,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*zb(i3,i2)*
     -        zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - 2*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*zb(i3,i2)*
     -        zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6)))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i5,i3)*zb(i7,i6)) + 
     -  (hb*za(i3,i5)*(-(s(i1,i7)*s(i3,i4)*s(i4,i5)*s(i6,i7)*
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i5,i4)*
     -          zb(i7,i2)**2*zb(i7,i3)) - 
     -       s(i1,i7)*s(i3,i4)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*
     -        za(i3,i5)*za(i6,i7)*zb(i5,i3)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i3) - s(i1,i7)*s(i3,i5)*s(i6,i7)*t(i1,i6,i7)*
     -        za(i1,i5)*za(i2,i6)*za(i3,i5)*za(i6,i7)*zb(i5,i3)*
     -        zb(i5,i4)*zb(i7,i2)**2*zb(i7,i3) - 
     -       2*s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i3,i5)*za(i6,i7)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i3) - 
     -       s(i1,i7)*s(i3,i4)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i6,i7)*zb(i4,i3)*zb(i7,i2)**2*zb(i7,i5) - 
     -       s(i1,i7)*s(i4,i5)*s(i6,i7)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i2)**2*
     -        zb(i7,i5) - s(i1,i7)*s(i3,i4)*s(i4,i5)*t(i1,i6,i7)*
     -        za(i1,i5)*za(i2,i6)**2*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*
     -        zb(i7,i2)**2*zb(i7,i6) - 
     -       s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i4,i3)*zb(i5,i2)*zb(i5,i3)*
     -        zb(i7,i2)**2*zb(i7,i6) - 
     -       s(i1,i7)*s(i3,i4)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)**2*za(i6,i7)*zb(i3,i2)*zb(i5,i4)*zb(i7,i2)**2*
     -        zb(i7,i6) - s(i1,i7)*s(i3,i4)*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)**2*za(i3,i5)*za(i6,i7)*zb(i3,i2)*zb(i5,i3)*
     -        zb(i5,i4)*zb(i7,i2)**2*zb(i7,i6) - 
     -       s(i1,i7)*s(i3,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i6) - 
     -       2*s(i1,i7)*s(i4,i5)*t(i1,i6,i7)*za(i1,i5)*za(i2,i6)**2*
     -        za(i3,i5)*za(i6,i7)*zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*
     -        zb(i7,i2)**2*zb(i7,i6) + 
     -       s(i3,i4)*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*zb(i4,i3)*
     -        zb(i5,i2)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*
     -        zb(i4,i3)*zb(i5,i2)*zb(i5,i3)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i3,i4)*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*
     -        zb(i3,i2)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i3,i4)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*
     -        zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i3,i5)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*
     -        zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + 2*s(i4,i5)*t(i2,i6,i7)*za(i1,i6)*za(i3,i5)*
     -        zb(i3,i2)*zb(i5,i3)*zb(i5,i4)*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) + s(i4,i5)**2*t(i2,i6,i7)*za(i1,i6)*
     -        (zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*zb(i7,i1)*
     -        (s(i2,i6)*za(i1,i6)*
     -           (s(i6,i7)*za(i5,i6) - za(i1,i5)*za(i6,i7)*zb(i7,i1)) + 
     -          t(i1,i6,i7)*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i7,i2))*
     -        zb(i7,i6) - s(i1,i7)*s(i4,i5)**2*t(i1,i6,i7)*za(i1,i5)*
     -        za(i2,i6)*za(i6,i7)*zb(i7,i2)**2*
     -        (s(i6,i7)*(zb(i5,i4)*zb(i7,i3) + zb(i4,i3)*zb(i7,i5)) + 
     -          za(i2,i6)*(zb(i4,i3)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*
     -           zb(i7,i6))))/
     -   (4._dp*s(i1,i7)*s(i2,i6)*s(i6,i7)*t(i1,i6,i7)*t(i2,i6,i7)*
     -     za(i6,i7)*zb(i5,i3)*zb(i7,i6))

      end function

c ======================================
c === "cleaned up" qq-type qqgg amplitudes
c === before use warning:
c === - don't rely on the helicity labels (check)
c === - please check the subroutine xzqqagg_qq on how to use them
c ======================================

c     q+,qb-, g+,g+,ga+, lb-,l+
      pure function zajj_tree_qqgg_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb)
        complex(dp) :: zajj_tree_qqgg_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: ta,t2a,t2b,t22a,t31a
        real(dp) :: s,t

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        !<j1|j2|j3]
        ta(i1,i2,i3)=za(i1,i2)*zb(i2,i3)
        !<i1|i2+i3|i4]
        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        !<i1|i2|i3|i4>
        t2b(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i3)*za(i3,i4)
        !<i1|i2+i3|i4+i5|i6>
        t22a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i4,i6)+t2b(i1,i2,i5,i6)
     &   +t2b(i1,i3,i4,i6)+t2b(i1,i3,i5,i6)
        !<i1|i2+i3+i4|i5|i6>
        t31a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i5,i6)+t2b(i1,i3,i5,i6)
     &   +t2b(i1,i4,i5,i6)

        ! amplitude by R. Mondini
        ! checked against hep-ph/98006317 by symmetrizing over gluons
        zajj_tree_qqgg_ppp =
     &   -(za(i1,i2)*za(i2,i6)**2)/(za(i1,i4)*
     &   za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4)*za(i6,i7))

      end function

      pure function zajj_tree_qqgg_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb)
        complex(dp) :: zajj_tree_qqgg_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: ta,t2a,t2b,t22a,t31a
        real(dp) :: s,t

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        !<j1|j2|j3]
        ta(i1,i2,i3)=za(i1,i2)*zb(i2,i3)
        !<i1|i2+i3|i4]
        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        !<i1|i2|i3|i4>
        t2b(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i3)*za(i3,i4)
        !<i1|i2+i3|i4+i5|i6>
        t22a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i4,i6)+t2b(i1,i2,i5,i6)
     &   +t2b(i1,i3,i4,i6)+t2b(i1,i3,i5,i6)
        !<i1|i2+i3+i4|i5|i6>
        t31a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i5,i6)+t2b(i1,i3,i5,i6)
     &   +t2b(i1,i4,i5,i6)

        ! amplitude by R. Mondini
        ! checked against hep-ph/98006317 by symmetrizing over gluons
        zajj_tree_qqgg_ppm=
     &   +1/(za(i1,i4)*za(i6,i7))*((za(i2,i6)**2*t2a(i5,i1,i4,i3)**2)/
     &   (t(i1,i4,i5)*t(i2,i6,i7)*t22a(i4,i1,i5,i6,i7,i2))+(za(i2,i5)*t2a(i6,i2,i5,i3)**2)/
     &   (t(i2,i3,i5)*zb(i5,i2)*t31a(i4,i1,i6,i7,i5,i2)))+1/(t(i1,i6,i7)*za(i2,i3)*za(i3,i4)*
     &   t22a(i4,i1,i5,i6,i7,i2))*(-(t22a(i2,i3,i4,i1,i7,i6))/(t2a(i4,i2,i3,i5))*
     &   ((za(i2,i6)*t2a(i4,i2,i3,i1)*t2a(i2,i3,i4,i1))/(zb(i5,i1)*za(i6,i7))-
     &   (za(i2,i4)*(ta(i5,i1,i7)*za(i2,i4)-za(i4,i5)*t2a(i2,i3,i4,i7)))/(za(i1,i4)))
     &   -(zb(i7,i1)*za(i2,i5)*za(i2,i6)*t2a(i2,i3,i4,i1))/(zb(i5,i1))-(zb(i7,i1)*
     &   za(i1,i2)*za(i2,i5)*(ta(i5,i1,i7)*za(i2,i4)-za(i4,i5)*t2a(i2,i3,i4,i7)))/
     &   (zb(i7,i6)*za(i1,i4)))
      end function

      pure function zajj_tree_qqgg_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb)
        complex(dp) :: zajj_tree_qqgg_pmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: ta,t2a,t2b,t22a,t31a
        real(dp) :: s,t

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        !<j1|j2|j3]
        ta(i1,i2,i3)=za(i1,i2)*zb(i2,i3)
        !<i1|i2+i3|i4]
        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        !<i1|i2|i3|i4>
        t2b(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i3)*za(i3,i4)
        !<i1|i2+i3|i4+i5|i6>
        t22a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i4,i6)+t2b(i1,i2,i5,i6)
     &   +t2b(i1,i3,i4,i6)+t2b(i1,i3,i5,i6)
        !<i1|i2+i3+i4|i5|i6>
        t31a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i5,i6)+t2b(i1,i3,i5,i6)
     &   +t2b(i1,i4,i5,i6)

        ! amplitude by R. Mondini
        ! checked against hep-ph/98006317 by symmetrizing over gluons
        zajj_tree_qqgg_pmp =
     &   za(i2,i6)**2*((t2a(i4,i1,i5,i3)**3)/(t(i1,i4,i5)*t(i2,i6,i7)*za(i1,i5)*
     &   t22a(i4,i1,i5,i6,i7,i2))-((zb(i3,i1))**3*t2a(i2,i1,i4,i3))/(t(i1,i3,i4)*zb(i4,i1)*
     &   zb(i4,i3)*za(i2,i5)*t2a(i2,i3,i4,i1)))/(za(i6,i7)*t2a(i5,i1,i4,i3))+
     &   (za(i2,i4)**3*((t2a(i6,i1,i7,i5)*(t2a(i4,i1,i5,i7)/za(i1,i5)-(zb(i5,i1)*
     &   za(i2,i6)*t2a(i4,i2,i3,i1))/(za(i6,i7)*t2a(i2,i3,i4,i1))))/t(i2,i3,i4)-(zb(i7,i1)*
     &   ((za(i1,i2)*t2a(i4,i1,i5,i7))/(zb(i7,i6)*za(i1,i5))+(zb(i5,i1)*za(i2,i4)*
     &   za(i2,i6))/t2a(i2,i3,i4,i1)))/za(i2,i5)))/(t(i1,i6,i7)*za(i2,i3)*za(i3,i4)*
     &   t22a(i4,i1,i5,i6,i7,i2))

      end function

      pure function zajj_tree_qqgg_mpp(i1,i2,i3,i4,i5,i6,i7,za,zb)
        complex(dp) :: zajj_tree_qqgg_mpp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: ta,t2a,t2b,t22a,t31a
        real(dp) :: s,t

        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        !<j1|j2|j3]
        ta(i1,i2,i3)=za(i1,i2)*zb(i2,i3)
        !<i1|i2+i3|i4]
        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        !<i1|i2|i3|i4>
        t2b(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i3)*za(i3,i4)
        !<i1|i2+i3|i4+i5|i6>
        t22a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i4,i6)+t2b(i1,i2,i5,i6)
     &   +t2b(i1,i3,i4,i6)+t2b(i1,i3,i5,i6)
        !<i1|i2+i3+i4|i5|i6>
        t31a(i1,i2,i3,i4,i5,i6)=t2b(i1,i2,i5,i6)+t2b(i1,i3,i5,i6)
     &   +t2b(i1,i4,i5,i6)

        ! amplitude by R. Mondini
        ! checked against hep-ph/98006317 by symmetrizing over gluons
        zajj_tree_qqgg_mpp =
     &   ((zb(i4,i2)*t2a(i1,i2,i3,i4)*t2a(i6,i2,i3,i4)**2)/(t(i2,i3,i4)*zb(i3,i2)*zb(i4,i3)*
     &   za(i1,i5)*t2a(i1,i3,i4,i2))+(za(i2,i3)**2*t2a(i6,i1,i7,i4)**2*t2a(i3,i2,i5,i4))/
     &   (t(i1,i6,i7)*t(i2,i3,i5)*za(i2,i5)*t22a(i3,i2,i5,i6,i7,i1)))/(za(i6,i7)*t2a(i5,i2,i3,i4))+
     &   (za(i1,i3)*((za(i2,i3)*t2a(i3,i1,i4,i7)*((za(i1,i2)*t2a(i3,i2,i6,i7))/
     &   (zb(i7,i6)*za(i1,i5))-(za(i2,i6)*t2a(i3,i1,i4,i5))/t(i1,i3,i4)))/za(i2,i5)+
     &   (1/t2a(i1,i3,i4,i2))*(za(i1,i3)*t2a(i6,i2,i3,i5)-ta(i1,i4,i5)*za(i3,i6))*
     &   ((za(i1,i3)*t2a(i3,i2,i6,i7))/za(i1,i5)-(za(i2,i6)*t2a(i3,i1,i4,i2)*
     &   t2a(i3,i1,i4,i5))/(t(i1,i3,i4)*za(i6,i7)))))/(t(i2,i6,i7)*za(i1,i4)*za(i3,i4)*
     &   t22a(i3,i2,i5,i6,i7,i1))

      end function


      end module
