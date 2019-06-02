      function Fcc_qpqmgmgp(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Fcc_qpqmgmgp
      
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: L0,Lsm1,Lsm1_2mh,Lsm1_2mht,I3m
      real(dp):: t  

      Fcc_qpqmgmgp=
     .((Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))+
     .Lsm1_2mht(s(j1,j4),t(j1,j2,j3),s(j2,j3),s(j5,j6)))*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))**2)/
     .(za(j5,j6)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j2,j3)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(2._dp*L0(-t(j2,j3,j4),-s(j5,j6))*za(j1,j5)*zb(j1,j4)*zb(j2,j4)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))/
     .(s(j5,j6)*za(j5,j6)*zb(j2,j3)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*zb(j3,j4))-
     .(I3m(s(j1,j4),s(j2,j3),s(j5,j6))*zb(j1,j4)*
     .(-(s(j1,j4)*s(j2,j3)*za(j2,j5)**2)+
     .(-(za(j1,j2)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3)))+
     .za(j2,j4)*(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))**2))/
     .(2._dp*za(j1,j4)*za(j5,j6)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j
     .4))*
     .zb(j2,j3)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))-
     .((Lsm1(-s(j1,j2),-t(j1,j2,j4),-s(j1,j4),-t(j1,j2,j4))+
     .Lsm1_2mht(s(j2,j3),t(j1,j2,j4),s(j1,j4),s(j5,j6)))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6))**2)/
     .(za(j1,j4)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j5,j6))-
     .(2._dp*L0(-t(j1,j3,j4),-s(j5,j6))*za(j1,j3)*za(j2,j3)*zb(j2,j6)*
     .(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     .(s(j5,j6)*za(j1,j4)*za(j3,j4)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*zb(j5,j6))+
     .(I3m(s(j2,j3),s(j1,j4),s(j5,j6))*za(j2,j3)*
     .(-(s(j1,j4)*s(j2,j3)*zb(j1,j6)**2)+
     .(zb(j1,j2)*(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6))+
     .zb(j1,j3)*(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))**2))/
     .(2._dp*za(j1,j4)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j2,j3)*zb(j5,j6))+
     .Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*
     .(-(((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j1,j5)*zb(j1,j3))-za(j2,j5)*zb(j2,j3))**2)/
     .(za(j5,j6)*zb(j2,j3)*(-(za(j1,j4)*zb(j1,j3))-
     .za(j2,j4)*zb(j2,j3))**3))+
     .(za(j4,j5)**2*zb(j1,j3)**2*t(j1,j2,j3)**2)/
     .(za(j5,j6)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j2,j3)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**3))
     
      Fcc_qpqmgmgp=Fcc_qpqmgmgp+
     .Lsm1_2mh(s(j3,j4),t(j1,j2,j4),s(j1,j2),s(j5,j6))*
     .(((-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6))**2*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))/
     .(za(j1,j4)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**3*
     .zb(j5,j6))-(za(j2,j4)**2*zb(j3,j6)**2*t(j1,j2,j4)**2)/
     .(za(j1,j4)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**3*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j5,j6)))+
     .(Lsm1(-s(j1,j4),-t(j1,j3,j4),-s(j3,j4),-t(j1,j3,j4))*za(j2,j5)**2*
     .zb(j1,j4)**2)/
     .(za(j5,j6)*zb(j1,j3)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .t(j1,j3,j4))-(2._dp*L0(-t(j1,j3,j4),-s(j3,j4))*za(j1,j3)*zb(j1,j4)*
     .(za(j1,j3)*zb(j3,j6)+za(j1,j4)*zb(j4,j6))*
     .(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     .(s(j3,j4)*za(j1,j4)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .zb(j5,j6)*t(j1,j3,j4))+Lsm1_2mh(s(j1,j2),t(j1,j3,j4),s(j3,j4),
     .s(j5,j6))*(-((za(j2,j5)**2*zb(j1,j4)**2*
     .(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4)))/
     .(za(j5,j6)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .zb(j3,j4)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .t(j1,j3,j4)))-(za(j1,j3)*
     .(-(za(j1,j3)*zb(j1,j2))-za(j3,j4)*zb(j2,j4))**2*
     .(za(j1,j3)*zb(j3,j6)+za(j1,j4)*zb(j4,j6))**2)/
     .(za(j1,j4)*za(j3,j4)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .*
     .3_dp*zb(j5,j6)*t(j1,j3,j4))+
     .(za(j1,j3)**3*zb(j2,j6)**2*t(j1,j3,j4))/
     .(za(j1,j4)*za(j3,j4)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .*3*
     .zb(j5,j6)))
     
      Fcc_qpqmgmgp=Fcc_qpqmgmgp+
     .I3m(s(j1,j2),s(j3,j4),s(j5,j6))*
     .((za(j1,j2)*za(j3,j5)*zb(j1,j4)*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))/
     .(za(j1,j4)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*t(j1,j2,j4
     .))
     .+(4._dp*za(j1,j2)*za(j2,j4)*za(j3,j5)*zb(j1,j4)*zb(j2,j6)*zb(j3,j4)-
     .s(j2,j3)*za(j1,j5)*zb(j1,j6)*
     .(-3._dp*za(j1,j2)*zb(j1,j4)+za(j2,j3)*zb(j3,j4))-
     .s(j3,j4)*za(j1,j2)*za(j3,j5)*zb(j1,j4)*zb(j3,j6)+
     .s(j2,j4)*za(j2,j3)*(-(za(j1,j5)*zb(j1,j4))-za(j2,j5)*zb(j2,j4))*
     .zb(j3,j6)+(za(j1,j5)*za(j2,j4)*zb(j1,j4)-
     .za(j2,j4)*za(j2,j5)*zb(j2,j4))*
     .(-(za(j1,j3)*zb(j1,j4)*zb(j3,j6))-
     .za(j2,j3)*zb(j2,j4)*zb(j3,j6))-
     .(za(j1,j2)*za(j3,j4)*zb(j1,j4)*zb(j2,j3)*
     .(-2._dp*za(j1,j2)*za(j3,j5)*zb(j3,j6)+
     .za(j2,j5)*(za(j1,j2)*zb(j2,j6)-za(j1,j3)*zb(j3,j6))))/
     .za(j1,j4)-za(j2,j4)*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j4))-za(j2,j3)*zb(j3,j4))*zb(j4,j6)-
     .za(j4,j5)*(-(s(j2,j4)*za(j1,j2)*zb(j1,j4))+
     .s(j1,j2)*za(j2,j3)*zb(j3,j4))*zb(j4,j6)-
     .s(j1,j4)*za(j3,j5)*zb(j3,j4)*
     .(za(j2,j3)*zb(j3,j6)+za(j2,j4)*zb(j4,j6))+
     .(-(za(j1,j3)*za(j2,j3)*zb(j1,j3))-
     .za(j1,j3)*za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j5)*zb(j1,j4)*zb(j3,j6))-
     .za(j2,j5)*zb(j2,j3)*zb(j4,j6))-
     .zb(j1,j4)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .(za(j1,j3)*za(j2,j5)*zb(j2,j6)+
     .2_dp*za(j1,j5)*(-(za(j1,j3)*zb(j1,j6))-za(j3,j4)*zb(j4,j6)))+
     .2_dp*za(j2,j3)*za(j2,j5)*zb(j3,j4)*
     .(-(za(j1,j4)*zb(j1,j6)*zb(j2,j4))-
     .zb(j2,j6)*(s(j1,j2)+t(j1,j2,j4))))/
     .(2._dp*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))+
     .(za(j2,j5)*zb(j1,j4)*(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))*
     .(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     .((-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*t(j1,j3,j4)**2)-
     .(2._dp*za(j1,j3)**2*za(j2,j5)*zb(j1,j4)*zb(j2,j6))/
     .(za(j1,j4)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*t(j1,j3,j4
     .))
     .+(za(j1,j3)*za(j2,j3)*za(j2,j5)*zb(j1,j4)*zb(j1,j6))/
     .(2._dp*za(j3,j4)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .t(j1,j3,j4))+(za(j2,j5)*zb(j1,j4)*
     .(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))/
     .(2._dp*za(j5,j6)*zb(j3,j4)*(-(za(j1,j2)*zb(j1,j3))-
     .za(j2,j4)*zb(j3,j4))*t(j1,j3,j4))+
     .(za(j1,j3)*(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6))**2*
     .(2._dp*s(j3,j4)*s(j5,j6)+(s(j1,j2)-s(j3,j4)-s(j5,j6))*t(j1,j3,j4)))/
     .(2._dp*za(j1,j4)*za(j3,j4)*(-(za(j1,j3)*zb(j2,j3))-
     .za(j1,j4)*zb(j2,j4))*zb(j5,j6)*t(j1,j3,j4)**2)+
     .(-((za(j1,j3)*za(j2,j3)*za(j2,j5)*zb(j1,j6)*zb(j2,j4))/za(j3,j4))+
     .za(j2,j3)*(-(za(j1,j5)*zb(j1,j4))-2._dp*za(j3,j5)*zb(j3,j4))*
     .zb(j4,j6)+(2._dp*za(j1,j2)*za(j5,j6)*zb(j4,j6)*
     .(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/za(j1,j4)-
     .(za(j2,j5)*zb(j1,j4)*
     .(-(za(j1,j5)*zb(j1,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))-
     .za(j1,j2)*za(j5,j6)*zb(j2,j4)*zb(j4,j6)))/
     .(za(j5,j6)*zb(j3,j4))+
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j3)*za(j2,j3)*zb(j3,j6)*
     .(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     .(za(j1,j4)*za(j3,j4)*zb(j5,j6))-
     .(za(j1,j3)*za(j2,j5)*(-(za(j1,j3)*zb(j1,j6)*
     .(s(j5,j6)-t(j1,j2,j4)))+
     .za(j3,j4)*zb(j4,j6)*(s(j5,j6)+t(j1,j3,j4))))/
     .(za(j1,j4)*za(j3,j4)))/
     .(2._dp*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))))
     
      Fcc_qpqmgmgp=Fcc_qpqmgmgp+
     .(2._dp*L0(-t(j2,j3,j4),-s(j3,j4))*za(j2,j3)*zb(j2,j4)*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))/
     .(s(j3,j4)*za(j5,j6)*zb(j2,j3)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*t(j2,j3,j4))-
     .(Lsm1(-s(j2,j3),-t(j2,j3,j4),-s(j3,j4),-t(j2,j3,j4))*za(j2,j3)**2*
     .zb(j1,j6)**2)/
     .(za(j2,j4)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j5,j6)*
     .t(j2,j3,j4))+Lsm1_2mh(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))*
     .(-((zb(j2,j4)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))**2*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))**2)/
     .(za(j5,j6)*zb(j2,j3)*(-(za(j1,j3)*zb(j2,j3))-
     .za(j1,j4)*zb(j2,j4))**3*zb(j3,j4)*t(j2,j3,j4)))-
     .(za(j2,j3)**2*(za(j2,j3)*zb(j1,j2)-za(j3,j4)*zb(j1,j4))*zb(j1,j6)*
     .*2)/
     .(za(j3,j4)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j5,j6)*t(j2,j3,j4)
     .)+(za(j1,j5)**2*zb(j2,j4)**3*t(j2,j3,j4))/
     .(za(j5,j6)*zb(j2,j3)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .*3*
     .zb(j3,j4)))
     
      Fcc_qpqmgmgp=Fcc_qpqmgmgp+
     .I3m(s(j1,j2),s(j3,j4),s(j5,j6))*
     .(-((za(j2,j3)*zb(j1,j2)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .zb(j4,j6))/
     .((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j2,j3)*t(j1,j2,j3)))
     .+(-(s(j1,j4)*za(j2,j5)*(3._dp*za(j2,j3)*zb(j1,j2)-za(j3,j4)*zb(j1,j4
     .))*
     .zb(j2,j6))+s(j1,j3)*za(j4,j5)*zb(j1,j4)*
     .(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))+
     .(-(za(j2,j3)*zb(j1,j3)*zb(j2,j4))-
     .za(j2,j4)*zb(j1,j4)*zb(j2,j4))*
     .(-(za(j1,j4)*za(j3,j5)*zb(j1,j6))-
     .za(j2,j3)*za(j4,j5)*zb(j2,j6))+
     .(-(za(j1,j3)*za(j4,j5)*zb(j1,j4))-
     .za(j2,j3)*za(j4,j5)*zb(j2,j4))*
     .(-(za(j1,j3)*zb(j1,j3)*zb(j1,j6))+
     .za(j2,j3)*zb(j1,j3)*zb(j2,j6))-
     .za(j2,j3)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(za(j1,j5)*zb(j1,j6)*zb(j2,j4)+
     .2_dp*zb(j2,j6)*(-(za(j2,j5)*zb(j2,j4))+za(j3,j5)*zb(j3,j4)))-
     .za(j3,j5)*(s(j1,j3)*za(j2,j3)*zb(j1,j2)-
     .s(j1,j2)*za(j3,j4)*zb(j1,j4))*zb(j3,j6)+
     .s(j3,j4)*za(j2,j3)*za(j4,j5)*zb(j1,j2)*zb(j4,j6)+
     .4_dp*za(j1,j5)*za(j2,j3)*za(j3,j4)*zb(j1,j2)*zb(j1,j3)*zb(j4,j6)+
     .s(j2,j3)*za(j3,j4)*(za(j3,j5)*zb(j1,j3)+za(j4,j5)*zb(j1,j4))*
     .zb(j4,j6)-za(j3,j5)*zb(j1,j3)*
     .(za(j2,j3)*zb(j1,j2)+za(j3,j4)*zb(j1,j4))*
     .(za(j1,j3)*zb(j3,j6)+za(j1,j4)*zb(j4,j6))-
     .(za(j1,j4)*za(j2,j3)*zb(j1,j2)*zb(j3,j4)*
     .(zb(j1,j6)*(-(za(j1,j5)*zb(j1,j2))-za(j4,j5)*zb(j2,j4))+
     .2_dp*za(j4,j5)*zb(j1,j2)*zb(j4,j6)))/zb(j2,j3)-
     .2_dp*za(j3,j4)*zb(j1,j4)*zb(j1,j6)*
     .(-(za(j1,j3)*za(j2,j5)*zb(j2,j3))-
     .za(j1,j5)*(s(j1,j2)+t(j1,j2,j3))))/
     .(2._dp*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4)))+
     .(za(j2,j3)*(za(j2,j3)*zb(j1,j2)-za(j3,j4)*zb(j1,j4))*zb(j1,j6)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))/
     .((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*t(j2,j3,j4)**2)-
     .(2._dp*za(j1,j5)*za(j2,j3)*zb(j1,j6)*zb(j2,j4)**2)/
     .(zb(j2,j3)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*t(j2,j3,j4
     .))
     .-(za(j2,j3)*za(j2,j5)*zb(j1,j4)*zb(j1,j6)*zb(j2,j4))/
     .(2._dp*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j3,j4)*t(j2,j3,j
     .4))+
     .(za(j2,j3)*(za(j2,j3)*zb(j1,j2)-za(j3,j4)*zb(j1,j4))*zb(j1,j6)*
     .(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     .(2._dp*za(j3,j4)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j5,j6)
     .*
     .t(j2,j3,j4))+(zb(j2,j4)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))**2*
     .(2._dp*s(j3,j4)*s(j5,j6)+(s(j1,j2)-s(j3,j4)-s(j5,j6))*t(j2,j3,j4)))/
     .(2._dp*za(j5,j6)*zb(j2,j3)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j
     .4))*
     .zb(j3,j4)*t(j2,j3,j4)**2)+
     .((za(j1,j3)*za(j2,j5)*zb(j1,j4)*zb(j1,j6)*zb(j2,j4))/zb(j3,j4)+
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j4,j5)*zb(j1,j4)*zb(j2,j4)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))/
     .(za(j5,j6)*zb(j2,j3)*zb(j3,j4))+
     .za(j3,j5)*zb(j1,j4)*(-(za(j2,j3)*zb(j2,j6))+
     .2_dp*za(j3,j4)*zb(j4,j6))+
     .(2._dp*za(j3,j5)*zb(j1,j2)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))*zb(j5,j6))/
     .zb(j2,j3)-(za(j2,j3)*zb(j1,j6)*
     .(-(za(j2,j3)*(-(za(j1,j3)*zb(j1,j2))-za(j3,j4)*zb(j2,j4))*
     .zb(j2,j6))-za(j1,j3)*za(j3,j5)*zb(j1,j2)*zb(j5,j6)))/
     .(za(j3,j4)*zb(j5,j6))+
     .(zb(j1,j6)*zb(j2,j4)*(-(za(j2,j5)*zb(j2,j4)*
     .(s(j5,j6)-t(j1,j2,j3)))-
     .za(j3,j5)*zb(j3,j4)*(s(j5,j6)+t(j2,j3,j4))))/
     .(zb(j2,j3)*zb(j3,j4)))/
     .(2._dp*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))))

      return
      end
      