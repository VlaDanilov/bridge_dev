      w(1)=Log(1 + xu*xv**(-1))
      w(2)=Log(xu)
      w(3)=Log(xv)
      w(4)=Log(1 - xv)
      w(5)=Li2( - xu*xv**(-1))
      w(6)=Li2(xv)
      w(7)=OneMu**(-1)
      w(8)=xu**(-1)
      w(9)=Log(1 - xu)
      w(10)=Log(1 - 1/(1 - xv)*xu)
      w(11)=Li2(1/(1 - xu)*xv)
      w(12)=OneMv**(-1)
      w(13)=w(3) + w(1)
      w(13)=w(13)*w(2)
      w(14)=Pi*im
      w(15)=w(14)*w(1)
      w(13)=w(13) + zeta2 - w(6) + w(15) + w(5)
      w(15)=w(1) + w(4) - w(14) + 1._dp/2._dp*w(3)
      w(16)= - w(3)*w(15)
      w(16)=w(16) + w(13)
      w(16)=CF*w(16)
      w(17)=w(14)*w(9)
      w(18)=w(10) - w(9)
      w(19)=w(18)*w(3)
      w(17)=w(17) + w(19)
      w(20)=w(11) - w(6)
      w(21)=w(9)**2
      w(22)=1._dp/2._dp*w(21) + w(20)
      w(23)=w(22) + w(17)
      w(23)=w(23)*w(8)
      w(24)=w(10) + w(4)
      w(25)=w(24) + w(14)
      w(26)=w(25) - w(3)
      w(27)=w(23) + w(26)
      w(27)=w(27)*w(8)
      w(28)=1._dp/2._dp*w(7)
      w(29)=w(28) - 1
      w(29)=w(7)*w(29)*w(26)
      w(27)=w(27) - w(29)
      w(27)=w(27)*xv
      w(17)=w(20) + w(17)
      w(17)= - w(21) - 2*w(17)
      w(17)=w(8)*w(17)
      w(17)=w(17) - 2*w(25) + w(3)
      w(17)=w(8)*w(17)
      w(20)= - 1._dp/2._dp - w(26)
      w(20)=w(7)*w(20)
      w(17)=w(27) + w(20) + w(17)
      w(17)=xv*w(17)
      w(20)=w(23) + w(25)
      w(20)=w(20)*w(8)
      w(21)= - 7 + w(25)
      w(17)=w(17) + 1._dp/2._dp*w(21) + w(20)
      w(17)=CF*w(17)
      w(21)=xu*w(12)
      w(29)=w(21) - 1
      w(30)=4*CF
      w(31)=w(29)*w(30)
      w(32)=w(12) - 1
      w(32)=w(32)*w(3)
      w(33)=w(32) + 1
      w(34)=w(33)*w(21)
      w(15)= - 3._dp/2._dp*w(12) + 3._dp/2._dp - w(15)
      w(15)=w(3)*w(15)
      w(13)=1._dp/2._dp*w(34) + w(15) - 9._dp/2._dp + w(13)
      w(13)=CF*w(13)
      w(15)=w(9) - 1
      w(14)=w(15)*w(14)
      w(14)=w(14) - w(23) + w(22) - w(24)
      w(15)=1 + w(18)
      w(15)=w(3)*w(15)
      w(15)=w(15) + w(14)
      w(15)=w(8)*w(15)
      w(18)=w(26)*w(28)
      w(15)=w(18) + w(15)
      w(15)=xv*w(15)
      w(18)=w(34) + 3*w(32) + 1 + w(25)
      w(15)=w(15) + 1._dp/2._dp*w(18) + w(20)
      w(15)=CF*w(15)
      w(18)=w(21)*w(30)
      w(20)= - w(34) + w(33)
      w(20)=1._dp/2._dp*CF*w(20)
      w(14)= - w(19) - w(14)
      w(14)=w(8)*w(14)
      w(19)=1 - w(26)
      w(19)=w(19)*w(28)
      w(14)= - w(27) + w(19) + w(14)
      w(14)=xv*w(14)
      w(19)= - w(32) - w(34)
      w(14)=1._dp/2._dp*w(19) + w(14)
      w(14)=CF*w(14)


      alphaG1WL2= 0

      alphaG2WL2= 0

      alphaG3WL2= 0

      alphaG4WL2= 0

      betaG1WL2= 0

      betaG2WL2= 0

      betaG3WL2= 0

      betaG4WL2= 0

      gammaG1WL2= 0

      gammaG2WL2= 0

      gammaG3WL2= 0

      gammaG4WL2= 0

      alphaG1WL1=w(16)

      alphaG2WL1=w(17)

      alphaG3WL1=w(31)

      betaG1WL1=w(13)

      betaG2WL1=w(15)

      betaG3WL1=w(18)

      gammaG1WL1=w(20)

      gammaG2WL1=w(14)

      gammaG3WL1= 0

      alphaG1WL0= 0

      alphaG2WL0=1

      alphaG3WL0= - w(29)

      betaG1WL0=1

      betaG2WL0= 0

      betaG3WL0= - w(21)

      gammaG1WL0= 0

      gammaG2WL0= 0

      gammaG3WL0= 0
