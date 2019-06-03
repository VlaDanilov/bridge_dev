      function pvSPPKKL(n1,n2,n3,n4,n5,P,K,L)
      implicit none
      include 'types.f'
      real(dp):: pvSPPKKL
      integer:: n1,n2,n3,n4,n5
      real(dp):: P(4),K(4),L(4)
      pvSPPKKL=
     . + P(n1)*P(n2)*K(n3)*K(n4)*L(n5)
     . + P(n1)*P(n2)*K(n3)*K(n5)*L(n4)
     . + P(n1)*P(n2)*K(n4)*K(n5)*L(n3)
     . + P(n1)*P(n3)*K(n2)*K(n4)*L(n5)
     . + P(n1)*P(n3)*K(n2)*K(n5)*L(n4)
     . + P(n1)*P(n3)*K(n4)*K(n5)*L(n2)
     . + P(n1)*P(n4)*K(n2)*K(n3)*L(n5)
     . + P(n1)*P(n4)*K(n2)*K(n5)*L(n3)
     . + P(n1)*P(n4)*K(n3)*K(n5)*L(n2)
     . + P(n1)*P(n5)*K(n2)*K(n3)*L(n4)
     . + P(n1)*P(n5)*K(n2)*K(n4)*L(n3)
     . + P(n1)*P(n5)*K(n3)*K(n4)*L(n2)
      pvSPPKKL=pvSPPKKL
     . + P(n2)*P(n3)*K(n1)*K(n4)*L(n5)
     . + P(n2)*P(n3)*K(n1)*K(n5)*L(n4)
     . + P(n2)*P(n3)*K(n4)*K(n5)*L(n1)
     . + P(n2)*P(n4)*K(n1)*K(n3)*L(n5)
     . + P(n2)*P(n4)*K(n1)*K(n5)*L(n3)
     . + P(n2)*P(n4)*K(n3)*K(n5)*L(n1)
     . + P(n2)*P(n5)*K(n1)*K(n3)*L(n4)
     . + P(n2)*P(n5)*K(n1)*K(n4)*L(n3)
     . + P(n2)*P(n5)*K(n3)*K(n4)*L(n1)
     . + P(n3)*P(n4)*K(n1)*K(n2)*L(n5)
     . + P(n3)*P(n4)*K(n1)*K(n5)*L(n2)
     . + P(n3)*P(n4)*K(n2)*K(n5)*L(n1)
     . + P(n3)*P(n5)*K(n1)*K(n2)*L(n4)
     . + P(n3)*P(n5)*K(n1)*K(n4)*L(n2)
     . + P(n3)*P(n5)*K(n2)*K(n4)*L(n1)
     . + P(n4)*P(n5)*K(n1)*K(n2)*L(n3)
     . + P(n4)*P(n5)*K(n1)*K(n3)*L(n2)
     . + P(n4)*P(n5)*K(n2)*K(n3)*L(n1)
      return
      end

      function pvSPPPKL(n1,n2,n3,n4,n5,P,K,L)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPKL
      integer n1,n2,n3,n4,n5
      real(dp):: P(4),K(4),L(4)
      pvSPPPKL=
     . + P(n1)*P(n2)*P(n3)*K(n4)*L(n5)
     . + P(n1)*P(n2)*P(n3)*K(n5)*L(n4)
     . + P(n1)*P(n2)*P(n4)*K(n3)*L(n5)
     . + P(n1)*P(n2)*P(n4)*K(n5)*L(n3)
     . + P(n1)*P(n2)*P(n5)*K(n3)*L(n4)
     . + P(n1)*P(n2)*P(n5)*K(n4)*L(n3)
     . + P(n1)*P(n3)*P(n4)*K(n2)*L(n5)
     . + P(n1)*P(n3)*P(n4)*K(n5)*L(n2)
     . + P(n1)*P(n3)*P(n5)*K(n2)*L(n4)
     . + P(n1)*P(n3)*P(n5)*K(n4)*L(n2)
     . + P(n1)*P(n4)*P(n5)*K(n2)*L(n3)
     . + P(n1)*P(n4)*P(n5)*K(n3)*L(n2)
     . + P(n2)*P(n3)*P(n4)*K(n1)*L(n5)
     . + P(n2)*P(n3)*P(n4)*K(n5)*L(n1)
     . + P(n2)*P(n3)*P(n5)*K(n1)*L(n4)
     . + P(n2)*P(n3)*P(n5)*K(n4)*L(n1)
     . + P(n2)*P(n4)*P(n5)*K(n1)*L(n3)
     . + P(n2)*P(n4)*P(n5)*K(n3)*L(n1)
     . + P(n3)*P(n4)*P(n5)*K(n1)*L(n2)
     . + P(n3)*P(n4)*P(n5)*K(n2)*L(n1)
      return
      end

      function pvSPPPKK(n1,n2,n3,n4,n5,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPKK
      integer n1,n2,n3,n4,n5
      real(dp):: P(4),K(4)
      pvSPPPKK=
     . + P(n1)*P(n2)*P(n3)*K(n4)*K(n5)
     . + P(n1)*P(n2)*P(n4)*K(n3)*K(n5)
     . + P(n1)*P(n2)*P(n5)*K(n3)*K(n4)
     . + P(n1)*P(n3)*P(n4)*K(n2)*K(n5)
     . + P(n1)*P(n3)*P(n5)*K(n2)*K(n4)
     . + P(n1)*P(n4)*P(n5)*K(n2)*K(n3)
     . + P(n2)*P(n3)*P(n4)*K(n1)*K(n5)
     . + P(n2)*P(n3)*P(n5)*K(n1)*K(n4)
     . + P(n2)*P(n4)*P(n5)*K(n1)*K(n3)
     . + P(n3)*P(n4)*P(n5)*K(n1)*K(n2)
      return
      end

      function pvSDDPKL(n1,n2,n3,n4,n5,P,K,L)
      implicit none 
      include 'types.f'
      real(dp):: pvSDDPKL
      integer n1,n2,n3,n4,n5
      real(dp):: P(4),K(4),L(4)
      include 'TRmetric.f'
      pvSDDPKL=
     . + g(n1,n2)*P(n3)*K(n4)*L(n5)
     . + g(n1,n2)*P(n3)*K(n5)*L(n4)
     . + g(n1,n2)*P(n4)*K(n3)*L(n5)
     . + g(n1,n2)*P(n4)*K(n5)*L(n3)
     . + g(n1,n2)*P(n5)*K(n3)*L(n4)
     . + g(n1,n2)*P(n5)*K(n4)*L(n3)
     . + g(n1,n3)*P(n2)*K(n4)*L(n5)
     . + g(n1,n3)*P(n2)*K(n5)*L(n4)
     . + g(n1,n3)*P(n4)*K(n2)*L(n5)
     . + g(n1,n3)*P(n4)*K(n5)*L(n2)
     . + g(n1,n3)*P(n5)*K(n2)*L(n4)
     . + g(n1,n3)*P(n5)*K(n4)*L(n2)
     . + g(n1,n4)*P(n2)*K(n3)*L(n5)
     . + g(n1,n4)*P(n2)*K(n5)*L(n3)
     . + g(n1,n4)*P(n3)*K(n2)*L(n5)
     . + g(n1,n4)*P(n3)*K(n5)*L(n2)
     . + g(n1,n4)*P(n5)*K(n2)*L(n3)
     . + g(n1,n4)*P(n5)*K(n3)*L(n2)
     . + g(n1,n5)*P(n2)*K(n3)*L(n4)
     . + g(n1,n5)*P(n2)*K(n4)*L(n3)
     . + g(n1,n5)*P(n3)*K(n2)*L(n4)
     . + g(n1,n5)*P(n3)*K(n4)*L(n2)
     . + g(n1,n5)*P(n4)*K(n2)*L(n3)
     . + g(n1,n5)*P(n4)*K(n3)*L(n2)
      pvSDDPKL=pvSDDPKL
     . + g(n2,n3)*P(n1)*K(n4)*L(n5)
     . + g(n2,n3)*P(n1)*K(n5)*L(n4)
     . + g(n2,n3)*P(n4)*K(n1)*L(n5)
     . + g(n2,n3)*P(n4)*K(n5)*L(n1)
     . + g(n2,n3)*P(n5)*K(n1)*L(n4)
     . + g(n2,n3)*P(n5)*K(n4)*L(n1)
     . + g(n2,n4)*P(n1)*K(n3)*L(n5)
     . + g(n2,n4)*P(n1)*K(n5)*L(n3)
     . + g(n2,n4)*P(n3)*K(n1)*L(n5)
     . + g(n2,n4)*P(n3)*K(n5)*L(n1)
     . + g(n2,n4)*P(n5)*K(n1)*L(n3)
     . + g(n2,n4)*P(n5)*K(n3)*L(n1)
     . + g(n2,n5)*P(n1)*K(n3)*L(n4)
     . + g(n2,n5)*P(n1)*K(n4)*L(n3)
     . + g(n2,n5)*P(n3)*K(n1)*L(n4)
     . + g(n2,n5)*P(n3)*K(n4)*L(n1)
     . + g(n2,n5)*P(n4)*K(n1)*L(n3)
     . + g(n2,n5)*P(n4)*K(n3)*L(n1)
     . + g(n3,n4)*P(n1)*K(n2)*L(n5)
     . + g(n3,n4)*P(n1)*K(n5)*L(n2)
     . + g(n3,n4)*P(n2)*K(n1)*L(n5)
     . + g(n3,n4)*P(n2)*K(n5)*L(n1)
     . + g(n3,n4)*P(n5)*K(n1)*L(n2)
     . + g(n3,n4)*P(n5)*K(n2)*L(n1)
     . + g(n3,n5)*P(n1)*K(n2)*L(n4)
     . + g(n3,n5)*P(n1)*K(n4)*L(n2)
     . + g(n3,n5)*P(n2)*K(n1)*L(n4)
     . + g(n3,n5)*P(n2)*K(n4)*L(n1)
     . + g(n3,n5)*P(n4)*K(n1)*L(n2)
     . + g(n3,n5)*P(n4)*K(n2)*L(n1)
     . + g(n4,n5)*P(n1)*K(n2)*L(n3)
     . + g(n4,n5)*P(n1)*K(n3)*L(n2)
     . + g(n4,n5)*P(n2)*K(n1)*L(n3)
     . + g(n4,n5)*P(n2)*K(n3)*L(n1)
     . + g(n4,n5)*P(n3)*K(n1)*L(n2)
     . + g(n4,n5)*P(n3)*K(n2)*L(n1)
      return
      end

      function pvSDDPPK(n1,n2,n3,n4,n5,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPPK
      integer n1,n2,n3,n4,n5
      real(dp):: P(4),K(4)
      include 'TRmetric.f'
      pvSDDPPK=
     . + g(n1,n2)*P(n3)*P(n4)*K(n5)
     . + g(n1,n2)*P(n3)*P(n5)*K(n4)
     . + g(n1,n2)*P(n4)*P(n5)*K(n3)
     . + g(n1,n3)*P(n2)*P(n4)*K(n5)
     . + g(n1,n3)*P(n2)*P(n5)*K(n4)
     . + g(n1,n3)*P(n4)*P(n5)*K(n2)
     . + g(n1,n4)*P(n2)*P(n3)*K(n5)
     . + g(n1,n4)*P(n2)*P(n5)*K(n3)
     . + g(n1,n4)*P(n3)*P(n5)*K(n2)
     . + g(n1,n5)*P(n2)*P(n3)*K(n4)
     . + g(n1,n5)*P(n2)*P(n4)*K(n3)
     . + g(n1,n5)*P(n3)*P(n4)*K(n2)
     . + g(n2,n3)*P(n1)*P(n4)*K(n5)
     . + g(n2,n3)*P(n1)*P(n5)*K(n4)
     . + g(n2,n3)*P(n4)*P(n5)*K(n1)
     . + g(n2,n4)*P(n1)*P(n3)*K(n5)
     . + g(n2,n4)*P(n1)*P(n5)*K(n3)
     . + g(n2,n4)*P(n3)*P(n5)*K(n1)
     . + g(n2,n5)*P(n1)*P(n3)*K(n4)
     . + g(n2,n5)*P(n1)*P(n4)*K(n3)
     . + g(n2,n5)*P(n3)*P(n4)*K(n1)
     . + g(n3,n4)*P(n1)*P(n2)*K(n5)
     . + g(n3,n4)*P(n1)*P(n5)*K(n2)
     . + g(n3,n4)*P(n2)*P(n5)*K(n1)
     . + g(n3,n5)*P(n1)*P(n2)*K(n4)
     . + g(n3,n5)*P(n1)*P(n4)*K(n2)
     . + g(n3,n5)*P(n2)*P(n4)*K(n1)
     . + g(n4,n5)*P(n1)*P(n2)*K(n3)
     . + g(n4,n5)*P(n1)*P(n3)*K(n2)
     . + g(n4,n5)*P(n2)*P(n3)*K(n1)
      return
      end

      function pvSDDDDP(n1,n2,n3,n4,n5,P)
      implicit none
      include 'types.f'
      real(dp):: pvSDDDDP
      integer n1,n2,n3,n4,n5
      real(dp):: P(4)
      include 'TRmetric.f'
      pvSDDDDP=
     . + g(n1,n2)*g(n3,n4)*P(n5)
     . + g(n1,n2)*g(n3,n5)*P(n4)
     . + g(n1,n2)*g(n4,n5)*P(n3)
     . + g(n1,n3)*g(n2,n4)*P(n5)
     . + g(n1,n3)*g(n2,n5)*P(n4)
     . + g(n1,n3)*g(n4,n5)*P(n2)
     . + g(n1,n4)*g(n2,n3)*P(n5)
     . + g(n1,n4)*g(n2,n5)*P(n3)
     . + g(n1,n4)*g(n3,n5)*P(n2)
     . + g(n1,n5)*g(n2,n3)*P(n4)
     . + g(n1,n5)*g(n2,n4)*P(n3)
     . + g(n1,n5)*g(n3,n4)*P(n2)
     . + g(n2,n3)*g(n4,n5)*P(n1)
     . + g(n2,n4)*g(n3,n5)*P(n1)
     . + g(n2,n5)*g(n3,n4)*P(n1)
      return
      end

      function pvSDDPPP(n1,n2,n3,n4,n5,P)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPPP
      integer n1,n2,n3,n4,n5
      real(dp):: P(4)
      include 'TRmetric.f'
      pvSDDPPP=
     . + g(n1,n2)*P(n3)*P(n4)*P(n5)
     . + g(n1,n3)*P(n2)*P(n4)*P(n5)
     . + g(n1,n4)*P(n2)*P(n3)*P(n5)
     . + g(n1,n5)*P(n2)*P(n3)*P(n4)
     . + g(n2,n3)*P(n1)*P(n4)*P(n5)
     . + g(n2,n4)*P(n1)*P(n3)*P(n5)
     . + g(n2,n5)*P(n1)*P(n3)*P(n4)
     . + g(n3,n4)*P(n1)*P(n2)*P(n5)
     . + g(n3,n5)*P(n1)*P(n2)*P(n4)
     . + g(n4,n5)*P(n1)*P(n2)*P(n3)
      return
      end

      function pvSPPPPK(n1,n2,n3,n4,n5,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPPK
      integer n1,n2,n3,n4,n5
      real(dp):: P(4),K(4)
      pvSPPPPK=
     . + P(n1)*P(n2)*P(n3)*P(n4)*K(n5)
     . + P(n1)*P(n2)*P(n3)*P(n5)*K(n4)
     . + P(n1)*P(n2)*P(n4)*P(n5)*K(n3)
     . + P(n1)*P(n3)*P(n4)*P(n5)*K(n2)
     . + P(n2)*P(n3)*P(n4)*P(n5)*K(n1)
      return
      end


      function pvSPK(n1,n2,q1,q2)
      implicit none
      include 'types.f'
      real(dp):: pvSPK
      integer n1,n2
      real(dp):: q1(4),q2(4)
      pvSPK=q1(n1)*q2(n2)+q2(n1)*q1(n2)
      return
      end

      function pvSPKL(n1,n2,n3,q1,q2,q3)
      implicit none
      include 'types.f'
      real(dp):: pvSPKL
      integer n1,n2,n3
      real(dp):: q1(4),q2(4),q3(4)
      pvSPKL=
     .   +q1(n1)*q2(n2)*q3(n3)+q2(n1)*q3(n2)*q1(n3)
     .   +q3(n1)*q1(n2)*q2(n3)+q3(n1)*q2(n2)*q1(n3)
     .   +q2(n1)*q1(n2)*q3(n3)+q1(n1)*q3(n2)*q2(n3)
      return
      end

      function pvSPKK(n1,n2,n3,q1,q2)
      implicit none
      include 'types.f'
      real(dp):: pvSPKK
      integer n1,n2,n3
      real(dp):: q1(4),q2(4)
      pvSPKK=
     .   +q1(n1)*q2(n2)*q2(n3)
     .   +q2(n1)*q2(n2)*q1(n3)
     .   +q2(n1)*q1(n2)*q2(n3)
      return
      end
      
      function pvSDDP(n1,n2,n3,q1)
      implicit none
      include 'types.f'
      real(dp):: pvSDDP
      integer n1,n2,n3
      real(dp):: q1(4)
      include 'TRmetric.f'
      pvSDDP=q1(n1)*g(n2,n3)+q1(n2)*g(n1,n3)+q1(n3)*g(n1,n2)
      return
      end

      function pvSPKKK(n1,n2,n3,n4,q1,q2)
      implicit none
      include 'types.f'
      real(dp):: pvSPKKK
      integer n1,n2,n3,n4
      real(dp):: q1(4),q2(4)
      pvSPKKK=
     . +q1(n1)*q2(n2)*q2(n3)*q2(n4)
     . +q1(n2)*q2(n3)*q2(n4)*q2(n1)
     . +q1(n3)*q2(n4)*q2(n1)*q2(n2)
     . +q1(n4)*q2(n1)*q2(n2)*q2(n3)
      return
      end


      function pvSPPKK(n1,n2,n3,n4,q1,q2)
      implicit none
      include 'types.f'
      real(dp):: pvSPPKK
      integer n1,n2,n3,n4
      real(dp):: q1(4),q2(4)
      pvSPPKK=
     . +q1(n1)*q1(n2)*q2(n3)*q2(n4)
     . +q1(n1)*q1(n3)*q2(n2)*q2(n4)
     . +q1(n1)*q1(n4)*q2(n2)*q2(n3)
     . +q1(n2)*q1(n3)*q2(n1)*q2(n4)
     . +q1(n2)*q1(n4)*q2(n1)*q2(n3)
     . +q1(n3)*q1(n4)*q2(n1)*q2(n2)
      return
      end

      function pvSPPKL(n1,n2,n3,n4,q1,q2,q3)
      implicit none
      include 'types.f'
      real(dp):: pvSPPKL
      integer n1,n2,n3,n4
      real(dp):: q1(4),q2(4),q3(4),pvSPK
      pvSPPKL=
     . +q1(n1)*q1(n2)*pvSPK(n3,n4,q2,q3)
     . +q1(n1)*q1(n3)*pvSPK(n2,n4,q2,q3)
     . +q1(n1)*q1(n4)*pvSPK(n2,n3,q2,q3)
     . +q1(n2)*q1(n3)*pvSPK(n1,n4,q2,q3)
     . +q1(n2)*q1(n4)*pvSPK(n1,n3,q2,q3)
     . +q1(n3)*q1(n4)*pvSPK(n1,n2,q2,q3)
      return
      end

      function pvSDDPP(n1,n2,n3,n4,q1)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPP
      integer n1,n2,n3,n4
      real(dp):: q1(4)
      include 'TRmetric.f'
      pvSDDPP=
     . +g(n1,n2)*q1(n3)*q1(n4)
     . +g(n1,n3)*q1(n2)*q1(n4)
     . +g(n1,n4)*q1(n2)*q1(n3)
     . +g(n2,n3)*q1(n1)*q1(n4)
     . +g(n2,n4)*q1(n1)*q1(n3)
     . +g(n3,n4)*q1(n1)*q1(n2)
      return
      end

      function pvSDDPK(n1,n2,n3,n4,q1,q2)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPK
      integer n1,n2,n3,n4
      real(dp):: q1(4),q2(4),pvSPK
      include 'TRmetric.f'
      pvSDDPK=
     . +g(n1,n2)*pvSPK(n3,n4,q1,q2)
     . +g(n1,n3)*pvSPK(n2,n4,q1,q2)
     . +g(n1,n4)*pvSPK(n2,n3,q1,q2)
     . +g(n2,n3)*pvSPK(n1,n4,q1,q2)
     . +g(n2,n4)*pvSPK(n1,n3,q1,q2)
     . +g(n3,n4)*pvSPK(n1,n2,q1,q2)
      return
      end

      function pvSDDDD(n1,n2,n3,n4)
      implicit none
      include 'types.f'
      real(dp):: pvSDDDD
      integer n1,n2,n3,n4
      include 'TRmetric.f'
      pvSDDDD=
     . g(n1,n2)*g(n3,n4)+g(n1,n3)*g(n2,n4)+g(n1,n4)*g(n2,n3)
      return
      end

      function pvSPKLM(n1,n2,n3,n4,q1,q2,q3,q4)
      implicit none
      include 'types.f'
      real(dp):: pvSPKLM
      integer n1,n2,n3,n4
      real(dp):: q1(4),q2(4),q3(4),q4(4)
      pvSPKLM=
     .   +q1(n1)*q2(n2)*q3(n3)*q4(n4)
     .   +q1(n1)*q2(n2)*q3(n4)*q4(n3)
     .   +q1(n1)*q2(n3)*q3(n2)*q4(n4)
     .   +q1(n1)*q2(n3)*q3(n4)*q4(n2)
     .   +q1(n1)*q2(n4)*q3(n2)*q4(n3)
     .   +q1(n1)*q2(n4)*q3(n3)*q4(n2)

     .   +q1(n2)*q2(n1)*q3(n3)*q4(n4)
     .   +q1(n2)*q2(n1)*q3(n4)*q4(n3)
     .   +q1(n2)*q2(n3)*q3(n1)*q4(n4)
     .   +q1(n2)*q2(n3)*q3(n4)*q4(n1)
     .   +q1(n2)*q2(n4)*q3(n1)*q4(n3)
     .   +q1(n2)*q2(n4)*q3(n3)*q4(n1)

     .   +q1(n3)*q2(n1)*q3(n2)*q4(n4)
     .   +q1(n3)*q2(n1)*q3(n4)*q4(n2)
     .   +q1(n3)*q2(n2)*q3(n1)*q4(n4)
     .   +q1(n3)*q2(n2)*q3(n4)*q4(n1)
     .   +q1(n3)*q2(n4)*q3(n1)*q4(n2)
     .   +q1(n3)*q2(n4)*q3(n2)*q4(n1)

     .   +q1(n4)*q2(n1)*q3(n2)*q4(n3)
     .   +q1(n4)*q2(n1)*q3(n3)*q4(n2)
     .   +q1(n4)*q2(n2)*q3(n1)*q4(n3)
     .   +q1(n4)*q2(n2)*q3(n3)*q4(n1)
     .   +q1(n4)*q2(n3)*q3(n1)*q4(n2)
     .   +q1(n4)*q2(n3)*q3(n2)*q4(n1)

      return
      end

      function pvSPPPPKK(n1,n2,n3,n4,n5,n6,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPPKK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      pvSPPPPKK=
     . + P(n1)*P(n2)*P(n3)*P(n4)*K(n5)*K(n6)
     . + P(n1)*P(n2)*P(n3)*P(n5)*K(n4)*K(n6)
     . + P(n1)*P(n2)*P(n3)*P(n6)*K(n4)*K(n5)
     . + P(n1)*P(n2)*P(n4)*P(n5)*K(n3)*K(n6)
     . + P(n1)*P(n2)*P(n4)*P(n6)*K(n3)*K(n5)
     . + P(n1)*P(n2)*P(n5)*P(n6)*K(n3)*K(n4)
     . + P(n1)*P(n3)*P(n4)*P(n5)*K(n2)*K(n6)
     . + P(n1)*P(n3)*P(n4)*P(n6)*K(n2)*K(n5)
     . + P(n1)*P(n3)*P(n5)*P(n6)*K(n2)*K(n4)
     . + P(n1)*P(n4)*P(n5)*P(n6)*K(n2)*K(n3)
     . + P(n2)*P(n3)*P(n4)*P(n5)*K(n1)*K(n6)
     . + P(n2)*P(n3)*P(n4)*P(n6)*K(n1)*K(n5)
     . + P(n2)*P(n3)*P(n5)*P(n6)*K(n1)*K(n4)
     . + P(n2)*P(n4)*P(n5)*P(n6)*K(n1)*K(n3)
     . + P(n3)*P(n4)*P(n5)*P(n6)*K(n1)*K(n2)
      end

      function pvSPPPKKL(n1,n2,n3,n4,n5,n6,P,K,L)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPKKL
      real(dp):: P(4),K(4),L(4)
      integer n1,n2,n3,n4,n5,n6
      pvSPPPKKL=
     . + P(n1)*P(n2)*P(n3)*K(n4)*K(n5)*L(n6)
     . + P(n1)*P(n2)*P(n3)*K(n4)*K(n6)*L(n5)
     . + P(n1)*P(n2)*P(n3)*K(n5)*K(n6)*L(n4)
     . + P(n1)*P(n2)*P(n4)*K(n3)*K(n5)*L(n6)
     . + P(n1)*P(n2)*P(n4)*K(n3)*K(n6)*L(n5)
     . + P(n1)*P(n2)*P(n4)*K(n5)*K(n6)*L(n3)
     . + P(n1)*P(n2)*P(n5)*K(n3)*K(n4)*L(n6)
     . + P(n1)*P(n2)*P(n5)*K(n3)*K(n6)*L(n4)
     . + P(n1)*P(n2)*P(n5)*K(n4)*K(n6)*L(n3)
     . + P(n1)*P(n2)*P(n6)*K(n3)*K(n4)*L(n5)
     . + P(n1)*P(n2)*P(n6)*K(n3)*K(n5)*L(n4)
     . + P(n1)*P(n2)*P(n6)*K(n4)*K(n5)*L(n3)
     . + P(n1)*P(n3)*P(n4)*K(n2)*K(n5)*L(n6)
     . + P(n1)*P(n3)*P(n4)*K(n2)*K(n6)*L(n5)
     . + P(n1)*P(n3)*P(n4)*K(n5)*K(n6)*L(n2)
     . + P(n1)*P(n3)*P(n5)*K(n2)*K(n4)*L(n6)
     . + P(n1)*P(n3)*P(n5)*K(n2)*K(n6)*L(n4)
     . + P(n1)*P(n3)*P(n5)*K(n4)*K(n6)*L(n2)
     . + P(n1)*P(n3)*P(n6)*K(n2)*K(n4)*L(n5)
     . + P(n1)*P(n3)*P(n6)*K(n2)*K(n5)*L(n4)
     . + P(n1)*P(n3)*P(n6)*K(n4)*K(n5)*L(n2)
     . + P(n1)*P(n4)*P(n5)*K(n2)*K(n3)*L(n6)
     . + P(n1)*P(n4)*P(n5)*K(n2)*K(n6)*L(n3)
     . + P(n1)*P(n4)*P(n5)*K(n3)*K(n6)*L(n2)
     . + P(n1)*P(n4)*P(n6)*K(n2)*K(n3)*L(n5)
     . + P(n1)*P(n4)*P(n6)*K(n2)*K(n5)*L(n3)
     . + P(n1)*P(n4)*P(n6)*K(n3)*K(n5)*L(n2)
     . + P(n1)*P(n5)*P(n6)*K(n2)*K(n3)*L(n4)
     . + P(n1)*P(n5)*P(n6)*K(n2)*K(n4)*L(n3)
     . + P(n1)*P(n5)*P(n6)*K(n3)*K(n4)*L(n2)
      pvSPPPKKL=pvSPPPKKL
     . + P(n2)*P(n3)*P(n4)*K(n1)*K(n5)*L(n6)
     . + P(n2)*P(n3)*P(n4)*K(n1)*K(n6)*L(n5)
     . + P(n2)*P(n3)*P(n4)*K(n5)*K(n6)*L(n1)
     . + P(n2)*P(n3)*P(n5)*K(n1)*K(n4)*L(n6)
     . + P(n2)*P(n3)*P(n5)*K(n1)*K(n6)*L(n4)
     . + P(n2)*P(n3)*P(n5)*K(n4)*K(n6)*L(n1)
     . + P(n2)*P(n3)*P(n6)*K(n1)*K(n4)*L(n5)
     . + P(n2)*P(n3)*P(n6)*K(n1)*K(n5)*L(n4)
     . + P(n2)*P(n3)*P(n6)*K(n4)*K(n5)*L(n1)
     . + P(n2)*P(n4)*P(n5)*K(n1)*K(n3)*L(n6)
     . + P(n2)*P(n4)*P(n5)*K(n1)*K(n6)*L(n3)
     . + P(n2)*P(n4)*P(n5)*K(n3)*K(n6)*L(n1)
     . + P(n2)*P(n4)*P(n6)*K(n1)*K(n3)*L(n5)
     . + P(n2)*P(n4)*P(n6)*K(n1)*K(n5)*L(n3)
     . + P(n2)*P(n4)*P(n6)*K(n3)*K(n5)*L(n1)
     . + P(n2)*P(n5)*P(n6)*K(n1)*K(n3)*L(n4)
     . + P(n2)*P(n5)*P(n6)*K(n1)*K(n4)*L(n3)
     . + P(n2)*P(n5)*P(n6)*K(n3)*K(n4)*L(n1)
     . + P(n3)*P(n4)*P(n5)*K(n1)*K(n2)*L(n6)
     . + P(n3)*P(n4)*P(n5)*K(n1)*K(n6)*L(n2)
     . + P(n3)*P(n4)*P(n5)*K(n2)*K(n6)*L(n1)
     . + P(n3)*P(n4)*P(n6)*K(n1)*K(n2)*L(n5)
     . + P(n3)*P(n4)*P(n6)*K(n1)*K(n5)*L(n2)
     . + P(n3)*P(n4)*P(n6)*K(n2)*K(n5)*L(n1)
     . + P(n3)*P(n5)*P(n6)*K(n1)*K(n2)*L(n4)
     . + P(n3)*P(n5)*P(n6)*K(n1)*K(n4)*L(n2)
     . + P(n3)*P(n5)*P(n6)*K(n2)*K(n4)*L(n1)
     . + P(n4)*P(n5)*P(n6)*K(n1)*K(n2)*L(n3)
     . + P(n4)*P(n5)*P(n6)*K(n1)*K(n3)*L(n2)
     . + P(n4)*P(n5)*P(n6)*K(n2)*K(n3)*L(n1)
      end

      function pvSPPKKLL(n1,n2,n3,n4,n5,n6,P,K,L)
      implicit none
      include 'types.f'
      real(dp):: pvSPPKKLL
      real(dp):: P(4),K(4),L(4)
      integer n1,n2,n3,n4,n5,n6
      pvSPPKKLL=
     . + P(n1)*P(n2)*K(n3)*K(n4)*L(n5)*L(n6)
     . + P(n1)*P(n2)*K(n3)*K(n5)*L(n4)*L(n6)
     . + P(n1)*P(n2)*K(n3)*K(n6)*L(n4)*L(n5)
     . + P(n1)*P(n2)*K(n4)*K(n5)*L(n3)*L(n6)
     . + P(n1)*P(n2)*K(n4)*K(n6)*L(n3)*L(n5)
     . + P(n1)*P(n2)*K(n5)*K(n6)*L(n3)*L(n4)
     . + P(n1)*P(n3)*K(n2)*K(n4)*L(n5)*L(n6)
     . + P(n1)*P(n3)*K(n2)*K(n5)*L(n4)*L(n6)
     . + P(n1)*P(n3)*K(n2)*K(n6)*L(n4)*L(n5)
     . + P(n1)*P(n3)*K(n4)*K(n5)*L(n2)*L(n6)
     . + P(n1)*P(n3)*K(n4)*K(n6)*L(n2)*L(n5)
     . + P(n1)*P(n3)*K(n5)*K(n6)*L(n2)*L(n4)
     . + P(n1)*P(n4)*K(n2)*K(n3)*L(n5)*L(n6)
     . + P(n1)*P(n4)*K(n2)*K(n5)*L(n3)*L(n6)
     . + P(n1)*P(n4)*K(n2)*K(n6)*L(n3)*L(n5)
     . + P(n1)*P(n4)*K(n3)*K(n5)*L(n2)*L(n6)
     . + P(n1)*P(n4)*K(n3)*K(n6)*L(n2)*L(n5)
     . + P(n1)*P(n4)*K(n5)*K(n6)*L(n2)*L(n3)
     . + P(n1)*P(n5)*K(n2)*K(n3)*L(n4)*L(n6)
     . + P(n1)*P(n5)*K(n2)*K(n4)*L(n3)*L(n6)
     . + P(n1)*P(n5)*K(n2)*K(n6)*L(n3)*L(n4)
     . + P(n1)*P(n5)*K(n3)*K(n4)*L(n2)*L(n6)
     . + P(n1)*P(n5)*K(n3)*K(n6)*L(n2)*L(n4)
     . + P(n1)*P(n5)*K(n4)*K(n6)*L(n2)*L(n3)
     . + P(n1)*P(n6)*K(n2)*K(n3)*L(n4)*L(n5)
     . + P(n1)*P(n6)*K(n2)*K(n4)*L(n3)*L(n5)
     . + P(n1)*P(n6)*K(n2)*K(n5)*L(n3)*L(n4)
     . + P(n1)*P(n6)*K(n3)*K(n4)*L(n2)*L(n5)
     . + P(n1)*P(n6)*K(n3)*K(n5)*L(n2)*L(n4)
     . + P(n1)*P(n6)*K(n4)*K(n5)*L(n2)*L(n3)
      pvSPPKKLL=pvSPPKKLL
     . + P(n2)*P(n3)*K(n1)*K(n4)*L(n5)*L(n6)
     . + P(n2)*P(n3)*K(n1)*K(n5)*L(n4)*L(n6)
     . + P(n2)*P(n3)*K(n1)*K(n6)*L(n4)*L(n5)
     . + P(n2)*P(n3)*K(n4)*K(n5)*L(n1)*L(n6)
     . + P(n2)*P(n3)*K(n4)*K(n6)*L(n1)*L(n5)
     . + P(n2)*P(n3)*K(n5)*K(n6)*L(n1)*L(n4)
     . + P(n2)*P(n4)*K(n1)*K(n3)*L(n5)*L(n6)
     . + P(n2)*P(n4)*K(n1)*K(n5)*L(n3)*L(n6)
     . + P(n2)*P(n4)*K(n1)*K(n6)*L(n3)*L(n5)
     . + P(n2)*P(n4)*K(n3)*K(n5)*L(n1)*L(n6)
     . + P(n2)*P(n4)*K(n3)*K(n6)*L(n1)*L(n5)
     . + P(n2)*P(n4)*K(n5)*K(n6)*L(n1)*L(n3)
     . + P(n2)*P(n5)*K(n1)*K(n3)*L(n4)*L(n6)
     . + P(n2)*P(n5)*K(n1)*K(n4)*L(n3)*L(n6)
     . + P(n2)*P(n5)*K(n1)*K(n6)*L(n3)*L(n4)
     . + P(n2)*P(n5)*K(n3)*K(n4)*L(n1)*L(n6)
     . + P(n2)*P(n5)*K(n3)*K(n6)*L(n1)*L(n4)
     . + P(n2)*P(n5)*K(n4)*K(n6)*L(n1)*L(n3)
     . + P(n2)*P(n6)*K(n1)*K(n3)*L(n4)*L(n5)
     . + P(n2)*P(n6)*K(n1)*K(n4)*L(n3)*L(n5)
     . + P(n2)*P(n6)*K(n1)*K(n5)*L(n3)*L(n4)
     . + P(n2)*P(n6)*K(n3)*K(n4)*L(n1)*L(n5)
     . + P(n2)*P(n6)*K(n3)*K(n5)*L(n1)*L(n4)
     . + P(n2)*P(n6)*K(n4)*K(n5)*L(n1)*L(n3)
      pvSPPKKLL=pvSPPKKLL
     . + P(n3)*P(n4)*K(n1)*K(n2)*L(n5)*L(n6)
     . + P(n3)*P(n4)*K(n1)*K(n5)*L(n2)*L(n6)
     . + P(n3)*P(n4)*K(n1)*K(n6)*L(n2)*L(n5)
     . + P(n3)*P(n4)*K(n2)*K(n5)*L(n1)*L(n6)
     . + P(n3)*P(n4)*K(n2)*K(n6)*L(n1)*L(n5)
     . + P(n3)*P(n4)*K(n5)*K(n6)*L(n1)*L(n2)
     . + P(n3)*P(n5)*K(n1)*K(n2)*L(n4)*L(n6)
     . + P(n3)*P(n5)*K(n1)*K(n4)*L(n2)*L(n6)
     . + P(n3)*P(n5)*K(n1)*K(n6)*L(n2)*L(n4)
     . + P(n3)*P(n5)*K(n2)*K(n4)*L(n1)*L(n6)
     . + P(n3)*P(n5)*K(n2)*K(n6)*L(n1)*L(n4)
     . + P(n3)*P(n5)*K(n4)*K(n6)*L(n1)*L(n2)
     . + P(n3)*P(n6)*K(n1)*K(n2)*L(n4)*L(n5)
     . + P(n3)*P(n6)*K(n1)*K(n4)*L(n2)*L(n5)
     . + P(n3)*P(n6)*K(n1)*K(n5)*L(n2)*L(n4)
     . + P(n3)*P(n6)*K(n2)*K(n4)*L(n1)*L(n5)
     . + P(n3)*P(n6)*K(n2)*K(n5)*L(n1)*L(n4)
     . + P(n3)*P(n6)*K(n4)*K(n5)*L(n1)*L(n2)
      pvSPPKKLL=pvSPPKKLL
     . + P(n4)*P(n5)*K(n1)*K(n2)*L(n3)*L(n6)
     . + P(n4)*P(n5)*K(n1)*K(n3)*L(n2)*L(n6)
     . + P(n4)*P(n5)*K(n1)*K(n6)*L(n2)*L(n3)
     . + P(n4)*P(n5)*K(n2)*K(n3)*L(n1)*L(n6)
     . + P(n4)*P(n5)*K(n2)*K(n6)*L(n1)*L(n3)
     . + P(n4)*P(n5)*K(n3)*K(n6)*L(n1)*L(n2)
     . + P(n4)*P(n6)*K(n1)*K(n2)*L(n3)*L(n5)
     . + P(n4)*P(n6)*K(n1)*K(n3)*L(n2)*L(n5)
     . + P(n4)*P(n6)*K(n1)*K(n5)*L(n2)*L(n3)
     . + P(n4)*P(n6)*K(n2)*K(n3)*L(n1)*L(n5)
     . + P(n4)*P(n6)*K(n2)*K(n5)*L(n1)*L(n3)
     . + P(n4)*P(n6)*K(n3)*K(n5)*L(n1)*L(n2)
     . + P(n5)*P(n6)*K(n1)*K(n2)*L(n3)*L(n4)
     . + P(n5)*P(n6)*K(n1)*K(n3)*L(n2)*L(n4)
     . + P(n5)*P(n6)*K(n1)*K(n4)*L(n2)*L(n3)
     . + P(n5)*P(n6)*K(n2)*K(n3)*L(n1)*L(n4)
     . + P(n5)*P(n6)*K(n2)*K(n4)*L(n1)*L(n3)
     . + P(n5)*P(n6)*K(n3)*K(n4)*L(n1)*L(n2)
      end

      function pvSDDPPKK(n1,n2,n3,n4,n5,n6,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPPKK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDPPKK=
     . + g(n1,n2)*P(n3)*P(n4)*K(n5)*K(n6)
     . + g(n1,n2)*P(n3)*P(n5)*K(n4)*K(n6)
     . + g(n1,n2)*P(n3)*P(n6)*K(n4)*K(n5)
     . + g(n1,n2)*P(n4)*P(n5)*K(n3)*K(n6)
     . + g(n1,n2)*P(n4)*P(n6)*K(n3)*K(n5)
     . + g(n1,n2)*P(n5)*P(n6)*K(n3)*K(n4)
     . + g(n1,n3)*P(n2)*P(n4)*K(n5)*K(n6)
     . + g(n1,n3)*P(n2)*P(n5)*K(n4)*K(n6)
     . + g(n1,n3)*P(n2)*P(n6)*K(n4)*K(n5)
     . + g(n1,n3)*P(n4)*P(n5)*K(n2)*K(n6)
     . + g(n1,n3)*P(n4)*P(n6)*K(n2)*K(n5)
     . + g(n1,n3)*P(n5)*P(n6)*K(n2)*K(n4)
     . + g(n1,n4)*P(n2)*P(n3)*K(n5)*K(n6)
     . + g(n1,n4)*P(n2)*P(n5)*K(n3)*K(n6)
     . + g(n1,n4)*P(n2)*P(n6)*K(n3)*K(n5)
     . + g(n1,n4)*P(n3)*P(n5)*K(n2)*K(n6)
     . + g(n1,n4)*P(n3)*P(n6)*K(n2)*K(n5)
     . + g(n1,n4)*P(n5)*P(n6)*K(n2)*K(n3)
     . + g(n1,n5)*P(n2)*P(n3)*K(n4)*K(n6)
     . + g(n1,n5)*P(n2)*P(n4)*K(n3)*K(n6)
     . + g(n1,n5)*P(n2)*P(n6)*K(n3)*K(n4)
     . + g(n1,n5)*P(n3)*P(n4)*K(n2)*K(n6)
     . + g(n1,n5)*P(n3)*P(n6)*K(n2)*K(n4)
     . + g(n1,n5)*P(n4)*P(n6)*K(n2)*K(n3)
     . + g(n1,n6)*P(n2)*P(n3)*K(n4)*K(n5)
     . + g(n1,n6)*P(n2)*P(n4)*K(n3)*K(n5)
     . + g(n1,n6)*P(n2)*P(n5)*K(n3)*K(n4)
     . + g(n1,n6)*P(n3)*P(n4)*K(n2)*K(n5)
     . + g(n1,n6)*P(n3)*P(n5)*K(n2)*K(n4)
     . + g(n1,n6)*P(n4)*P(n5)*K(n2)*K(n3)
      pvSDDPPKK=pvSDDPPKK
     . + g(n2,n3)*P(n1)*P(n4)*K(n5)*K(n6)
     . + g(n2,n3)*P(n1)*P(n5)*K(n4)*K(n6)
     . + g(n2,n3)*P(n1)*P(n6)*K(n4)*K(n5)
     . + g(n2,n3)*P(n4)*P(n5)*K(n1)*K(n6)
     . + g(n2,n3)*P(n4)*P(n6)*K(n1)*K(n5)
     . + g(n2,n3)*P(n5)*P(n6)*K(n1)*K(n4)
     . + g(n2,n4)*P(n1)*P(n3)*K(n5)*K(n6)
     . + g(n2,n4)*P(n1)*P(n5)*K(n3)*K(n6)
     . + g(n2,n4)*P(n1)*P(n6)*K(n3)*K(n5)
     . + g(n2,n4)*P(n3)*P(n5)*K(n1)*K(n6)
     . + g(n2,n4)*P(n3)*P(n6)*K(n1)*K(n5)
     . + g(n2,n4)*P(n5)*P(n6)*K(n1)*K(n3)
     . + g(n2,n5)*P(n1)*P(n3)*K(n4)*K(n6)
     . + g(n2,n5)*P(n1)*P(n4)*K(n3)*K(n6)
     . + g(n2,n5)*P(n1)*P(n6)*K(n3)*K(n4)
     . + g(n2,n5)*P(n3)*P(n4)*K(n1)*K(n6)
     . + g(n2,n5)*P(n3)*P(n6)*K(n1)*K(n4)
     . + g(n2,n5)*P(n4)*P(n6)*K(n1)*K(n3)
     . + g(n2,n6)*P(n1)*P(n3)*K(n4)*K(n5)
     . + g(n2,n6)*P(n1)*P(n4)*K(n3)*K(n5)
     . + g(n2,n6)*P(n1)*P(n5)*K(n3)*K(n4)
     . + g(n2,n6)*P(n3)*P(n4)*K(n1)*K(n5)
     . + g(n2,n6)*P(n3)*P(n5)*K(n1)*K(n4)
     . + g(n2,n6)*P(n4)*P(n5)*K(n1)*K(n3)
     . + g(n3,n4)*P(n1)*P(n2)*K(n5)*K(n6)
     . + g(n3,n4)*P(n1)*P(n5)*K(n2)*K(n6)
     . + g(n3,n4)*P(n1)*P(n6)*K(n2)*K(n5)
     . + g(n3,n4)*P(n2)*P(n5)*K(n1)*K(n6)
     . + g(n3,n4)*P(n2)*P(n6)*K(n1)*K(n5)
     . + g(n3,n4)*P(n5)*P(n6)*K(n1)*K(n2)
     . + g(n3,n5)*P(n1)*P(n2)*K(n4)*K(n6)
     . + g(n3,n5)*P(n1)*P(n4)*K(n2)*K(n6)
     . + g(n3,n5)*P(n1)*P(n6)*K(n2)*K(n4)
     . + g(n3,n5)*P(n2)*P(n4)*K(n1)*K(n6)
     . + g(n3,n5)*P(n2)*P(n6)*K(n1)*K(n4)
     . + g(n3,n5)*P(n4)*P(n6)*K(n1)*K(n2)
     . + g(n3,n6)*P(n1)*P(n2)*K(n4)*K(n5)
     . + g(n3,n6)*P(n1)*P(n4)*K(n2)*K(n5)
     . + g(n3,n6)*P(n1)*P(n5)*K(n2)*K(n4)
     . + g(n3,n6)*P(n2)*P(n4)*K(n1)*K(n5)
     . + g(n3,n6)*P(n2)*P(n5)*K(n1)*K(n4)
     . + g(n3,n6)*P(n4)*P(n5)*K(n1)*K(n2)
      pvSDDPPKK=pvSDDPPKK
     . + g(n4,n5)*P(n1)*P(n2)*K(n3)*K(n6)
     . + g(n4,n5)*P(n1)*P(n3)*K(n2)*K(n6)
     . + g(n4,n5)*P(n1)*P(n6)*K(n2)*K(n3)
     . + g(n4,n5)*P(n2)*P(n3)*K(n1)*K(n6)
     . + g(n4,n5)*P(n2)*P(n6)*K(n1)*K(n3)
     . + g(n4,n5)*P(n3)*P(n6)*K(n1)*K(n2)
     . + g(n4,n6)*P(n1)*P(n2)*K(n3)*K(n5)
     . + g(n4,n6)*P(n1)*P(n3)*K(n2)*K(n5)
     . + g(n4,n6)*P(n1)*P(n5)*K(n2)*K(n3)
     . + g(n4,n6)*P(n2)*P(n3)*K(n1)*K(n5)
     . + g(n4,n6)*P(n2)*P(n5)*K(n1)*K(n3)
     . + g(n4,n6)*P(n3)*P(n5)*K(n1)*K(n2)
     . + g(n5,n6)*P(n1)*P(n2)*K(n3)*K(n4)
     . + g(n5,n6)*P(n1)*P(n3)*K(n2)*K(n4)
     . + g(n5,n6)*P(n1)*P(n4)*K(n2)*K(n3)
     . + g(n5,n6)*P(n2)*P(n3)*K(n1)*K(n4)
     . + g(n5,n6)*P(n2)*P(n4)*K(n1)*K(n3)
     . + g(n5,n6)*P(n3)*P(n4)*K(n1)*K(n2)
      end

      function pvSDDPPKL(n1,n2,n3,n4,n5,n6,P,K,L)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPPKL
      real(dp):: P(4),K(4),L(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDPPKL=
     . + g(n1,n2)*P(n3)*P(n4)*K(n5)*L(n6)
     . + g(n1,n2)*P(n3)*P(n4)*K(n6)*L(n5)
     . + g(n1,n2)*P(n3)*P(n5)*K(n4)*L(n6)
     . + g(n1,n2)*P(n3)*P(n5)*K(n6)*L(n4)
     . + g(n1,n2)*P(n3)*P(n6)*K(n4)*L(n5)
     . + g(n1,n2)*P(n3)*P(n6)*K(n5)*L(n4)
     . + g(n1,n2)*P(n4)*P(n5)*K(n3)*L(n6)
     . + g(n1,n2)*P(n4)*P(n5)*K(n6)*L(n3)
     . + g(n1,n2)*P(n4)*P(n6)*K(n3)*L(n5)
     . + g(n1,n2)*P(n4)*P(n6)*K(n5)*L(n3)
     . + g(n1,n2)*P(n5)*P(n6)*K(n3)*L(n4)
     . + g(n1,n2)*P(n5)*P(n6)*K(n4)*L(n3)
     . + g(n1,n3)*P(n2)*P(n4)*K(n5)*L(n6)
     . + g(n1,n3)*P(n2)*P(n4)*K(n6)*L(n5)
     . + g(n1,n3)*P(n2)*P(n5)*K(n4)*L(n6)
     . + g(n1,n3)*P(n2)*P(n5)*K(n6)*L(n4)
     . + g(n1,n3)*P(n2)*P(n6)*K(n4)*L(n5)
     . + g(n1,n3)*P(n2)*P(n6)*K(n5)*L(n4)
     . + g(n1,n3)*P(n4)*P(n5)*K(n2)*L(n6)
     . + g(n1,n3)*P(n4)*P(n5)*K(n6)*L(n2)
     . + g(n1,n3)*P(n4)*P(n6)*K(n2)*L(n5)
     . + g(n1,n3)*P(n4)*P(n6)*K(n5)*L(n2)
     . + g(n1,n3)*P(n5)*P(n6)*K(n2)*L(n4)
     . + g(n1,n3)*P(n5)*P(n6)*K(n4)*L(n2)
     . + g(n1,n4)*P(n2)*P(n3)*K(n5)*L(n6)
     . + g(n1,n4)*P(n2)*P(n3)*K(n6)*L(n5)
     . + g(n1,n4)*P(n2)*P(n5)*K(n3)*L(n6)
     . + g(n1,n4)*P(n2)*P(n5)*K(n6)*L(n3)
     . + g(n1,n4)*P(n2)*P(n6)*K(n3)*L(n5)
     . + g(n1,n4)*P(n2)*P(n6)*K(n5)*L(n3)
     . + g(n1,n4)*P(n3)*P(n5)*K(n2)*L(n6)
     . + g(n1,n4)*P(n3)*P(n5)*K(n6)*L(n2)
     . + g(n1,n4)*P(n3)*P(n6)*K(n2)*L(n5)
     . + g(n1,n4)*P(n3)*P(n6)*K(n5)*L(n2)
     . + g(n1,n4)*P(n5)*P(n6)*K(n2)*L(n3)
     . + g(n1,n4)*P(n5)*P(n6)*K(n3)*L(n2)
      pvSDDPPKL=pvSDDPPKL
     . + g(n1,n5)*P(n2)*P(n3)*K(n4)*L(n6)
     . + g(n1,n5)*P(n2)*P(n3)*K(n6)*L(n4)
     . + g(n1,n5)*P(n2)*P(n4)*K(n3)*L(n6)
     . + g(n1,n5)*P(n2)*P(n4)*K(n6)*L(n3)
     . + g(n1,n5)*P(n2)*P(n6)*K(n3)*L(n4)
     . + g(n1,n5)*P(n2)*P(n6)*K(n4)*L(n3)
     . + g(n1,n5)*P(n3)*P(n4)*K(n2)*L(n6)
     . + g(n1,n5)*P(n3)*P(n4)*K(n6)*L(n2)
     . + g(n1,n5)*P(n3)*P(n6)*K(n2)*L(n4)
     . + g(n1,n5)*P(n3)*P(n6)*K(n4)*L(n2)
     . + g(n1,n5)*P(n4)*P(n6)*K(n2)*L(n3)
     . + g(n1,n5)*P(n4)*P(n6)*K(n3)*L(n2)
     . + g(n1,n6)*P(n2)*P(n3)*K(n4)*L(n5)
     . + g(n1,n6)*P(n2)*P(n3)*K(n5)*L(n4)
     . + g(n1,n6)*P(n2)*P(n4)*K(n3)*L(n5)
     . + g(n1,n6)*P(n2)*P(n4)*K(n5)*L(n3)
     . + g(n1,n6)*P(n2)*P(n5)*K(n3)*L(n4)
     . + g(n1,n6)*P(n2)*P(n5)*K(n4)*L(n3)
     . + g(n1,n6)*P(n3)*P(n4)*K(n2)*L(n5)
     . + g(n1,n6)*P(n3)*P(n4)*K(n5)*L(n2)
     . + g(n1,n6)*P(n3)*P(n5)*K(n2)*L(n4)
     . + g(n1,n6)*P(n3)*P(n5)*K(n4)*L(n2)
     . + g(n1,n6)*P(n4)*P(n5)*K(n2)*L(n3)
     . + g(n1,n6)*P(n4)*P(n5)*K(n3)*L(n2)
     . + g(n2,n3)*P(n1)*P(n4)*K(n5)*L(n6)
     . + g(n2,n3)*P(n1)*P(n4)*K(n6)*L(n5)
     . + g(n2,n3)*P(n1)*P(n5)*K(n4)*L(n6)
     . + g(n2,n3)*P(n1)*P(n5)*K(n6)*L(n4)
     . + g(n2,n3)*P(n1)*P(n6)*K(n4)*L(n5)
     . + g(n2,n3)*P(n1)*P(n6)*K(n5)*L(n4)
     . + g(n2,n3)*P(n4)*P(n5)*K(n1)*L(n6)
     . + g(n2,n3)*P(n4)*P(n5)*K(n6)*L(n1)
     . + g(n2,n3)*P(n4)*P(n6)*K(n1)*L(n5)
     . + g(n2,n3)*P(n4)*P(n6)*K(n5)*L(n1)
     . + g(n2,n3)*P(n5)*P(n6)*K(n1)*L(n4)
     . + g(n2,n3)*P(n5)*P(n6)*K(n4)*L(n1)
      pvSDDPPKL=pvSDDPPKL
     . + g(n2,n4)*P(n1)*P(n3)*K(n5)*L(n6)
     . + g(n2,n4)*P(n1)*P(n3)*K(n6)*L(n5)
     . + g(n2,n4)*P(n1)*P(n5)*K(n3)*L(n6)
     . + g(n2,n4)*P(n1)*P(n5)*K(n6)*L(n3)
     . + g(n2,n4)*P(n1)*P(n6)*K(n3)*L(n5)
     . + g(n2,n4)*P(n1)*P(n6)*K(n5)*L(n3)
     . + g(n2,n4)*P(n3)*P(n5)*K(n1)*L(n6)
     . + g(n2,n4)*P(n3)*P(n5)*K(n6)*L(n1)
     . + g(n2,n4)*P(n3)*P(n6)*K(n1)*L(n5)
     . + g(n2,n4)*P(n3)*P(n6)*K(n5)*L(n1)
     . + g(n2,n4)*P(n5)*P(n6)*K(n1)*L(n3)
     . + g(n2,n4)*P(n5)*P(n6)*K(n3)*L(n1)
     . + g(n2,n5)*P(n1)*P(n3)*K(n4)*L(n6)
     . + g(n2,n5)*P(n1)*P(n3)*K(n6)*L(n4)
     . + g(n2,n5)*P(n1)*P(n4)*K(n3)*L(n6)
     . + g(n2,n5)*P(n1)*P(n4)*K(n6)*L(n3)
     . + g(n2,n5)*P(n1)*P(n6)*K(n3)*L(n4)
     . + g(n2,n5)*P(n1)*P(n6)*K(n4)*L(n3)
     . + g(n2,n5)*P(n3)*P(n4)*K(n1)*L(n6)
     . + g(n2,n5)*P(n3)*P(n4)*K(n6)*L(n1)
     . + g(n2,n5)*P(n3)*P(n6)*K(n1)*L(n4)
     . + g(n2,n5)*P(n3)*P(n6)*K(n4)*L(n1)
     . + g(n2,n5)*P(n4)*P(n6)*K(n1)*L(n3)
     . + g(n2,n5)*P(n4)*P(n6)*K(n3)*L(n1)
     . + g(n2,n6)*P(n1)*P(n3)*K(n4)*L(n5)
     . + g(n2,n6)*P(n1)*P(n3)*K(n5)*L(n4)
     . + g(n2,n6)*P(n1)*P(n4)*K(n3)*L(n5)
     . + g(n2,n6)*P(n1)*P(n4)*K(n5)*L(n3)
     . + g(n2,n6)*P(n1)*P(n5)*K(n3)*L(n4)
     . + g(n2,n6)*P(n1)*P(n5)*K(n4)*L(n3)
     . + g(n2,n6)*P(n3)*P(n4)*K(n1)*L(n5)
     . + g(n2,n6)*P(n3)*P(n4)*K(n5)*L(n1)
     . + g(n2,n6)*P(n3)*P(n5)*K(n1)*L(n4)
     . + g(n2,n6)*P(n3)*P(n5)*K(n4)*L(n1)
     . + g(n2,n6)*P(n4)*P(n5)*K(n1)*L(n3)
     . + g(n2,n6)*P(n4)*P(n5)*K(n3)*L(n1)
      pvSDDPPKL=pvSDDPPKL
     . + g(n3,n4)*P(n1)*P(n2)*K(n5)*L(n6)
     . + g(n3,n4)*P(n1)*P(n2)*K(n6)*L(n5)
     . + g(n3,n4)*P(n1)*P(n5)*K(n2)*L(n6)
     . + g(n3,n4)*P(n1)*P(n5)*K(n6)*L(n2)
     . + g(n3,n4)*P(n1)*P(n6)*K(n2)*L(n5)
     . + g(n3,n4)*P(n1)*P(n6)*K(n5)*L(n2)
     . + g(n3,n4)*P(n2)*P(n5)*K(n1)*L(n6)
     . + g(n3,n4)*P(n2)*P(n5)*K(n6)*L(n1)
     . + g(n3,n4)*P(n2)*P(n6)*K(n1)*L(n5)
     . + g(n3,n4)*P(n2)*P(n6)*K(n5)*L(n1)
     . + g(n3,n4)*P(n5)*P(n6)*K(n1)*L(n2)
     . + g(n3,n4)*P(n5)*P(n6)*K(n2)*L(n1)
     . + g(n3,n5)*P(n1)*P(n2)*K(n4)*L(n6)
     . + g(n3,n5)*P(n1)*P(n2)*K(n6)*L(n4)
     . + g(n3,n5)*P(n1)*P(n4)*K(n2)*L(n6)
     . + g(n3,n5)*P(n1)*P(n4)*K(n6)*L(n2)
     . + g(n3,n5)*P(n1)*P(n6)*K(n2)*L(n4)
     . + g(n3,n5)*P(n1)*P(n6)*K(n4)*L(n2)
     . + g(n3,n5)*P(n2)*P(n4)*K(n1)*L(n6)
     . + g(n3,n5)*P(n2)*P(n4)*K(n6)*L(n1)
     . + g(n3,n5)*P(n2)*P(n6)*K(n1)*L(n4)
     . + g(n3,n5)*P(n2)*P(n6)*K(n4)*L(n1)
     . + g(n3,n5)*P(n4)*P(n6)*K(n1)*L(n2)
     . + g(n3,n5)*P(n4)*P(n6)*K(n2)*L(n1)
     . + g(n3,n6)*P(n1)*P(n2)*K(n4)*L(n5)
     . + g(n3,n6)*P(n1)*P(n2)*K(n5)*L(n4)
     . + g(n3,n6)*P(n1)*P(n4)*K(n2)*L(n5)
     . + g(n3,n6)*P(n1)*P(n4)*K(n5)*L(n2)
     . + g(n3,n6)*P(n1)*P(n5)*K(n2)*L(n4)
     . + g(n3,n6)*P(n1)*P(n5)*K(n4)*L(n2)
     . + g(n3,n6)*P(n2)*P(n4)*K(n1)*L(n5)
     . + g(n3,n6)*P(n2)*P(n4)*K(n5)*L(n1)
     . + g(n3,n6)*P(n2)*P(n5)*K(n1)*L(n4)
     . + g(n3,n6)*P(n2)*P(n5)*K(n4)*L(n1)
     . + g(n3,n6)*P(n4)*P(n5)*K(n1)*L(n2)
     . + g(n3,n6)*P(n4)*P(n5)*K(n2)*L(n1)
      pvSDDPPKL=pvSDDPPKL
     . + g(n4,n5)*P(n1)*P(n2)*K(n3)*L(n6)
     . + g(n4,n5)*P(n1)*P(n2)*K(n6)*L(n3)
     . + g(n4,n5)*P(n1)*P(n3)*K(n2)*L(n6)
     . + g(n4,n5)*P(n1)*P(n3)*K(n6)*L(n2)
     . + g(n4,n5)*P(n1)*P(n6)*K(n2)*L(n3)
     . + g(n4,n5)*P(n1)*P(n6)*K(n3)*L(n2)
     . + g(n4,n5)*P(n2)*P(n3)*K(n1)*L(n6)
     . + g(n4,n5)*P(n2)*P(n3)*K(n6)*L(n1)
     . + g(n4,n5)*P(n2)*P(n6)*K(n1)*L(n3)
     . + g(n4,n5)*P(n2)*P(n6)*K(n3)*L(n1)
     . + g(n4,n5)*P(n3)*P(n6)*K(n1)*L(n2)
     . + g(n4,n5)*P(n3)*P(n6)*K(n2)*L(n1)
     . + g(n4,n6)*P(n1)*P(n2)*K(n3)*L(n5)
     . + g(n4,n6)*P(n1)*P(n2)*K(n5)*L(n3)
     . + g(n4,n6)*P(n1)*P(n3)*K(n2)*L(n5)
     . + g(n4,n6)*P(n1)*P(n3)*K(n5)*L(n2)
     . + g(n4,n6)*P(n1)*P(n5)*K(n2)*L(n3)
     . + g(n4,n6)*P(n1)*P(n5)*K(n3)*L(n2)
     . + g(n4,n6)*P(n2)*P(n3)*K(n1)*L(n5)
     . + g(n4,n6)*P(n2)*P(n3)*K(n5)*L(n1)
     . + g(n4,n6)*P(n2)*P(n5)*K(n1)*L(n3)
     . + g(n4,n6)*P(n2)*P(n5)*K(n3)*L(n1)
     . + g(n4,n6)*P(n3)*P(n5)*K(n1)*L(n2)
     . + g(n4,n6)*P(n3)*P(n5)*K(n2)*L(n1)
     . + g(n5,n6)*P(n1)*P(n2)*K(n3)*L(n4)
     . + g(n5,n6)*P(n1)*P(n2)*K(n4)*L(n3)
     . + g(n5,n6)*P(n1)*P(n3)*K(n2)*L(n4)
     . + g(n5,n6)*P(n1)*P(n3)*K(n4)*L(n2)
     . + g(n5,n6)*P(n1)*P(n4)*K(n2)*L(n3)
     . + g(n5,n6)*P(n1)*P(n4)*K(n3)*L(n2)
     . + g(n5,n6)*P(n2)*P(n3)*K(n1)*L(n4)
     . + g(n5,n6)*P(n2)*P(n3)*K(n4)*L(n1)
     . + g(n5,n6)*P(n2)*P(n4)*K(n1)*L(n3)
     . + g(n5,n6)*P(n2)*P(n4)*K(n3)*L(n1)
     . + g(n5,n6)*P(n3)*P(n4)*K(n1)*L(n2)
     . + g(n5,n6)*P(n3)*P(n4)*K(n2)*L(n1)
      end

      function pvSPPPPPK(n1,n2,n3,n4,n5,n6,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPPPK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      pvSPPPPPK=
     . + P(n1)*P(n2)*P(n3)*P(n4)*P(n5)*K(n6)
     . + P(n1)*P(n2)*P(n3)*P(n4)*P(n6)*K(n5)
     . + P(n1)*P(n2)*P(n3)*P(n5)*P(n6)*K(n4)
     . + P(n1)*P(n2)*P(n4)*P(n5)*P(n6)*K(n3)
     . + P(n1)*P(n3)*P(n4)*P(n5)*P(n6)*K(n2)
     . + P(n2)*P(n3)*P(n4)*P(n5)*P(n6)*K(n1)
      end

      function pvSPPPPKL(n1,n2,n3,n4,n5,n6,P,K,L)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPPKL
      real(dp):: P(4),K(4),L(4)
      integer n1,n2,n3,n4,n5,n6
      pvSPPPPKL=
     . + P(n1)*P(n2)*P(n3)*P(n4)*K(n5)*L(n6)
     . + P(n1)*P(n2)*P(n3)*P(n4)*K(n6)*L(n5)
     . + P(n1)*P(n2)*P(n3)*P(n5)*K(n4)*L(n6)
     . + P(n1)*P(n2)*P(n3)*P(n5)*K(n6)*L(n4)
     . + P(n1)*P(n2)*P(n3)*P(n6)*K(n4)*L(n5)
     . + P(n1)*P(n2)*P(n3)*P(n6)*K(n5)*L(n4)
     . + P(n1)*P(n2)*P(n4)*P(n5)*K(n3)*L(n6)
     . + P(n1)*P(n2)*P(n4)*P(n5)*K(n6)*L(n3)
     . + P(n1)*P(n2)*P(n4)*P(n6)*K(n3)*L(n5)
     . + P(n1)*P(n2)*P(n4)*P(n6)*K(n5)*L(n3)
     . + P(n1)*P(n2)*P(n5)*P(n6)*K(n3)*L(n4)
     . + P(n1)*P(n2)*P(n5)*P(n6)*K(n4)*L(n3)
     . + P(n1)*P(n3)*P(n4)*P(n5)*K(n2)*L(n6)
     . + P(n1)*P(n3)*P(n4)*P(n5)*K(n6)*L(n2)
     . + P(n1)*P(n3)*P(n4)*P(n6)*K(n2)*L(n5)
     . + P(n1)*P(n3)*P(n4)*P(n6)*K(n5)*L(n2)
     . + P(n1)*P(n3)*P(n5)*P(n6)*K(n2)*L(n4)
     . + P(n1)*P(n3)*P(n5)*P(n6)*K(n4)*L(n2)
     . + P(n1)*P(n4)*P(n5)*P(n6)*K(n2)*L(n3)
     . + P(n1)*P(n4)*P(n5)*P(n6)*K(n3)*L(n2)
     . + P(n2)*P(n3)*P(n4)*P(n5)*K(n1)*L(n6)
     . + P(n2)*P(n3)*P(n4)*P(n5)*K(n6)*L(n1)
     . + P(n2)*P(n3)*P(n4)*P(n6)*K(n1)*L(n5)
     . + P(n2)*P(n3)*P(n4)*P(n6)*K(n5)*L(n1)
     . + P(n2)*P(n3)*P(n5)*P(n6)*K(n1)*L(n4)
     . + P(n2)*P(n3)*P(n5)*P(n6)*K(n4)*L(n1)
     . + P(n2)*P(n4)*P(n5)*P(n6)*K(n1)*L(n3)
     . + P(n2)*P(n4)*P(n5)*P(n6)*K(n3)*L(n1)
     . + P(n3)*P(n4)*P(n5)*P(n6)*K(n1)*L(n2)
     . + P(n3)*P(n4)*P(n5)*P(n6)*K(n2)*L(n1)
      end

      function pvSPPPKKK(n1,n2,n3,n4,n5,n6,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSPPPKKK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      pvSPPPKKK=
     . + P(n1)*P(n2)*P(n3)*K(n4)*K(n5)*K(n6)
     . + P(n1)*P(n2)*P(n4)*K(n3)*K(n5)*K(n6)
     . + P(n1)*P(n2)*P(n5)*K(n3)*K(n4)*K(n6)
     . + P(n1)*P(n2)*P(n6)*K(n3)*K(n4)*K(n5)
     . + P(n1)*P(n3)*P(n4)*K(n2)*K(n5)*K(n6)
     . + P(n1)*P(n3)*P(n5)*K(n2)*K(n4)*K(n6)
     . + P(n1)*P(n3)*P(n6)*K(n2)*K(n4)*K(n5)
     . + P(n1)*P(n4)*P(n5)*K(n2)*K(n3)*K(n6)
     . + P(n1)*P(n4)*P(n6)*K(n2)*K(n3)*K(n5)
     . + P(n1)*P(n5)*P(n6)*K(n2)*K(n3)*K(n4)
     . + P(n2)*P(n3)*P(n4)*K(n1)*K(n5)*K(n6)
     . + P(n2)*P(n3)*P(n5)*K(n1)*K(n4)*K(n6)
     . + P(n2)*P(n3)*P(n6)*K(n1)*K(n4)*K(n5)
     . + P(n2)*P(n4)*P(n5)*K(n1)*K(n3)*K(n6)
     . + P(n2)*P(n4)*P(n6)*K(n1)*K(n3)*K(n5)
     . + P(n2)*P(n5)*P(n6)*K(n1)*K(n3)*K(n4)
     . + P(n3)*P(n4)*P(n5)*K(n1)*K(n2)*K(n6)
     . + P(n3)*P(n4)*P(n6)*K(n1)*K(n2)*K(n5)
     . + P(n3)*P(n5)*P(n6)*K(n1)*K(n2)*K(n4)
     . + P(n4)*P(n5)*P(n6)*K(n1)*K(n2)*K(n3)
      end

      function pvSDDDDDD(n1,n2,n3,n4,n5,n6)
      implicit none
      include 'types.f'
      real(dp):: pvSDDDDDD
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDDDDD=
     . + g(n1,n2)*g(n3,n4)*g(n5,n6)
     . + g(n1,n2)*g(n3,n5)*g(n4,n6)
     . + g(n1,n2)*g(n3,n6)*g(n4,n5)
     . + g(n1,n3)*g(n2,n4)*g(n5,n6)
     . + g(n1,n3)*g(n2,n5)*g(n4,n6)
     . + g(n1,n3)*g(n2,n6)*g(n4,n5)
     . + g(n1,n4)*g(n2,n3)*g(n5,n6)
     . + g(n1,n4)*g(n2,n5)*g(n3,n6)
     . + g(n1,n4)*g(n2,n6)*g(n3,n5)
     . + g(n1,n5)*g(n2,n3)*g(n4,n6)
     . + g(n1,n5)*g(n2,n4)*g(n3,n6)
     . + g(n1,n5)*g(n2,n6)*g(n3,n4)
     . + g(n1,n6)*g(n2,n3)*g(n4,n5)
     . + g(n1,n6)*g(n2,n4)*g(n3,n5)
     . + g(n1,n6)*g(n2,n5)*g(n3,n4)
      end

      function pvSDDDDPP(n1,n2,n3,n4,n5,n6,P)
      implicit none
      include 'types.f'
      real(dp):: pvSDDDDPP
      real(dp):: P(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDDDPP=
     . + g(n1,n2)*g(n3,n4)*P(n5)*P(n6)
     . + g(n1,n2)*g(n3,n5)*P(n4)*P(n6)
     . + g(n1,n2)*g(n3,n6)*P(n4)*P(n5)
     . + g(n1,n2)*g(n4,n5)*P(n3)*P(n6)
     . + g(n1,n2)*g(n4,n6)*P(n3)*P(n5)
     . + g(n1,n2)*g(n5,n6)*P(n3)*P(n4)
     . + g(n1,n3)*g(n2,n4)*P(n5)*P(n6)
     . + g(n1,n3)*g(n2,n5)*P(n4)*P(n6)
     . + g(n1,n3)*g(n2,n6)*P(n4)*P(n5)
     . + g(n1,n3)*g(n4,n5)*P(n2)*P(n6)
     . + g(n1,n3)*g(n4,n6)*P(n2)*P(n5)
     . + g(n1,n3)*g(n5,n6)*P(n2)*P(n4)
     . + g(n1,n4)*g(n2,n3)*P(n5)*P(n6)
     . + g(n1,n4)*g(n2,n5)*P(n3)*P(n6)
     . + g(n1,n4)*g(n2,n6)*P(n3)*P(n5)
     . + g(n1,n4)*g(n3,n5)*P(n2)*P(n6)
     . + g(n1,n4)*g(n3,n6)*P(n2)*P(n5)
     . + g(n1,n4)*g(n5,n6)*P(n2)*P(n3)
     . + g(n1,n5)*g(n2,n3)*P(n4)*P(n6)
     . + g(n1,n5)*g(n2,n4)*P(n3)*P(n6)
     . + g(n1,n5)*g(n2,n6)*P(n3)*P(n4)
     . + g(n1,n5)*g(n3,n4)*P(n2)*P(n6)
     . + g(n1,n5)*g(n3,n6)*P(n2)*P(n4)
     . + g(n1,n5)*g(n4,n6)*P(n2)*P(n3)
     . + g(n1,n6)*g(n2,n3)*P(n4)*P(n5)
     . + g(n1,n6)*g(n2,n4)*P(n3)*P(n5)
     . + g(n1,n6)*g(n2,n5)*P(n3)*P(n4)
     . + g(n1,n6)*g(n3,n4)*P(n2)*P(n5)
     . + g(n1,n6)*g(n3,n5)*P(n2)*P(n4)
     . + g(n1,n6)*g(n4,n5)*P(n2)*P(n3)
      pvSDDDDPP=pvSDDDDPP
     . + g(n2,n3)*g(n4,n5)*P(n1)*P(n6)
     . + g(n2,n3)*g(n4,n6)*P(n1)*P(n5)
     . + g(n2,n3)*g(n5,n6)*P(n1)*P(n4)
     . + g(n2,n4)*g(n3,n5)*P(n1)*P(n6)
     . + g(n2,n4)*g(n3,n6)*P(n1)*P(n5)
     . + g(n2,n4)*g(n5,n6)*P(n1)*P(n3)
     . + g(n2,n5)*g(n3,n4)*P(n1)*P(n6)
     . + g(n2,n5)*g(n3,n6)*P(n1)*P(n4)
     . + g(n2,n5)*g(n4,n6)*P(n1)*P(n3)
     . + g(n2,n6)*g(n3,n4)*P(n1)*P(n5)
     . + g(n2,n6)*g(n3,n5)*P(n1)*P(n4)
     . + g(n2,n6)*g(n4,n5)*P(n1)*P(n3)
     . + g(n3,n4)*g(n5,n6)*P(n1)*P(n2)
     . + g(n3,n5)*g(n4,n6)*P(n1)*P(n2)
     . + g(n3,n6)*g(n4,n5)*P(n1)*P(n2)
      end

      function pvSDDPKKK(n1,n2,n3,n4,n5,n6,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPKKK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDPKKK=
     . + g(n1,n2)*P(n3)*K(n4)*K(n5)*K(n6)
     . + g(n1,n2)*P(n4)*K(n3)*K(n5)*K(n6)
     . + g(n1,n2)*P(n5)*K(n3)*K(n4)*K(n6)
     . + g(n1,n2)*P(n6)*K(n3)*K(n4)*K(n5)
     . + g(n1,n3)*P(n2)*K(n4)*K(n5)*K(n6)
     . + g(n1,n3)*P(n4)*K(n2)*K(n5)*K(n6)
     . + g(n1,n3)*P(n5)*K(n2)*K(n4)*K(n6)
     . + g(n1,n3)*P(n6)*K(n2)*K(n4)*K(n5)
     . + g(n1,n4)*P(n2)*K(n3)*K(n5)*K(n6)
     . + g(n1,n4)*P(n3)*K(n2)*K(n5)*K(n6)
     . + g(n1,n4)*P(n5)*K(n2)*K(n3)*K(n6)
     . + g(n1,n4)*P(n6)*K(n2)*K(n3)*K(n5)
     . + g(n1,n5)*P(n2)*K(n3)*K(n4)*K(n6)
     . + g(n1,n5)*P(n3)*K(n2)*K(n4)*K(n6)
     . + g(n1,n5)*P(n4)*K(n2)*K(n3)*K(n6)
     . + g(n1,n5)*P(n6)*K(n2)*K(n3)*K(n4)
     . + g(n1,n6)*P(n2)*K(n3)*K(n4)*K(n5)
     . + g(n1,n6)*P(n3)*K(n2)*K(n4)*K(n5)
     . + g(n1,n6)*P(n4)*K(n2)*K(n3)*K(n5)
     . + g(n1,n6)*P(n5)*K(n2)*K(n3)*K(n4)
      pvSDDPKKK=pvSDDPKKK
     . + g(n2,n3)*P(n1)*K(n4)*K(n5)*K(n6)
     . + g(n2,n3)*P(n4)*K(n1)*K(n5)*K(n6)
     . + g(n2,n3)*P(n5)*K(n1)*K(n4)*K(n6)
     . + g(n2,n3)*P(n6)*K(n1)*K(n4)*K(n5)
     . + g(n2,n4)*P(n1)*K(n3)*K(n5)*K(n6)
     . + g(n2,n4)*P(n3)*K(n1)*K(n5)*K(n6)
     . + g(n2,n4)*P(n5)*K(n1)*K(n3)*K(n6)
     . + g(n2,n4)*P(n6)*K(n1)*K(n3)*K(n5)
     . + g(n2,n5)*P(n1)*K(n3)*K(n4)*K(n6)
     . + g(n2,n5)*P(n3)*K(n1)*K(n4)*K(n6)
     . + g(n2,n5)*P(n4)*K(n1)*K(n3)*K(n6)
     . + g(n2,n5)*P(n6)*K(n1)*K(n3)*K(n4)
     . + g(n2,n6)*P(n1)*K(n3)*K(n4)*K(n5)
     . + g(n2,n6)*P(n3)*K(n1)*K(n4)*K(n5)
     . + g(n2,n6)*P(n4)*K(n1)*K(n3)*K(n5)
     . + g(n2,n6)*P(n5)*K(n1)*K(n3)*K(n4)
     . + g(n3,n4)*P(n1)*K(n2)*K(n5)*K(n6)
     . + g(n3,n4)*P(n2)*K(n1)*K(n5)*K(n6)
     . + g(n3,n4)*P(n5)*K(n1)*K(n2)*K(n6)
     . + g(n3,n4)*P(n6)*K(n1)*K(n2)*K(n5)
     . + g(n3,n5)*P(n1)*K(n2)*K(n4)*K(n6)
     . + g(n3,n5)*P(n2)*K(n1)*K(n4)*K(n6)
     . + g(n3,n5)*P(n4)*K(n1)*K(n2)*K(n6)
     . + g(n3,n5)*P(n6)*K(n1)*K(n2)*K(n4)
     . + g(n3,n6)*P(n1)*K(n2)*K(n4)*K(n5)
     . + g(n3,n6)*P(n2)*K(n1)*K(n4)*K(n5)
     . + g(n3,n6)*P(n4)*K(n1)*K(n2)*K(n5)
     . + g(n3,n6)*P(n5)*K(n1)*K(n2)*K(n4)
     . + g(n4,n5)*P(n1)*K(n2)*K(n3)*K(n6)
     . + g(n4,n5)*P(n2)*K(n1)*K(n3)*K(n6)
     . + g(n4,n5)*P(n3)*K(n1)*K(n2)*K(n6)
     . + g(n4,n5)*P(n6)*K(n1)*K(n2)*K(n3)
     . + g(n4,n6)*P(n1)*K(n2)*K(n3)*K(n5)
     . + g(n4,n6)*P(n2)*K(n1)*K(n3)*K(n5)
     . + g(n4,n6)*P(n3)*K(n1)*K(n2)*K(n5)
     . + g(n4,n6)*P(n5)*K(n1)*K(n2)*K(n3)
     . + g(n5,n6)*P(n1)*K(n2)*K(n3)*K(n4)
     . + g(n5,n6)*P(n2)*K(n1)*K(n3)*K(n4)
     . + g(n5,n6)*P(n3)*K(n1)*K(n2)*K(n4)
     . + g(n5,n6)*P(n4)*K(n1)*K(n2)*K(n3)
      end

      function pvSDDDDPK(n1,n2,n3,n4,n5,n6,P,K)
      implicit none
      include 'types.f'
      real(dp):: pvSDDDDPK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDDDPK=
     . + g(n1,n2)*g(n3,n4)*P(n5)*K(n6)
     . + g(n1,n2)*g(n3,n4)*P(n6)*K(n5)
     . + g(n1,n2)*g(n3,n5)*P(n4)*K(n6)
     . + g(n1,n2)*g(n3,n5)*P(n6)*K(n4)
     . + g(n1,n2)*g(n3,n6)*P(n4)*K(n5)
     . + g(n1,n2)*g(n3,n6)*P(n5)*K(n4)
     . + g(n1,n2)*g(n4,n5)*P(n3)*K(n6)
     . + g(n1,n2)*g(n4,n5)*P(n6)*K(n3)
     . + g(n1,n2)*g(n4,n6)*P(n3)*K(n5)
     . + g(n1,n2)*g(n4,n6)*P(n5)*K(n3)
     . + g(n1,n2)*g(n5,n6)*P(n3)*K(n4)
     . + g(n1,n2)*g(n5,n6)*P(n4)*K(n3)
     . + g(n1,n3)*g(n2,n4)*P(n5)*K(n6)
     . + g(n1,n3)*g(n2,n4)*P(n6)*K(n5)
     . + g(n1,n3)*g(n2,n5)*P(n4)*K(n6)
     . + g(n1,n3)*g(n2,n5)*P(n6)*K(n4)
     . + g(n1,n3)*g(n2,n6)*P(n4)*K(n5)
     . + g(n1,n3)*g(n2,n6)*P(n5)*K(n4)
     . + g(n1,n3)*g(n4,n5)*P(n2)*K(n6)
     . + g(n1,n3)*g(n4,n5)*P(n6)*K(n2)
     . + g(n1,n3)*g(n4,n6)*P(n2)*K(n5)
     . + g(n1,n3)*g(n4,n6)*P(n5)*K(n2)
     . + g(n1,n3)*g(n5,n6)*P(n2)*K(n4)
     . + g(n1,n3)*g(n5,n6)*P(n4)*K(n2)
     . + g(n1,n4)*g(n2,n3)*P(n5)*K(n6)
     . + g(n1,n4)*g(n2,n3)*P(n6)*K(n5)
     . + g(n1,n4)*g(n2,n5)*P(n3)*K(n6)
     . + g(n1,n4)*g(n2,n5)*P(n6)*K(n3)
     . + g(n1,n4)*g(n2,n6)*P(n3)*K(n5)
     . + g(n1,n4)*g(n2,n6)*P(n5)*K(n3)
     . + g(n1,n4)*g(n3,n5)*P(n2)*K(n6)
     . + g(n1,n4)*g(n3,n5)*P(n6)*K(n2)
     . + g(n1,n4)*g(n3,n6)*P(n2)*K(n5)
     . + g(n1,n4)*g(n3,n6)*P(n5)*K(n2)
     . + g(n1,n4)*g(n5,n6)*P(n2)*K(n3)
     . + g(n1,n4)*g(n5,n6)*P(n3)*K(n2)
     . + g(n1,n5)*g(n2,n3)*P(n4)*K(n6)
     . + g(n1,n5)*g(n2,n3)*P(n6)*K(n4)
     . + g(n1,n5)*g(n2,n4)*P(n3)*K(n6)
     . + g(n1,n5)*g(n2,n4)*P(n6)*K(n3)
     . + g(n1,n5)*g(n2,n6)*P(n3)*K(n4)
     . + g(n1,n5)*g(n2,n6)*P(n4)*K(n3)
     . + g(n1,n5)*g(n3,n4)*P(n2)*K(n6)
     . + g(n1,n5)*g(n3,n4)*P(n6)*K(n2)
     . + g(n1,n5)*g(n3,n6)*P(n2)*K(n4)
     . + g(n1,n5)*g(n3,n6)*P(n4)*K(n2)
     . + g(n1,n5)*g(n4,n6)*P(n2)*K(n3)
     . + g(n1,n5)*g(n4,n6)*P(n3)*K(n2)
     . + g(n1,n6)*g(n2,n3)*P(n4)*K(n5)
     . + g(n1,n6)*g(n2,n3)*P(n5)*K(n4)
     . + g(n1,n6)*g(n2,n4)*P(n3)*K(n5)
     . + g(n1,n6)*g(n2,n4)*P(n5)*K(n3)
     . + g(n1,n6)*g(n2,n5)*P(n3)*K(n4)
     . + g(n1,n6)*g(n2,n5)*P(n4)*K(n3)
     . + g(n1,n6)*g(n3,n4)*P(n2)*K(n5)
     . + g(n1,n6)*g(n3,n4)*P(n5)*K(n2)
     . + g(n1,n6)*g(n3,n5)*P(n2)*K(n4)
     . + g(n1,n6)*g(n3,n5)*P(n4)*K(n2)
     . + g(n1,n6)*g(n4,n5)*P(n2)*K(n3)
     . + g(n1,n6)*g(n4,n5)*P(n3)*K(n2)
      pvSDDDDPK=pvSDDDDPK
     . + g(n2,n3)*g(n4,n5)*P(n1)*K(n6)
     . + g(n2,n3)*g(n4,n5)*P(n6)*K(n1)
     . + g(n2,n3)*g(n4,n6)*P(n1)*K(n5)
     . + g(n2,n3)*g(n4,n6)*P(n5)*K(n1)
     . + g(n2,n3)*g(n5,n6)*P(n1)*K(n4)
     . + g(n2,n3)*g(n5,n6)*P(n4)*K(n1)
     . + g(n2,n4)*g(n3,n5)*P(n1)*K(n6)
     . + g(n2,n4)*g(n3,n5)*P(n6)*K(n1)
     . + g(n2,n4)*g(n3,n6)*P(n1)*K(n5)
     . + g(n2,n4)*g(n3,n6)*P(n5)*K(n1)
     . + g(n2,n4)*g(n5,n6)*P(n1)*K(n3)
     . + g(n2,n4)*g(n5,n6)*P(n3)*K(n1)
     . + g(n2,n5)*g(n3,n4)*P(n1)*K(n6)
     . + g(n2,n5)*g(n3,n4)*P(n6)*K(n1)
     . + g(n2,n5)*g(n3,n6)*P(n1)*K(n4)
     . + g(n2,n5)*g(n3,n6)*P(n4)*K(n1)
     . + g(n2,n5)*g(n4,n6)*P(n1)*K(n3)
     . + g(n2,n5)*g(n4,n6)*P(n3)*K(n1)
     . + g(n2,n6)*g(n3,n4)*P(n1)*K(n5)
     . + g(n2,n6)*g(n3,n4)*P(n5)*K(n1)
     . + g(n2,n6)*g(n3,n5)*P(n1)*K(n4)
     . + g(n2,n6)*g(n3,n5)*P(n4)*K(n1)
     . + g(n2,n6)*g(n4,n5)*P(n1)*K(n3)
     . + g(n2,n6)*g(n4,n5)*P(n3)*K(n1)
      pvSDDDDPK=pvSDDDDPK
     . + g(n3,n4)*g(n5,n6)*P(n1)*K(n2)
     . + g(n3,n4)*g(n5,n6)*P(n2)*K(n1)
     . + g(n3,n5)*g(n4,n6)*P(n1)*K(n2)
     . + g(n3,n5)*g(n4,n6)*P(n2)*K(n1)
     . + g(n3,n6)*g(n4,n5)*P(n1)*K(n2)
     . + g(n3,n6)*g(n4,n5)*P(n2)*K(n1)
      end

      function pvSDDPPPP(n1,n2,n3,n4,n5,n6,P)
      implicit none
      include 'types.f'
      real(dp):: pvSDDPPPP
      real(dp):: P(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDPPPP=
     . + g(n1,n2)*P(n3)*P(n4)*P(n5)*P(n6)
     . + g(n1,n3)*P(n2)*P(n4)*P(n5)*P(n6)
     . + g(n1,n4)*P(n2)*P(n3)*P(n5)*P(n6)
     . + g(n1,n5)*P(n2)*P(n3)*P(n4)*P(n6)
     . + g(n1,n6)*P(n2)*P(n3)*P(n4)*P(n5)
     . + g(n2,n3)*P(n1)*P(n4)*P(n5)*P(n6)
     . + g(n2,n4)*P(n1)*P(n3)*P(n5)*P(n6)
     . + g(n2,n5)*P(n1)*P(n3)*P(n4)*P(n6)
     . + g(n2,n6)*P(n1)*P(n3)*P(n4)*P(n5)
     . + g(n3,n4)*P(n1)*P(n2)*P(n5)*P(n6)
     . + g(n3,n5)*P(n1)*P(n2)*P(n4)*P(n6)
     . + g(n3,n6)*P(n1)*P(n2)*P(n4)*P(n5)
     . + g(n4,n5)*P(n1)*P(n2)*P(n3)*P(n6)
     . + g(n4,n6)*P(n1)*P(n2)*P(n3)*P(n5)
     . + g(n5,n6)*P(n1)*P(n2)*P(n3)*P(n4)
      end
    
      function pvSDDPPPK(n1,n2,n3,n4,n5,n6,P,K) 
      implicit none
      include 'types.f'
      real(dp):: pvSDDPPPK
      real(dp):: P(4),K(4)
      integer n1,n2,n3,n4,n5,n6
      include 'TRmetric.f'
      pvSDDPPPK=
     . + g(n1,n2)*P(n3)*P(n4)*P(n5)*K(n6)
     . + g(n1,n2)*P(n3)*P(n4)*P(n6)*K(n5)
     . + g(n1,n2)*P(n3)*P(n5)*P(n6)*K(n4)
     . + g(n1,n2)*P(n4)*P(n5)*P(n6)*K(n3)
     . + g(n1,n3)*P(n2)*P(n4)*P(n5)*K(n6)
     . + g(n1,n3)*P(n2)*P(n4)*P(n6)*K(n5)
     . + g(n1,n3)*P(n2)*P(n5)*P(n6)*K(n4)
     . + g(n1,n3)*P(n4)*P(n5)*P(n6)*K(n2)
     . + g(n1,n4)*P(n2)*P(n3)*P(n5)*K(n6)
     . + g(n1,n4)*P(n2)*P(n3)*P(n6)*K(n5)
     . + g(n1,n4)*P(n2)*P(n5)*P(n6)*K(n3)
     . + g(n1,n4)*P(n3)*P(n5)*P(n6)*K(n2)
     . + g(n1,n5)*P(n2)*P(n3)*P(n4)*K(n6)
     . + g(n1,n5)*P(n2)*P(n3)*P(n6)*K(n4)
     . + g(n1,n5)*P(n2)*P(n4)*P(n6)*K(n3)
     . + g(n1,n5)*P(n3)*P(n4)*P(n6)*K(n2)
     . + g(n1,n6)*P(n2)*P(n3)*P(n4)*K(n5)
     . + g(n1,n6)*P(n2)*P(n3)*P(n5)*K(n4)
     . + g(n1,n6)*P(n2)*P(n4)*P(n5)*K(n3)
     . + g(n1,n6)*P(n3)*P(n4)*P(n5)*K(n2)
      pvSDDPPPK=pvSDDPPPK
     . + g(n2,n3)*P(n1)*P(n4)*P(n5)*K(n6)
     . + g(n2,n3)*P(n1)*P(n4)*P(n6)*K(n5)
     . + g(n2,n3)*P(n1)*P(n5)*P(n6)*K(n4)
     . + g(n2,n3)*P(n4)*P(n5)*P(n6)*K(n1)
     . + g(n2,n4)*P(n1)*P(n3)*P(n5)*K(n6)
     . + g(n2,n4)*P(n1)*P(n3)*P(n6)*K(n5)
     . + g(n2,n4)*P(n1)*P(n5)*P(n6)*K(n3)
     . + g(n2,n4)*P(n3)*P(n5)*P(n6)*K(n1)
     . + g(n2,n5)*P(n1)*P(n3)*P(n4)*K(n6)
     . + g(n2,n5)*P(n1)*P(n3)*P(n6)*K(n4)
     . + g(n2,n5)*P(n1)*P(n4)*P(n6)*K(n3)
     . + g(n2,n5)*P(n3)*P(n4)*P(n6)*K(n1)
     . + g(n2,n6)*P(n1)*P(n3)*P(n4)*K(n5)
     . + g(n2,n6)*P(n1)*P(n3)*P(n5)*K(n4)
     . + g(n2,n6)*P(n1)*P(n4)*P(n5)*K(n3)
     . + g(n2,n6)*P(n3)*P(n4)*P(n5)*K(n1)
      pvSDDPPPK=pvSDDPPPK
     . + g(n3,n4)*P(n1)*P(n2)*P(n5)*K(n6)
     . + g(n3,n4)*P(n1)*P(n2)*P(n6)*K(n5)
     . + g(n3,n4)*P(n1)*P(n5)*P(n6)*K(n2)
     . + g(n3,n4)*P(n2)*P(n5)*P(n6)*K(n1)
     . + g(n3,n5)*P(n1)*P(n2)*P(n4)*K(n6)
     . + g(n3,n5)*P(n1)*P(n2)*P(n6)*K(n4)
     . + g(n3,n5)*P(n1)*P(n4)*P(n6)*K(n2)
     . + g(n3,n5)*P(n2)*P(n4)*P(n6)*K(n1)
     . + g(n3,n6)*P(n1)*P(n2)*P(n4)*K(n5)
     . + g(n3,n6)*P(n1)*P(n2)*P(n5)*K(n4)
     . + g(n3,n6)*P(n1)*P(n4)*P(n5)*K(n2)
     . + g(n3,n6)*P(n2)*P(n4)*P(n5)*K(n1)
     . + g(n4,n5)*P(n1)*P(n2)*P(n3)*K(n6)
     . + g(n4,n5)*P(n1)*P(n2)*P(n6)*K(n3)
     . + g(n4,n5)*P(n1)*P(n3)*P(n6)*K(n2)
     . + g(n4,n5)*P(n2)*P(n3)*P(n6)*K(n1)
     . + g(n4,n6)*P(n1)*P(n2)*P(n3)*K(n5)
     . + g(n4,n6)*P(n1)*P(n2)*P(n5)*K(n3)
     . + g(n4,n6)*P(n1)*P(n3)*P(n5)*K(n2)
     . + g(n4,n6)*P(n2)*P(n3)*P(n5)*K(n1)
     . + g(n5,n6)*P(n1)*P(n2)*P(n3)*K(n4)
     . + g(n5,n6)*P(n1)*P(n2)*P(n4)*K(n3)
     . + g(n5,n6)*P(n1)*P(n3)*P(n4)*K(n2)
     . + g(n5,n6)*P(n2)*P(n3)*P(n4)*K(n1)
      end


