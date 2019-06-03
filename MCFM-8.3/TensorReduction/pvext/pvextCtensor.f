      subroutine pvextCtensor(p1,p2,m1s,m2s,m3s,FC0,FC1,FC2,FC3)
      implicit none
C     p1,p2 are the external momenta
C     m1s,m2s,m3s are the squares of the masses in the propagators
      include 'types.f'
      include 'pvCnames.f'
      include 'pvextCv.f'
      include 'TRydef.f'
      include 'TRmetric.f'
      complex(dp):: FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     & FC3(y3max,-2:0)
      real(dp):: p1Dp1,p2Dp2,p3Dp3,p1(4),p2(4),p3(4),m1s,m2s,m3s,
     & pvSPK,pvSPKK,pvSDDP
      integer::n1,n2,n3,C0i,pvextCcache
      logical,save:: first=.true.
!$omp threadprivate(first)

      if (first) then
      first=.false.
      call pvarraysetup
      endif
      
      p3(:)=-p1(:)-p2(:)
      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2

      C0i=pvextCcache(p1Dp1,p2Dp2,p3Dp3,m1s,m2s,m3s)
      
      FC0(:)=Cv(C0i+cc0,:)

      do n1=1,4
      FC1(n1,:)=Cv(C0i+cc1,:)*p1(n1)+Cv(C0i+cc2,:)*p2(n1)
      enddo


      do n1=1,4
      do n2=n1,4
      FC2(y2(n1,n2),:)=
     & +p1(n1)*p1(n2)*Cv(C0i+cc11,:)
     & +p2(n1)*p2(n2)*Cv(C0i+cc22,:)
     & +pvSPK(n1,n2,p1,p2)*Cv(C0i+cc12,:)
     & +g(n1,n2)*Cv(C0i+cc00,:)
      enddo
      enddo

      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      FC3(y3(n1,n2,n3),:)=
     & +p1(n1)*p1(n2)*p1(n3)*Cv(C0i+cc111,:)
     & +p2(n1)*p2(n2)*p2(n3)*Cv(C0i+cc222,:)
     & +pvSPKK(n1,n2,n3,p2,p1)*Cv(C0i+cc112,:)
     & +pvSPKK(n1,n2,n3,p1,p2)*Cv(C0i+cc122,:)
     & +pvSDDP(n1,n2,n3,p1)*Cv(C0i+cc001,:)
     & +pvSDDP(n1,n2,n3,p2)*Cv(C0i+cc002,:)
      enddo
      enddo
      enddo

      return
      end

