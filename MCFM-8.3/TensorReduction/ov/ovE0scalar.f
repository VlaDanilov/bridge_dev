      subroutine ovE0scalar(E0,p1,p2,p3,p4,p5,p1p2,p2p3,p3p4,p4p5,p5p1,
     . m1s,m2s,m3s,m4s,m5s)
      implicit none
      include 'types.f'
      include 'TRscale.f'
      complex(dp):: E0(-2:0),Det,
     . D01(-2:0),D02(-2:0),D03(-2:0),D04(-2:0),D05(-2:0)
      real(dp):: p1,p2,p3,p4,p5,p1p2,p2p3,p3p4,p4p5,p5p1,
     . m1s,m2s,m3s,m4s,m5s
      integer i, j,ep
      complex(dp):: Y(5,5), Yi(5,5),detYi(5)
      complex(dp):: Yflat(25), Yiflat(25)
      equivalence (Y, Yflat)
      equivalence (Yi, Yiflat)

      complex(dp):: pvXDet,trI4
      external pvXDet
      include 'cplx.h'
c      Complex(dp):: Det5



      Y(1,1) = cplx1(2*m1s)
      Y(1,2) = cplx1(m1s + m2s - p1)
      Y(2,1) = Y(1,2)
      Y(1,3) = cplx1(m1s + m3s - p1p2)
      Y(3,1) = Y(1,3)
      Y(1,4) = cplx1(m1s + m4s - p4p5)
      Y(4,1) = Y(1,4)
      Y(1,5) = cplx1(m1s + m5s - p5)
      Y(5,1) = Y(1,5)
      Y(2,2) = cplx1(2*m2s)
      Y(2,3) = cplx1(m2s + m3s - p2)
      Y(3,2) = Y(2,3)
      Y(2,4) = cplx1(m2s + m4s - p2p3)
      Y(4,2) = Y(2,4)
      Y(2,5) = cplx1(m2s + m5s - p5p1)
      Y(5,2) = Y(2,5)
      Y(3,3) = cplx1(2*m3s)
      Y(3,4) = cplx1(m3s + m4s - p3)
      Y(4,3) = Y(3,4)
      Y(3,5) = cplx1(m3s + m5s - p3p4)
      Y(5,3) = Y(3,5)
      Y(4,4) = cplx1(2*m4s)
      Y(4,5) = cplx1(m4s + m5s - p4)
      Y(5,4) = Y(4,5)
      Y(5,5) = cplx1(2*m5s)
      
      do i = 1, 5
        do j = 1, 25
          Yiflat(j) = Yflat(j)
        enddo
        do j = 1, 5
          Yi(j,i) = cplx2(1d0,0d0)
        enddo
        detYi(i) = pvXDet(Yi, 5)
      enddo

        do ep=-2,0
        D01(ep)=trI4(p2,p3,p4,p5p1,p2p3,p3p4,m2s,m3s,m4s,m5s,musq,ep)
        enddo

        do ep=-2,0
        D02(ep)=trI4(p1p2,p3,p4,p5,p4p5,p3p4,m1s,m3s,m4s,m5s,musq,ep)
        enddo
        do ep=-2,0
        D03(ep)=trI4(p1,p2p3,p4,p5,p4p5,p5p1,m1s,m2s,m4s,m5s,musq,ep)
        enddo
        do ep=-2,0
        D04(ep)=trI4(p1,p2,p3p4,p5,p1p2,p5p1,m1s,m2s,m3s,m5s,musq,ep)
        enddo
        do ep=-2,0
        D05(ep)=trI4(p1,p2,p3,p4p5,p1p2,p2p3,m1s,m2s,m3s,m4s,musq,ep)
        enddo

      Det=pvXDet(Y,5)

        do ep=-2,0
        E0(ep) = -(
     &    detYi(1)*D01(ep)+detYi(2)*D02(ep)+detYi(3)*D03(ep)
     &   +detYi(4)*D04(ep)+detYi(5)*D05(ep))/Det

        
        enddo
      end

