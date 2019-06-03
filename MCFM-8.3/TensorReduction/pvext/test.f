      program test
      include 'types.f'
      include 'pvCnames.f'
      include 'pvextCv.f'
      include 'TRydef.f'
      include 'pvverbose.f'
      integer r,rmaxC,rmaxD,rmaxE
      real(dp):: p1(4),p2(4),p3(4),p4(4),q1(4),q2(4),q3(4),q4(4),
     & m1s,m2s,m3s,m4s,m5s
      complex(dp):: FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     & FC3(y3max,-2:0),FC4(y4max,-2:0),FC5(y5max,-2:0),FC6(y6max,-2:0),
     & FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     & FD3(y3max,-2:0),FD4(y4max,-2:0),FD5(y5max,-2:0),FD6(y6max,-2:0),
     & FE0(-2:0),FE1(y1max,-2:0),FE2(y2max,-2:0),
     & FE3(y3max,-2:0),FE4(y4max,-2:0),FE5(y5max,-2:0)
      logical failed,dopvext

c--- flag to check pvext or not
      dopvext=.false.
      
      p1(:)=0._dp
      p1(4)=3._dp
      p1(3)=3._dp
      
      p2(:)=0._dp
      p2(4)=2._dp
      p2(3)=-2._dp
      
      p3(4)=3.5_dp
      p3(1)=1.1_dp
      p3(2)=0._dp
      p3(3)=1.3_dp
      
      p4(4)=4.5_dp
      p4(1)=-1.1_dp
      p4(2)=0.4_dp
      p4(3)=2.7_dp
      
      m1s=0.123_dp
      m2s=0.419_dp
      m3s=0.391_dp
      m4s=0.0273_dp
      m5s=0.0773_dp
      m1s=0._dp
      m2s=0._dp
      m3s=0._dp
      m4s=0._dp
      m5s=0._dp
      
      q1=p1
      q2=p1+p2
      q3=p1+p2+p3
      q4=p1+p2+p3+p4
      
      call qlinit
      call pvsetmudim(1._dp)
      call TRsetmaxindex(6,6,4)

      if (dopvext) then
        call pvextCtensor(p1,p2,m1s,m2s,m3s,FC0,FC1,FC2,FC3)
        call pvextDtensor(p1,p2,p3,m1s,m2s,m3s,m4s,FD0,FD1,FD2,FD3,FD4)
        rmaxC=3
        rmaxD=4
        rmaxE=0
      else
        call pvCtensor(q1,q2,m1s,m2s,m3s,FC0,FC1,FC2,FC3,FC4,FC5,FC6)
        call pvDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     &                 FD0,FD1,FD2,FD3,FD4,FD5,FD6)
        call pvEtensor(q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     &                 FE0,FE1,FE2,FE3,FE4,FE5)
        rmaxC=6
        rmaxD=6
        rmaxE=5
      endif

c--- triangles
      do r=1,rmaxC
        call pvCcheck(r,q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
        write(6,*) 'pvCcheck: rank,failed=',r,failed
      enddo
            

c--- boxes
      do r=1,rmaxD
        call pvDcheck(r,q1,q2,q3,m1s,m2s,m3s,m4s,
     & FD0,FD1,FD2,FD3,FD4,FD5,FD6,failed)
        write(6,*) 'pvDcheck: rank,failed=',r,failed
      enddo
      
c--- pentagons
      do r=1,rmaxE
        call pvEcheck(r,q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     & FE0,FE1,FE2,FE3,FE4,FE5,failed)
        write(6,*) 'pvEcheck: rank,failed=',r,failed
      enddo
      
      stop
      end
      
