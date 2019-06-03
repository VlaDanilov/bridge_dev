      program test
      include 'types.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'TRydef.f'
      include 'pvverbose.f'
      real(dp) p1(4),p2(4),p3(4),p4(4),q1(4),q2(4),q3(4),q4(4),
     & m1s,m2s,m3s,m4s
      complex(dp):: FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     & FC3(y3max,-2:0),FC4(y4max,-2:0),FC5(y5max,-2:0),FC6(y6max,-2:0),
     & FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     & FD3(y3max,-2:0),FD4(y4max,-2:0),FD5(y5max,-2:0),FD6(y6max,-2:0)
      logical failed
    
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
      p4(2)=0._dp
      p4(3)=2.7_dp
      
      m1s=0.123_dp
      m2s=0.419_dp
      m3s=0.391_dp
      m4s=0.0273_dp
      m5s=0.0773_dp
!      m1s=0._dp
!      m2s=0._dp
!      m3s=0._dp
!      m4s=0._dp
      
      q1=p1
      q2=p1+p2
      q3=p1+p2+p3
      q4=p1+p2+p3+p4
      
      call qlinit
      call pvsetmudim(1._dp)
      call TRsetmaxindex(6,6,5)

c--- triangles
      call pvCtensor(q1,q2,m1s,m2s,m3s,
     . FC0,FC1,FC2,FC3,FC4,FC5,FC6)
      call pvCcheck(6,q1,q2,m1s,m2s,m3s,
     & FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      write(6,*) 'pvCcheck: failed=',failed
            
c--- boxes
      call pvDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     . FD0,FD1,FD2,FD3,FD4,FD5,FD6)
      call pvDcheck(6,q1,q2,q3,m1s,m2s,m3s,m4s,
     & FD0,FD1,FD2,FD3,FD4,FD5,FD6,failed)
      write(6,*) 'pvDcheck: failed=',failed
      
      call pvEtensor(q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     &  FE0,FE1,FE2,FE3,FE4,FE5)
      call pvEcheck(4,q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     & FE0,FE1,FE2,FE3,FE4,FE5,failed)
      write(6,*) 'pvEcheck: failed=',failed
      stop
      end
      
