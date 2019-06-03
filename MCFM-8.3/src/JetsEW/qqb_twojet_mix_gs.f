      subroutine qqb_twojet_mix_gs(p,msqc)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c--- for mixed QCD/Z matrix element for dijet production,
c     q(-p1)+qbar(-p2) --> q(p3)+q~(p4) + g(p5)

c--- all momenta are incoming
c
      implicit none 
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd
c--- slightly obtuse notation, to simplify declaration lines
      real(dp):: p(mxpart,4),msqc(maxd,fn:nf,fn:nf)
      double precision
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq15_3(fn:nf,fn:nf),msq15_4(fn:nf,fn:nf),
     & msq35_1(fn:nf,fn:nf),msq45_1(fn:nf,fn:nf),
     & msq25_3(fn:nf,fn:nf),msq25_4(fn:nf,fn:nf),
     & msq35_2(fn:nf,fn:nf),msq45_2(fn:nf,fn:nf),
     & msq35_4(fn:nf,fn:nf),msq45_3(fn:nf,fn:nf),
     & msq13_2(fn:nf,fn:nf),msq24_1(fn:nf,fn:nf),
     & msq14_2(fn:nf,fn:nf),msq23_1(fn:nf,fn:nf),

     & msq15_2v(fn:nf,fn:nf),msq25_1v(fn:nf,fn:nf),
     & msq15_3v(fn:nf,fn:nf),msq15_4v(fn:nf,fn:nf),
     & msq35_1v(fn:nf,fn:nf),msq45_1v(fn:nf,fn:nf),
     & msq25_3v(fn:nf,fn:nf),msq25_4v(fn:nf,fn:nf),
     & msq35_2v(fn:nf,fn:nf),msq45_2v(fn:nf,fn:nf),
     & msq35_4v(fn:nf,fn:nf),msq45_3v(fn:nf,fn:nf),
     & msq13_2v(fn:nf,fn:nf),msq24_1v(fn:nf,fn:nf),
     & msq14_2v(fn:nf,fn:nf),msq23_1v(fn:nf,fn:nf),

     & sub15_2(4),sub25_1(4),
     & sub15_3(4),sub15_4(4),
     & sub35_1(4),sub45_1(4),
     & sub25_3(4),sub25_4(4),
     & sub35_2(4),sub45_2(4),
     & sub35_4(4),sub45_3(4),
     & sub13_2(4),sub24_1(4),
     & sub14_2(4),sub23_1(4),

     & sub15_2v,sub25_1v,
     & sub15_3v,sub15_4v,
     & sub35_1v,sub45_1v,
     & sub25_3v,sub25_4v,
     & sub35_2v,sub45_2v,
     & sub35_4v,sub45_3v,
     & sub13_2v,sub24_1v,
     & sub14_2v,sub23_1v

      double precision
     & msq_mix15_2(0:3,fn:nf,fn:nf),msq_mix25_1(0:3,fn:nf,fn:nf),
     & msq_mix15_3(0:3,fn:nf,fn:nf),msq_mix35_1(0:3,fn:nf,fn:nf),
     & msq_mix15_4(0:3,fn:nf,fn:nf),msq_mix45_1(0:3,fn:nf,fn:nf),
     & msq_mix25_3(0:3,fn:nf,fn:nf),msq_mix35_2(0:3,fn:nf,fn:nf),
     & msq_mix25_4(0:3,fn:nf,fn:nf),msq_mix45_2(0:3,fn:nf,fn:nf),
     & msq_mix35_4(0:3,fn:nf,fn:nf),msq_mix45_3(0:3,fn:nf,fn:nf),
     & msq_mix13_2(0:3,fn:nf,fn:nf),msq_mix24_1(0:3,fn:nf,fn:nf),
     & msq_mix14_2(0:3,fn:nf,fn:nf),msq_mix23_1(0:3,fn:nf,fn:nf)

      double precision
     & msq_mix15_2v(0:3,fn:nf,fn:nf),msq_mix25_1v(0:3,fn:nf,fn:nf),
     & msq_mix15_3v(0:3,fn:nf,fn:nf),msq_mix35_1v(0:3,fn:nf,fn:nf),
     & msq_mix15_4v(0:3,fn:nf,fn:nf),msq_mix45_1v(0:3,fn:nf,fn:nf),
     & msq_mix25_3v(0:3,fn:nf,fn:nf),msq_mix35_2v(0:3,fn:nf,fn:nf),
     & msq_mix25_4v(0:3,fn:nf,fn:nf),msq_mix45_2v(0:3,fn:nf,fn:nf),
     & msq_mix35_4v(0:3,fn:nf,fn:nf),msq_mix45_3v(0:3,fn:nf,fn:nf),
     & msq_mix13_2v(0:3,fn:nf,fn:nf),msq_mix24_1v(0:3,fn:nf,fn:nf),
     & msq_mix14_2v(0:3,fn:nf,fn:nf),msq_mix23_1v(0:3,fn:nf,fn:nf)

      real(dp):: bit(8)

      external qqb_twojet_mix,qqb_twojet_mix_gqsub,
     & qqb_twojet_mix_gqsub2,donothing_gvec

      ndmax=12

c-- initialize the matrix elements to zero
      msqc(:,:,:)=zip
      
c--initial-final/final-initial
      call dips(1,p,1,5,3,sub15_3,sub15_3v,msq15_3,msq15_3v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix15_3)
      call dips(1,p,3,5,1,sub35_1,sub35_1v,msq35_1,msq35_1v,
     . qqb_twojet_mix,donothing_gvec)
      
      call dips(2,p,2,5,3,sub25_3,sub25_3v,msq25_3,msq25_3v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix25_3)
      call dips(2,p,3,5,2,sub35_2,sub35_2v,msq35_2,msq35_2v,
     . qqb_twojet_mix,donothing_gvec)
      
      call dips(3,p,1,5,4,sub15_4,sub15_4v,msq15_4,msq15_4v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix15_4)
      call dips(3,p,4,5,1,sub45_1,sub45_1v,msq45_1,msq45_1v,
     . qqb_twojet_mix,donothing_gvec)

      call dips(4,p,2,5,4,sub25_4,sub25_4v,msq25_4,msq25_4v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix25_4)
      call dips(4,p,4,5,2,sub45_2,sub45_2v,msq45_2,msq45_2v,
     . qqb_twojet_mix,donothing_gvec)
 
c--initial-initial
      call dips(5,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix15_2)
      
      call dips(6,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix25_1)
      
c--final-final
      call dips(7,p,3,5,4,sub35_4,sub35_4v,msq35_4,msq35_4v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix35_4)

      call dips(8,p,4,5,3,sub45_3,sub45_3v,msq45_3,msq45_3v,
     . qqb_twojet_mix,donothing_gvec)
      call store_msq_mix(msq_mix45_3)
 
c--initial-initial for gluon-initiated pieces: use special gqsub
      call dips(9,p,1,3,2,sub13_2,sub13_2v,msq13_2,msq13_2v,
     . qqb_twojet_mix_gqsub,donothing_gvec)
      call store_msq_mix(msq_mix13_2)
            
      call dips(10,p,2,4,1,sub24_1,sub24_1v,msq24_1,msq24_1v,
     . qqb_twojet_mix_gqsub,donothing_gvec)
      call store_msq_mix(msq_mix24_1)
      
c--initial-initial for gluon-initiated pieces: use special gqsub2
      call dips(11,p,1,4,2,sub14_2,sub14_2v,msq14_2,msq14_2v,
     . qqb_twojet_mix_gqsub2,donothing_gvec)
      call store_msq_mix(msq_mix14_2)
            
      call dips(12,p,2,3,1,sub23_1,sub23_1v,msq23_1,msq23_1v,
     . qqb_twojet_mix_gqsub2,donothing_gvec)
      call store_msq_mix(msq_mix23_1)
      
c--- fill the dipole contributions
      do j=-nf,nf
      do k=-nf,nf

c--- quark-quark or antiquark-antiquark
      if (j*k .gt. 0) then

      bit(:)=0._dp
      bit(2)=bit(2)+msq_mix25_3(0,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)+msq_mix15_4(0,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(5)=bit(5)-msq_mix15_2(0,j,k)*sub15_2(qq)
      bit(6)=bit(6)-msq_mix25_1(0,j,k)*sub25_1(qq)
      bit(7)=bit(7)-msq_mix35_4(0,j,k)*sub35_4(qq)
      bit(8)=bit(8)-msq_mix45_3(0,j,k)*sub45_3(qq)

      bit(1)=bit(1)+two*Cf*msq_mix15_3(1,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)-one/xn*msq_mix25_3(1,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)-one/xn*msq_mix15_4(1,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)+two*Cf*msq_mix25_4(1,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)+one/xn*msq_mix15_2(1,j,k)*sub15_2(qq)
      bit(6)=bit(6)+one/xn*msq_mix25_1(1,j,k)*sub25_1(qq)
      bit(7)=bit(7)+one/xn*msq_mix35_4(1,j,k)*sub35_4(qq)
      bit(8)=bit(8)+one/xn*msq_mix45_3(1,j,k)*sub45_3(qq)

      bit(1)=bit(1)-one/xn*msq_mix15_3(2,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)+two*Cf*msq_mix25_3(2,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)+two*Cf*msq_mix15_4(2,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)-one/xn*msq_mix25_4(2,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)+one/xn*msq_mix15_2(2,j,k)*sub15_2(qq)
      bit(6)=bit(6)+one/xn*msq_mix25_1(2,j,k)*sub25_1(qq)
      bit(7)=bit(7)+one/xn*msq_mix35_4(2,j,k)*sub35_4(qq)
      bit(8)=bit(8)+one/xn*msq_mix45_3(2,j,k)*sub45_3(qq)

      bit(1)=bit(1)+msq_mix15_3(3,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(4)=bit(4)+msq_mix25_4(3,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)-msq_mix15_2(3,j,k)*sub15_2(qq)
      bit(6)=bit(6)-msq_mix25_1(3,j,k)*sub25_1(qq)
      bit(7)=bit(7)-msq_mix35_4(3,j,k)*sub35_4(qq)
      bit(8)=bit(8)-msq_mix45_3(3,j,k)*sub45_3(qq)

      msqc(1:8,j,k)=bit(:)
      
      endif

      if (j*k .lt. 0) then
c--- quark-antiquark or antiquark-quark

      if (j .eq. -k) then
c------ quarks of the same flavor
      bit(:)=0._dp
      bit(1)=bit(1)+msq_mix15_3(0,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)-msq_mix25_3(0,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)-msq_mix15_4(0,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)+msq_mix25_4(0,j,k)*(sub25_4(qq)+sub45_2(qq))

      bit(1)=bit(1)-one/xn*msq_mix15_3(1,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)+one/xn*msq_mix25_3(1,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)+one/xn*msq_mix15_4(1,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)-one/xn*msq_mix25_4(1,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)+two*Cf*msq_mix15_2(1,j,k)*sub15_2(qq)
      bit(6)=bit(6)+two*Cf*msq_mix25_1(1,j,k)*sub25_1(qq)
      bit(7)=bit(7)+two*Cf*msq_mix35_4(1,j,k)*sub35_4(qq)
      bit(8)=bit(8)+two*Cf*msq_mix45_3(1,j,k)*sub45_3(qq)

      bit(1)=bit(1)+two*Cf*msq_mix15_3(2,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)+one/xn*msq_mix25_3(2,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)+one/xn*msq_mix15_4(2,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)+two*Cf*msq_mix25_4(2,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)-one/xn*msq_mix15_2(2,j,k)*sub15_2(qq)
      bit(6)=bit(6)-one/xn*msq_mix25_1(2,j,k)*sub25_1(qq)
      bit(7)=bit(7)-one/xn*msq_mix35_4(2,j,k)*sub35_4(qq)
      bit(8)=bit(8)-one/xn*msq_mix45_3(2,j,k)*sub45_3(qq)

      bit(2)=bit(2)-msq_mix25_3(3,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)-msq_mix15_4(3,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(5)=bit(5)+msq_mix15_2(3,j,k)*sub15_2(qq)
      bit(6)=bit(6)+msq_mix25_1(3,j,k)*sub25_1(qq)
      bit(7)=bit(7)+msq_mix35_4(3,j,k)*sub35_4(qq)
      bit(8)=bit(8)+msq_mix45_3(3,j,k)*sub45_3(qq)

      else
      
      bit(:)=0._dp
      bit(2)=bit(2)-msq_mix25_3(0,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)-msq_mix15_4(0,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(5)=bit(5)+msq_mix15_2(0,j,k)*sub15_2(qq)
      bit(6)=bit(6)+msq_mix25_1(0,j,k)*sub25_1(qq)
      bit(7)=bit(7)+msq_mix35_4(0,j,k)*sub35_4(qq)
      bit(8)=bit(8)+msq_mix45_3(0,j,k)*sub45_3(qq)

      bit(1)=bit(1)+two*Cf*msq_mix15_3(1,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)+one/xn*msq_mix25_3(1,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)+one/xn*msq_mix15_4(1,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)+two*Cf*msq_mix25_4(1,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)-one/xn*msq_mix15_2(1,j,k)*sub15_2(qq)
      bit(6)=bit(6)-one/xn*msq_mix25_1(1,j,k)*sub25_1(qq)
      bit(7)=bit(7)-one/xn*msq_mix35_4(1,j,k)*sub35_4(qq)
      bit(8)=bit(8)-one/xn*msq_mix45_3(1,j,k)*sub45_3(qq)

      bit(1)=bit(1)-one/xn*msq_mix15_3(2,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)+one/xn*msq_mix25_3(2,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)+one/xn*msq_mix15_4(2,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)-one/xn*msq_mix25_4(2,j,k)*(sub25_4(qq)+sub45_2(qq))
      bit(5)=bit(5)+two*Cf*msq_mix15_2(2,j,k)*sub15_2(qq)
      bit(6)=bit(6)+two*Cf*msq_mix25_1(2,j,k)*sub25_1(qq)
      bit(7)=bit(7)+two*Cf*msq_mix35_4(2,j,k)*sub35_4(qq)
      bit(8)=bit(8)+two*Cf*msq_mix45_3(2,j,k)*sub45_3(qq)

      bit(1)=bit(1)+msq_mix15_3(3,j,k)*(sub15_3(qq)+sub35_1(qq))
      bit(2)=bit(2)-msq_mix25_3(3,j,k)*(sub25_3(qq)+sub35_2(qq))
      bit(3)=bit(3)-msq_mix15_4(3,j,k)*(sub15_4(qq)+sub45_1(qq))
      bit(4)=bit(4)+msq_mix25_4(3,j,k)*(sub25_4(qq)+sub45_2(qq))

      endif
      
      msqc(1:8,j,k)=bit(:)
      
      endif

c--- quark-gluon cases: only trivial collinear singularities
      if ((j .lt. 0) .and. (k .eq. 0)) then
         msq25_1(:,:)=msq_mix25_1(1,:,:)+msq_mix25_1(2,:,:)
         msq24_1(:,:)=msq_mix24_1(1,:,:)+msq_mix24_1(2,:,:)
         msq23_1(:,:)=msq_mix23_1(1,:,:)+msq_mix23_1(2,:,:)
         msqc(6,j,k)=sub25_1(qg)*(
     &     +msq25_1(j,-1)+msq25_1(j,-2)
     &     +msq25_1(j,-3)+msq25_1(j,-4))
         msqc(10,j,k)=sub24_1(qg)*(
     &      msq24_1(j,+1)+msq24_1(j,+2)
     &     +msq24_1(j,+3)+msq24_1(j,+4)
     &     -msq24_1(j,-j)+0.5_dp*msq24_1(j,-j))
         msqc(12,j,k)=sub23_1(qg)*msq23_1(j,-j)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
         msq25_1(:,:)=msq_mix25_1(1,:,:)+msq_mix25_1(2,:,:)
         msq24_1(:,:)=msq_mix24_1(1,:,:)+msq_mix24_1(2,:,:)
         msq23_1(:,:)=msq_mix23_1(1,:,:)+msq_mix23_1(2,:,:)
         msqc(6,j,k)=sub25_1(qg)*(
     &      msq25_1(j,+1)+msq25_1(j,+2)
     &     +msq25_1(j,+3)+msq25_1(j,+4))
         msqc(10,j,k)=sub24_1(qg)*(
     &      msq24_1(j,-1)+msq24_1(j,-2)
     &     +msq24_1(j,-3)+msq24_1(j,-4)
     &     -msq24_1(j,-j)+0.5_dp*msq24_1(j,-j))
         msqc(12,j,k)=sub23_1(qg)*msq23_1(j,-j)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
         msq15_2(:,:)=msq_mix15_2(1,:,:)+msq_mix15_2(2,:,:)
         msq13_2(:,:)=msq_mix13_2(1,:,:)+msq_mix13_2(2,:,:)
         msq14_2(:,:)=msq_mix14_2(1,:,:)+msq_mix14_2(2,:,:)
         msqc(5,j,k)=sub15_2(qg)*(
     &      msq15_2(-1,k)+msq15_2(-2,k)
     &     +msq15_2(-3,k)+msq15_2(-4,k))
         msqc(9,j,k)=sub13_2(qg)*(
     &      msq13_2(+1,k)+msq13_2(+2,k)
     &     +msq13_2(+3,k)+msq13_2(+4,k)
     &     -msq13_2(-k,k)+0.5_dp*msq13_2(-k,k))
         msqc(11,j,k)=sub14_2(qg)*msq14_2(-k,k)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
         msq15_2(:,:)=msq_mix15_2(1,:,:)+msq_mix15_2(2,:,:)
         msq13_2(:,:)=msq_mix13_2(1,:,:)+msq_mix13_2(2,:,:)
         msq14_2(:,:)=msq_mix14_2(1,:,:)+msq_mix14_2(2,:,:)
         msqc(5,j,k)=sub15_2(qg)*(
     &      msq15_2(+1,k)+msq15_2(+2,k)
     &     +msq15_2(+3,k)+msq15_2(+4,k))
         msqc(9,j,k)=sub13_2(qg)*(
     &      msq13_2(-1,k)+msq13_2(-2,k)
     &     +msq13_2(-3,k)+msq13_2(-4,k)
     &     -msq13_2(-k,k)+0.5_dp*msq13_2(-k,k))
         msqc(11,j,k)=sub14_2(qg)*msq14_2(-k,k)
      endif

      enddo
      enddo

      return
      end
      


      subroutine store_msq_mix(msqto)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'msq_mix.f'
      integer i,j,k
      real(dp):: msqto(0:3,-nf:nf,-nf:nf)
      
      do i=0,3
      do j=-nf,nf
      do k=-nf,nf
      msqto(i,j,k)=msq_mix(i,j,k)
      enddo
      enddo
      enddo
      
      return
      end
      
