      subroutine qqb_QQb_mix_gs(p,msqc)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c--- for mixed QCD/Z matrix element for top pair production,
c     q(-p1)+qbar(-p2) -->  Q(p3)+Q~(P4) + g(p5)
c--- all momenta are incoming
c
      implicit none 
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'    
      integer j,k
c --- remember: nd will count the dipoles
      integer nd
c--- slightly obtuse notation, to simplify declaration lines
      real(dp):: p(mxpart,4),msqc(maxd,fn:nf,fn:nf)
      real(dp)::
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq15_3(fn:nf,fn:nf),msq15_4(fn:nf,fn:nf),
     & msq35_1(fn:nf,fn:nf),msq45_1(fn:nf,fn:nf),
     & msq25_3(fn:nf,fn:nf),msq25_4(fn:nf,fn:nf),
     & msq35_2(fn:nf,fn:nf),msq45_2(fn:nf,fn:nf),

     & msq15_2v(fn:nf,fn:nf),msq25_1v(fn:nf,fn:nf),
     & msq15_3v(fn:nf,fn:nf),msq15_4v(fn:nf,fn:nf),
     & msq35_1v(fn:nf,fn:nf),msq45_1v(fn:nf,fn:nf),
     & msq25_3v(fn:nf,fn:nf),msq25_4v(fn:nf,fn:nf),
     & msq35_2v(fn:nf,fn:nf),msq45_2v(fn:nf,fn:nf),

     & sub15_2(4),sub25_1(4),
     & sub15_3(4),sub15_4(4),
     & sub35_1(4),sub45_1(4),
     & sub25_3(4),sub25_4(4),
     & sub35_2(4),sub45_2(4),

     & sub15_2v,sub25_1v,
     & sub15_3v,sub15_4v,
     & sub35_1v,sub45_1v,
     & sub25_3v,sub25_4v,
     & sub35_2v,sub45_2v


      real(dp)::
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m15_3(0:2,fn:nf,fn:nf),m35_1(0:2,fn:nf,fn:nf),
     & m15_4(0:2,fn:nf,fn:nf),m45_1(0:2,fn:nf,fn:nf),
     & m25_3(0:2,fn:nf,fn:nf),m35_2(0:2,fn:nf,fn:nf),
     & m25_4(0:2,fn:nf,fn:nf),m45_2(0:2,fn:nf,fn:nf)

      real(dp)::
     & m15_2v(0:2,fn:nf,fn:nf),m25_1v(0:2,fn:nf,fn:nf),
     & m15_3v(0:2,fn:nf,fn:nf),m35_1v(0:2,fn:nf,fn:nf),
     & m15_4v(0:2,fn:nf,fn:nf),m45_1v(0:2,fn:nf,fn:nf),
     & m25_3v(0:2,fn:nf,fn:nf),m35_2v(0:2,fn:nf,fn:nf),
     & m25_4v(0:2,fn:nf,fn:nf),m45_2v(0:2,fn:nf,fn:nf)

      external qqb_QQb_mix,donothing_gvec

      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

      ndmax=8

c-- initialize the matrix elements to zero
      msqc(:,:,:)=zip
      
c--initial final
      call dips_mass(1,p,1,5,3,sub15_3,sub15_3v,msq15_3,msq15_3v,
     . qqb_QQb_mix,donothing_gvec)
      call dips_mass(2,p,2,5,3,sub25_3,sub25_3v,msq25_3,msq25_3v,
     . qqb_QQb_mix,donothing_gvec)
      call dips_mass(3,p,1,5,4,sub15_4,sub15_4v,msq15_4,msq15_4v,
     . qqb_QQb_mix,donothing_gvec)
      call dips_mass(4,p,2,5,4,sub25_4,sub25_4v,msq25_4,msq25_4v,
     . qqb_QQb_mix,donothing_gvec)

c--final-initial 
      call dips_mass(5,p,3,5,1,sub35_1,sub35_1v,msq35_1,msq35_1v,
     . qqb_QQb_mix,donothing_gvec)
      call dips_mass(6,p,3,5,2,sub35_2,sub35_2v,msq35_2,msq35_2v,
     . qqb_QQb_mix,donothing_gvec)
      call dips_mass(7,p,4,5,1,sub45_1,sub45_1v,msq45_1,msq45_1v,
     . qqb_QQb_mix,donothing_gvec)
      call dips_mass(8,p,4,5,2,sub45_2,sub45_2v,msq45_2,msq45_2v,
     . qqb_QQb_mix,donothing_gvec)

c      write(6,*) 'gs:sub15_2(gg)',sub15_2(gg)
c      write(6,*) 'gs:sub25_1(gg)',sub25_1(gg)

c      write(6,*) 'gs:sub15_3(gg)',sub15_3(gg)
c      write(6,*) 'gs:sub15_4(gg)',sub15_4(gg)
c      write(6,*) 'gs:sub25_3(gg)',sub25_3(gg)
c      write(6,*) 'gs:sub25_4(gg)',sub25_4(gg)

c      write(6,*) 'gs:sub35_1(qq)',sub35_1(qq)
c      write(6,*) 'gs:sub35_2(qq)',sub35_2(qq)
c      write(6,*) 'gs:sub45_1(qq)',sub45_1(qq)
c      write(6,*) 'gs:sub45_2(qq)',sub45_2(qq)

c      write(6,*) 'gs:sub35_4(qq)',sub35_4(qq)
c      write(6,*) 'gs:sub45_3(qq)',sub45_3(qq)


c--- fill the dipole contributions
      do j=-nf,nf
      k=-j

c--- note relative signs, so no collinear contribution only soft
      if (j > 0) then
        msqc(1,j,k)= +four*msq15_3(j,k)*sub15_3(qq)
        msqc(2,j,k)= -four*msq25_3(j,k)*sub25_3(qq)
        msqc(3,j,k)= -four*msq15_4(j,k)*sub15_4(qq)
        msqc(4,j,k)= +four*msq25_4(j,k)*sub25_4(qq)

        msqc(5,j,k)= +four*msq35_1(j,k)*sub35_1(qq)
        msqc(6,j,k)= -four*msq35_2(j,k)*sub35_2(qq)
        msqc(7,j,k)= -four*msq45_1(j,k)*sub45_1(qq)
        msqc(8,j,k)= +four*msq45_2(j,k)*sub45_2(qq)
      elseif (j < 0) then
        msqc(1,j,k)= -four*msq15_3(j,k)*sub15_3(qq)
        msqc(2,j,k)= +four*msq25_3(j,k)*sub25_3(qq)
        msqc(3,j,k)= +four*msq15_4(j,k)*sub15_4(qq)
        msqc(4,j,k)= -four*msq25_4(j,k)*sub25_4(qq)

        msqc(5,j,k)= -four*msq35_1(j,k)*sub35_1(qq)
        msqc(6,j,k)= +four*msq35_2(j,k)*sub35_2(qq)
        msqc(7,j,k)= +four*msq45_1(j,k)*sub45_1(qq)
        msqc(8,j,k)= -four*msq45_2(j,k)*sub45_2(qq)
      endif

      enddo

      return
      end
      
