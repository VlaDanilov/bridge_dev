      subroutine pvDfill(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N)
      implicit none
C     N is the offset in the common block
C     p1,p2,p3,p4 are the invariant masses sqaured of external lines
C     m1,m2,m3,m4 are the masses squared of internal lines
      include 'types.f'
      include 'TRconstants.f'
      include 'TRscale.f'
      include 'pvCnames.f'
      include 'pvDnames.f'
      include 'pvCv.f'
      include 'pvDv.f'
      include 'pvverbose.f'
      include 'TRmaxindex.f'
      include 'pvRespectmaxcindex.f'
      include 'pvrecurflags.f'
      include 'pvforcerecalc.f'
      integer:: C234,C134,C124,C123,ep,epmj,N,j,perm(3),pvCcache
      integer,parameter::np=3
      real(dp):: p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,f1,f2,f3
      complex(dp) G(np,np),csum(-2:0),c1sum(-2:0),c2sum(-2:0),
     . c0sum(-2:0),in(3,-2:0),trI4
      complex(dp):: c11sum(-2:0),c00sum(-2:0),c12sum(-2:0),c22sum(-2:0)
      logical:: exceptional
      logical,save:: first=.true.
      real(dp):: q1save(4),q2save(4),q3save(4)
      real(dp):: triq1save(4),triq2save(4)
      common/q123save/q1save,q2save,q3save  
      common/q12save/triq1save,triq2save  
!$omp threadprivate(first,/q123save/,/q12save/)
      real(dp),save:: idp2(0:2),idp1(0:2),id(0:2),idm1(0:2),
     & idm2(0:2),idm3(0:2)
!$omp threadprivate(idp2,idp1,id,idm1,idm2,idm3)
      integer,save:: icall,irecur,irecur2,irecur3,irecur4
!$omp threadprivate(icall,irecur,irecur2,irecur3,irecur4)
      include 'cplx.h'
     
      if (first) then
      first=.false.
C--idp2=1/[D+2]
      idp2(0)=1._dp/six
      idp2(1)=idp2(0)/3._dp
      idp2(2)=idp2(1)/3._dp
C--idp1=1/[D+1]
      idp1(0)=0.2_dp
      idp1(1)=idp1(0)*0.4_dp
      idp1(2)=idp1(1)*0.4_dp
C--id=1/D
      id(0)=0.25_dp
      id(1)=id(0)*half
      id(2)=id(1)*half
C--idm1=1/[D-1]
      idm1(0)=1._dp/3._dp
      idm1(1)=idm1(0)*2._dp/3._dp
      idm1(2)=idm1(1)*2._dp/3._dp
C--idm2=1/[D-2]
      idm2(0)=half
      idm2(1)=idm2(0)
      idm2(2)=idm2(1)
C--idm3=1/[D-3]
      idm3(0)=1._dp
      idm3(1)=2._dp*idm3(0)
      idm3(2)=4._dp*idm3(1)
c--- variables for statistics reporting
      irecur=0
      irecur2=0
      irecur3=0
      irecur4=0
      icall=0      
c--- print out flags for recursion
c      write(6,*) 'pvDfill recursion flags:'
c      write(6,*) '  doGsing  ',doGsing
c      write(6,*) '  doGYsing ',doGYsing
c      write(6,*) '  doPsing  ',doPsing
c      write(6,*) '  doPFsing ',doPFsing
      endif

c--- statistics accounting and reporting      
      icall=icall+1
      if (pvverbose) then
      if (mod(icall,50000) .eq. 0) then
        write(6,77) icall,
     & 1d2*dfloat(icall-irecur-irecur2-irecur3-irecur4)/dfloat(icall),
     & 1d2*dfloat(irecur)/dfloat(icall),
     & 1d2*dfloat(irecur2)/dfloat(icall),
     & 1d2*dfloat(irecur3)/dfloat(icall),
     & 1d2*dfloat(irecur4)/dfloat(icall)
      endif
      endif
   77 format(' +++ Dfill ',i9,': ',5(f6.2,'% : '))
      
      f1 = m2 - m1 - p1
      f2 = m3 - m1 - p1p2
      f3 = m4 - m1 - p4

!---double up the Gram matrix to remove factors of 1/2 in Eqs.
      G(1,1) = cplx1(two*p1)
      G(2,2) = cplx1(two*p1p2)
      G(3,3) = cplx1(two*p4)
      G(1,2) = cplx1(p1+p1p2 - p2)
      G(2,1) = G(1,2)
      G(1,3) = cplx1(p1+p4 - p2p3)
      G(3,1) = G(1,3)
      G(2,3) = cplx1(p1p2 - p3+p4)
      G(3,2) = G(2,3)

c--- Check for small kinematic quantities requiring alternate recursion
c      if (pvverbose) write(6,*) 'Check box Gsing'
c      Gsing=pvGramsing(G,3)

C     Y(i,j)=mi^2+mj^2-(q_i-q_j)^2
C     where q_1=0,  q_2=p1,  q_3=p_1+p_2, q_4=p_1+p_2+p_3;

c      Y(1,1) = cplx1(two*m1)
c      Y(1,2) = cplx1(m1 + m2 - p1)
c      Y(2,1) = Y(1,2)
c      Y(1,3) = cplx1(m1 + m3 - p1p2)
c      Y(3,1) = Y(1,3)
c      Y(1,4) = cplx1(m1 + m4 - p4)
c      Y(4,1) = Y(1,4)
c      Y(2,2) = cplx1(two*m2)
c      Y(2,3) = cplx1(m2 + m3 - p2)
c      Y(3,2) = Y(2,3)
c      Y(2,4) = cplx1(m2 + m4 - p2p3)
c      Y(4,2) = Y(2,4)
c      Y(3,3) = cplx1(two*m3)
c      Y(3,4) = cplx1(m3 + m4 - p3)
c      Y(4,3) = Y(3,4)
c      Y(4,4) = cplx1(two*m4)
      
c      if (pvverbose) write(6,*) 'Check box Ysing'
c      Ysing=pvGramsing(Y,4)

c-- find maximum entry in Gram matrix
c      Gmax=zip
c      do j=1,3
c      do k=j,3
c      if (abs(G(j,k)) .gt. Gmax) Gmax=abs(G(j,k))
c      enddo
c      enddo
      
c      if (pvverbose) write(6,*) 'Gmax=',Gmax
c      Psing=.false.
c--- criterion for small momenta recursion
c      if (Gmax .lt. weenumber) Psing=.true.

c-- find maximum of f1, f2 and f3
c      fmax=max(abs(f1),abs(f2),abs(f3))
      
c      if (pvverbose) write(6,*) 'fmax=',fmax
c      Fsing=.false.
c--- criterion for small momenta and small f(k) recursion
c      if (fmax .lt. weenumber) Fsing=.true.

c      if (pvverbose) write(6,*) 'pvDfill: Gsing,Ysing,Psing,Fsing',
c     & Gsing,Ysing,Psing,Fsing

c--- If alternative recursion required, ensure triangle integrals
c--- are computed to maximum available rank
      if (doPFsing) pvRespectmaxcindex=.false.
      if (doPsing)  pvRespectmaxcindex=.false.
      if (doGYsing) pvRespectmaxcindex=.false.
      if (doGsing)  pvRespectmaxcindex=.false.
            
      pvRespectmaxcindex=.true.

      pvforcerecalc=.false.
c--- Set up relevant triangle pinchings
      triq1save(:)=q2save(:)-q1save(:)
      triq2save(:)=q3save(:)-q1save(:)
      C234=pvCcache(p2,p3,p2p3,m2,m3,m4)
      triq1save(:)=q2save(:)
      triq2save(:)=q3save(:)
      C134=pvCcache(p1p2,p3,p4,m1,m3,m4)
      triq1save(:)=q1save(:)
      triq2save(:)=q3save(:)
      C124=pvCcache(p1,p2p3,p4,m1,m2,m4)
      triq1save(:)=q1save(:)
      triq2save(:)=q2save(:)
      C123=pvCcache(p1,p2,p1p2,m1,m2,m3)

c--- Return to default behaviour
      pvRespectmaxcindex=.true.

c--- Make call to alternate recursion routines if required       
      exceptional=.false.
      
      if     (doPFsing) then
c--- for small momenta and small f(k)
        if (pvverbose) then
          write(6,*) 'USING BOX SMALL MOMENTA AND f(k) RECURSION'
      endif
        call Dfill_recur4(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N)
      irecur4=irecur4+1
      return
      elseif (doPsing) then
c--- for small momenta
        if (pvverbose) then
          write(6,*) 'USING BOX SMALL MOMENTA RECURSION'
      endif
        call Dfill_recur3(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N)
      irecur3=irecur3+1
      return
      elseif (doGYsing) then
c--- for small Gram and small Y
        if (pvverbose) then
          write(6,*) 'USING BOX SMALL Y AND SMALL G RECURSION'
      endif
        call Dfill_recur2(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N,
     .                     exceptional)
      irecur2=irecur2+1
        if (exceptional) then
c------ for exceptional configurations, fall through to normal PV
        continue
      else
c------ otherwise, we're done
        return
      endif
      elseif (doGsing) then
c--- for small Gram only  
        if (pvverbose) then
          write(6,*) 'USING BOX SMALL G RECURSION'
      endif
        call Dfill_recur (p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N)
      irecur=irecur+1
      return
      endif
c--- otherwise, usual PV is fine

c      if (exceptional) write(6,*) 'WARNING: EXCEPTIONAL POINT'

c--- initialize integrals      
      do ep=-2,0
      do j=1,Ndd
      Dv(N+j,ep)=cplx2(1d5,-1d5)
      enddo
      enddo

      call XLUDecomp(G, 3, perm)
      
      do j=1,Ndd
      Dv(N+j,ep)=cplx1(10000._dp)
      enddo

      do ep=-2,0
      Dv(N+dd0,ep) =trI4(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,musq,ep)
      enddo

      c0sum(:)=Cv(cc0+C234,:)+Cv(cc1+C234,:)+Cv(cc2+C234,:)
      c1sum(:)=Cv(cc1+C234,:)+Cv(cc11+C234,:)+Cv(cc12+C234,:)
      c2sum(:)=Cv(cc2+C234,:)+Cv(cc12+C234,:)+Cv(cc22+C234,:)
      csum(:)=c0sum(:)+c1sum(:)+c2sum(:)

      c00sum(:) = Cv(cc00+C234,:) +
     &    Cv(cc001+C234,:)+Cv(cc002+C234,:)
      c11sum(:) = Cv(cc11+C234,:) +
     &    Cv(cc111+C234,:)+Cv(cc112+C234,:)
      c12sum(:) = Cv(cc12+C234,:) +
     &    Cv(cc112+C234,:)+Cv(cc122+C234,:)
      c22sum(:) = Cv(cc22+C234,:) +
     &    Cv(cc122+C234,:)+Cv(cc222+C234,:)


      in(1,:) = f1*Dv(N+dd0,:) - Cv(cc0+C234,:)+Cv(cc0+C134,:)
      in(2,:) = f2*Dv(N+dd0,:) - Cv(cc0+C234,:)+Cv(cc0+C124,:)
      in(3,:) = f3*Dv(N+dd0,:) - Cv(cc0+C234,:)+Cv(cc0+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1,:) = in(1,:)
      Dv(N+dd2,:) = in(2,:)
      Dv(N+dd3,:) = in(3,:)

      do ep=-2,0
      Dv(N+dd00,ep) = czip
      if (ep .eq. -2) goto 20
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd00,ep) = Dv(N+dd00,ep)+idm3(j)*(m1*Dv(N+dd0,epmj)
     &  -half*(Dv(N+dd1,epmj)*f1 +Dv(N+dd2,epmj)*f2 +Dv(N+dd3,epmj)*f3 
     & - Cv(cc0+C234,epmj)))
      enddo
 20   continue
      enddo 

      in(1,:) = f1*Dv(N+dd1,:)+c0sum(:) - two*Dv(N+dd00,:)
      in(2,:) = f2*Dv(N+dd1,:)+c0sum(:)+Cv(cc1+C124,:)
      in(3,:) = f3*Dv(N+dd1,:)+c0sum(:)+Cv(cc1+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd11,:) = in(1,:)
      Dv(N+dd12,:) = in(2,:)
      Dv(N+dd13,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd2,:) - Cv(cc1+C234,:)+Cv(cc1+C134,:)
      in(2,:) = f2*Dv(N+dd2,:) - Cv(cc1+C234,:) - two*Dv(N+dd00,:)
      in(3,:) = f3*Dv(N+dd2,:) - Cv(cc1+C234,:)+Cv(cc2+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd12,:) = half*(Dv(N+dd12,:)+in(1,:))
      Dv(N+dd22,:) = in(2,:)
      Dv(N+dd23,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd3,:) - Cv(cc2+C234,:)+Cv(cc2+C134,:)
      in(2,:) = f2*Dv(N+dd3,:) - Cv(cc2+C234,:)+Cv(cc2+C124,:)
      in(3,:) = f3*Dv(N+dd3,:) - Cv(cc2+C234,:) - two*Dv(N+dd00,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd13,:) = half*(Dv(N+dd13,:)+in(1,:))
      Dv(N+dd23,:) = half*(Dv(N+dd23,:)+in(2,:))
      Dv(N+dd33,:) = in(3,:)

C --- end of two index tensors

      if (maxdindex .eq. 2) return

C--- three index tensors
      do ep=-2,0
      do j=dd001,dd003
      Dv(N+j,ep) = czip
      enddo
      if (ep .eq. -2) goto 30
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd001,ep)=Dv(N+dd001,ep)+idm2(j)*half*(-f1*Dv(N +
     &    dd11,epmj)-f2*Dv(N+dd12,epmj)-f3*Dv(N+dd13,epmj
     &    )-Cv(cc0+C234,epmj)-Cv(cc1+C234,epmj)-Cv(cc2 +
     &    C234,epmj)+two*m1*Dv(N+dd1,epmj))
      Dv(N+dd002,ep)=Dv(N+dd002,ep)+idm2(j)*half*(-f1*Dv(N +
     &    dd12,epmj)-f2*Dv(N+dd22,epmj)-f3*Dv(N+dd23,epmj
     &    )+Cv(cc1+C234,epmj)+two*m1*Dv(N+dd2,epmj))
      Dv(N+dd003,ep)=Dv(N+dd003,ep)+idm2(j)*half*(-f1*Dv(N +
     &    dd13,epmj)-f2*Dv(N+dd23,epmj)-f3*Dv(N+dd33,epmj
     &    )+Cv(cc2+C234,epmj)+two*m1*Dv(N+dd3,epmj))
      enddo
 30   continue
      enddo

      in(1,:) = f1*Dv(N+dd11,:) - csum(:) - four*Dv(N+dd001,:)
      in(2,:) = f2*Dv(N+dd11,:) - csum(:)+Cv(cc11+C124,:)
      in(3,:) = f3*Dv(N+dd11,:) - csum(:)+Cv(cc11+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd111,:) = in(1,:)
      Dv(N+dd112,:) = in(2,:)
      Dv(N+dd113,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd22,:) - Cv(cc11+C234,:)+Cv(cc11+C134,:)
      in(2,:) = f2*Dv(N+dd22,:) - Cv(cc11+C234,:)
     . - four*Dv(N+dd002,:)
      in(3,:) = f3*Dv(N+dd22,:) - Cv(cc11+C234,:)+Cv(cc22+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd122,:) = in(1,:)
      Dv(N+dd222,:) = in(2,:)
      Dv(N+dd223,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd33,:) - Cv(cc22+C234,:)+Cv(cc22+C134,:)
      in(2,:) = f2*Dv(N+dd33,:) - Cv(cc22+C234,:)+Cv(cc22+C124,:)
      in(3,:) = f3*Dv(N+dd33,:) - Cv(cc22+C234,:)
     . - four*Dv(N+dd003,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd133,:) = in(1,:)
      Dv(N+dd233,:) = in(2,:)
      Dv(N+dd333,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd13,:)+c2sum(:) - two*Dv(N+dd003,:)
      in(2,:) = f2*Dv(N+dd13,:)+c2sum(:)+Cv(cc12+C124,:)
      in(3,:) = f3*Dv(N+dd13,:)+c2sum(:) - two*Dv(N+dd001,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd113,:) = half*(Dv(N+dd113,:)+in(1,:))
      Dv(N+dd123,:) = in(2,:)
      Dv(N+dd133,:) = half*(Dv(N+dd133,:)+in(3,:))

c--- check the contents of box array    
c      write(6,*) 'PV: D array'
c      do j=1,24
c        write(6,'(i3,2e20.12)') j,Dv(j+N,0)
c      enddo

      if (maxdindex .eq. 3) return

C--- four index tensors
      do ep=-2,0
      do j=dd0000,dd0033
      Dv(N+j,ep) = czip
      enddo
      if (ep .eq. -2) goto 40
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd0000,ep) = Dv(N+dd0000,ep)+idm1(j)*(m1*Dv(N+dd00,epmj) -
     &  half
     & *(f1*Dv(N+dd001,epmj)+f2*Dv(N+dd002,epmj)+f3*Dv(N+dd003,epmj)-
     &    Cv(cc00+C234,epmj)))
      Dv(N+dd0011,ep) = Dv(N+dd0011,ep)+idm1(j)*(m1*Dv(N+dd11,epmj)-
     &  half
     & *(f1*Dv(N+dd111,epmj)+f2*Dv(N+dd112,epmj)+f3*Dv(N+dd113,epmj)
     & -csum(epmj)))
      Dv(N+dd0012,ep) = Dv(N+dd0012,ep)+idm1(j)*(m1*Dv(N+dd12,epmj)-
     &  half
     & *(f1*Dv(N+dd112,epmj)+f2*Dv(N+dd122,epmj)+f3*Dv(N+dd123,epmj)
     & +c1sum(epmj)))
      Dv(N+dd0013,ep) = Dv(N+dd0013,ep)+idm1(j)*(m1*Dv(N+dd13,epmj)-
     &  half
     & *(f1*Dv(N+dd113,epmj)+f2*Dv(N+dd123,epmj)+f3*Dv(N+dd133,epmj)
     & +c2sum(epmj)))
      Dv(N+dd0022,ep) = Dv(N+dd0022,ep)+idm1(j)*(m1*Dv(N+dd22,epmj)-
     &  half
     & *(f1*Dv(N+dd122,epmj)+f2*Dv(N+dd222,epmj)+f3*Dv(N+dd223,epmj)-
     &    Cv(cc11+C234,epmj)))
      Dv(N+dd0023,ep) = Dv(N+dd0023,ep)+idm1(j)*(m1*Dv(N+dd23,epmj)-
     &  half
     & *(f1*Dv(N+dd123,epmj)+f2*Dv(N+dd223,epmj)+f3*Dv(N+dd233,epmj)-
     &    Cv(cc12+C234,epmj)))
      Dv(N+dd0033,ep) = Dv(N+dd0033,ep)+idm1(j)*(m1*Dv(N+dd33,epmj)-
     &  half
     & *(f1*Dv(N+dd133,epmj)+f2*Dv(N+dd233,epmj)+f3*Dv(N+dd333,epmj)-
     &    Cv(cc22+C234,epmj)))
      enddo
 40   continue
      enddo

      c1sum(:) = c1sum(:)+c11sum(:)+c12sum(:)
      c2sum(:) = c2sum(:)+c12sum(:)+c22sum(:)
      csum(:) = csum(:)+c1sum(:)+c2sum(:)

      in(1,:) = f1*Dv(N+dd111,:)+csum(:) - 6d0*Dv(N+dd0011,:)
      in(2,:) = f2*Dv(N+dd111,:)+csum(:)+Cv(cc111+C124,:)
      in(3,:) = f3*Dv(N+dd111,:)+csum(:)+Cv(cc111+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1111,:) = in(1,:)
      Dv(N+dd1112,:) = in(2,:)
      Dv(N+dd1113,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd113,:) - c2sum(:) - four*Dv(N+dd0013,:)
      in(2,:) = f2*Dv(N+dd113,:) - c2sum(:)+Cv(cc112+C124,:)
      in(3,:) = f3*Dv(N+dd113,:) - c2sum(:) - two*Dv(N+dd0011,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1113,:) = half*(Dv(N+dd1113,:)+in(1,:))
      Dv(N+dd1123,:) = in(2,:)
      Dv(N+dd1133,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd122,:)+c11sum(:) - two*Dv(N+dd0022,:)
      in(2,:) = f2*Dv(N+dd122,:)+c11sum(:) - four*Dv(N+dd0012,:)
      in(3,:) = f3*Dv(N+dd122,:)+c11sum(:)+Cv(cc122+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1122,:) = in(1,:)
      Dv(N+dd1222,:) = in(2,:)
      Dv(N+dd1223,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd222,:)-Cv(cc111+C234,:)+Cv(cc111+C134,:)
      in(2,:) = f2*Dv(N+dd222,:)-Cv(cc111+C234,:)-6*Dv(N+dd0022,:)
      in(3,:) = f3*Dv(N+dd222,:)-Cv(cc111+C234,:)+Cv(cc222+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1222,:) = half*(Dv(N+dd1222,:)+in(1,:))
      Dv(N+dd2222,:) = in(2,:)
      Dv(N+dd2223,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd233,:)-Cv(cc122+C234,:)+Cv(cc122+C134,:)
      in(2,:) = f2*Dv(N+dd233,:)-Cv(cc122+C234,:)-2*Dv(N+dd0033,:)
      in(3,:) = f3*Dv(N+dd233,:)-Cv(cc122+C234,:)-4*Dv(N+dd0023,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1233,:) = in(1,:)
      Dv(N+dd2233,:) = in(2,:)
      Dv(N+dd2333,:) = in(3,:)

      in(1,:) = f1*Dv(N+dd333,:)-Cv(cc222+C234,:)+Cv(cc222+C134,:)
      in(2,:) = f2*Dv(N+dd333,:)-Cv(cc222+C234,:)+Cv(cc222+C124,:)
      in(3,:) = f3*Dv(N+dd333,:)-Cv(cc222+C234,:)-6*Dv(N+dd0033,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1333,:) = in(1,:)
      Dv(N+dd2333,:) = half*(Dv(N+dd2333,:)+in(2,:))
      Dv(N+dd3333,:) = in(3,:)
      c00sum(:) = c00sum(:) +
     &    Cv(cc001+C234,:)+Cv(cc0011+C234,:)+Cv(cc0012+C234,:) +
     &    Cv(cc002+C234,:)+Cv(cc0012+C234,:)+Cv(cc0022+C234,:)
      c11sum(:) = c11sum(:) +
     &    Cv(cc111+C234,:)+Cv(cc1111+C234,:)+Cv(cc1112+C234,:) +
     &    Cv(cc112+C234,:)+Cv(cc1112+C234,:)+Cv(cc1122+C234,:)
      c12sum(:) = c12sum(:) +
     &    Cv(cc112+C234,:)+Cv(cc1112+C234,:)+Cv(cc1122+C234,:) +
     &    Cv(cc122+C234,:)+Cv(cc1122+C234,:)+Cv(cc1222+C234,:)
      c22sum(:) = c22sum(:) +
     &    Cv(cc122+C234,:)+Cv(cc1122+C234,:)+Cv(cc1222+C234,:) +
     &    Cv(cc222+C234,:)+Cv(cc1222+C234,:)+Cv(cc2222+C234,:)
      c1sum(:) = c1sum(:)+c11sum(:)+c12sum(:)
      c2sum(:) = c2sum(:)+c12sum(:)+c22sum(:)
      csum(:) = csum(:)+c1sum(:)+c2sum(:)


      if (maxdindex .eq. 4) return

C------begin of five index tensors
      do ep=-2,0
      do j=dd00001,dd00333
      Dv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 50
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd00001,ep)=Dv(N+dd00001,ep)+id(j)*half*(-f1*Dv(N+
     &    dd0011,epmj)-f2*Dv(N+dd0012,epmj)-f3*Dv(N+dd0013,
     &    epmj)-Cv(cc00+C234,epmj)-Cv(cc001+C234,epmj)-
     &    Cv(cc002+C234,epmj)+two*m1*Dv(N+dd001,epmj))

      Dv(N+dd00002,ep)=Dv(N+dd00002,ep)+id(j)*half*(-f1*Dv(N+
     &    dd0012,epmj)-f2*Dv(N+dd0022,epmj)-f3*Dv(N+dd0023,
     &    epmj)+Cv(cc001+C234,epmj)+two*m1*Dv(N+dd002,epmj))

      Dv(N+dd00003,ep)=Dv(N+dd00003,ep)+id(j)*half*(-f1*Dv(N+
     &    dd0013,epmj)-f2*Dv(N+dd0023,epmj)-f3*Dv(N+dd0033,
     &    epmj)+Cv(cc002+C234,epmj)+two*m1*Dv(N+dd003,epmj))

      Dv(N+dd00111,ep)=Dv(N+dd00111,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1111,epmj)-f2*Dv(N+dd1112,epmj)-f3*Dv(N+dd1113,
     &    epmj)-Cv(cc0+C234,epmj)-three*Cv(cc1+C234,epmj)
     &    -three*Cv(cc11+C234,epmj)-Cv(cc111+C234,epmj)-three
     &    *Cv(cc112+C234,epmj)-six*Cv(cc12+C234,epmj)-three
     &    *Cv(cc122+C234,epmj)-three*Cv(cc2+C234,epmj)-three*
     &    Cv(cc22+C234,epmj)-Cv(cc222+C234,epmj)
     &    +two*m1*Dv(N+dd111,epmj))

      Dv(N+dd00112,ep)=Dv(N+dd00112,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1112,epmj)-f2*Dv(N+dd1122,epmj)-f3*Dv(N+dd1123,
     &    epmj)+Cv(cc1+C234,epmj)+two*Cv(cc11+C234,epmj)
     &    +Cv(cc111+C234,epmj)+two*Cv(cc112+C234,epmj)+two
     &    *Cv(cc12+C234,epmj)+Cv(cc122+C234,epmj)
     &    +two*m1*Dv(N+dd112,epmj))

      Dv(N+dd00113,ep)=Dv(N+dd00113,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1113,epmj)-f2*Dv(N+dd1123,epmj)-f3*Dv(N+dd1133,
     &    epmj)+Cv(cc112+C234,epmj)+two*Cv(cc12+C234,ep-
     &    j)+two*Cv(cc122+C234,epmj)+Cv(cc2+C234,epmj)+two
     &    *Cv(cc22+C234,epmj)+Cv(cc222+C234,epmj)
     &    +two*m1*Dv(N+dd113,epmj))

      Dv(N+dd00122,ep)=Dv(N+dd00122,ep)+id(j)*half *(-f1*Dv(N+
     &    dd1122,epmj)-f2*Dv(N+dd1222,epmj)-f3*Dv(N+dd1223,
     &    epmj)-Cv(cc11+C234,epmj)-Cv(cc111+C234,epmj)-
     &    Cv(cc112+C234,epmj)
     &    +two*m1*Dv(N+dd122,epmj))

      Dv(N+dd00123,ep)=Dv(N+dd00123,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1123,epmj)-f2*Dv(N+dd1223,epmj)-f3*Dv(N+dd1233,
     &    epmj)-Cv(cc112+C234,epmj)-Cv(cc12+C234,epmj)-
     &    Cv(cc122+C234,epmj)
     &    +two*m1*Dv(N+dd123,epmj))

      Dv(N+dd00133,ep)=Dv(N+dd00133,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1133,epmj)-f2*Dv(N+dd1233,epmj)-f3*Dv(N+dd1333,
     &    epmj)-Cv(cc122+C234,epmj)-Cv(cc22+C234,epmj)-
     &    Cv(cc222+C234,epmj)
     &    +two*m1*Dv(N+dd133,epmj))

      Dv(N+dd00222,ep)=Dv(N+dd00222,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1222,epmj)-f2*Dv(N+dd2222,epmj)-f3*Dv(N+dd2223,
     &    epmj)+Cv(cc111+C234,epmj)
     &    +two*m1*Dv(N+dd222,epmj))

      Dv(N+dd00223,ep)=Dv(N+dd00223,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1223,epmj)-f2*Dv(N+dd2223,epmj)-f3*Dv(N+dd2233,
     &    epmj)+Cv(cc112+C234,epmj)
     &    +two*m1*Dv(N+dd223,epmj))

      Dv(N+dd00233,ep)=Dv(N+dd00233,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1233,epmj)-f2*Dv(N+dd2233,epmj)-f3*Dv(N+dd2333,
     &    epmj)+Cv(cc122+C234,epmj)
     &    +two*m1*Dv(N+dd233,epmj))

      Dv(N+dd00333,ep)=Dv(N+dd00333,ep)+id(j)*half*(-f1*Dv(N+
     &    dd1333,epmj)-f2*Dv(N+dd2333,epmj)-f3*Dv(N+dd3333,
     &    epmj)+Cv(cc222+C234,epmj)
     &    +two*m1*Dv(N+dd333,epmj))

      enddo
 50   continue
      enddo

      in(1,:)=f1*Dv(N+dd1111,:)-csum(:)-8*Dv(N+dd00111,:)
      in(2,:)=f2*Dv(N+dd1111,:)-csum(:)+Cv(cc1111+C124,:)
      in(3,:)=f3*Dv(N+dd1111,:)-csum(:)+Cv(cc1111+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd11111,:)=in(1,:)
      Dv(N+dd11112,:)=in(2,:)
      Dv(N+dd11113,:)=in(3,:)

      in(1,:)=f1*Dv(N+dd2222,:)-Cv(cc1111+C234,:)+Cv(cc1111+C134,:)
      in(2,:)=f2*Dv(N+dd2222,:)-Cv(cc1111+C234,:)-8*Dv(N+dd00222,:)
      in(3,:)=f3*Dv(N+dd2222,:)-Cv(cc1111+C234,:)+Cv(cc2222+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd12222,:)=in(1,:)
      Dv(N+dd22222,:)=in(2,:)
      Dv(N+dd22223,:)=in(3,:)

      in(1,:)=f1*Dv(N+dd3333,:)-Cv(cc2222+C234,:)+Cv(cc2222+C134,:)
      in(2,:)=f2*Dv(N+dd3333,:)-Cv(cc2222+C234,:)+Cv(cc2222+C124,:)
      in(3,:)=f3*Dv(N+dd3333,:)-Cv(cc2222+C234,:)-8*Dv(N+dd00333,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd13333,:)=in(1,:)
      Dv(N+dd23333,:)=in(2,:)
      Dv(N+dd33333,:)=in(3,:)

      in(1,:)=f1*Dv(N+dd1122,:)-c11sum(:)-4*Dv(N+dd00122,:)
      in(2,:)=f2*Dv(N+dd1122,:)-c11sum(:)-4*Dv(N+dd00112,:)
      in(3,:)=f3*Dv(N+dd1122,:)-c11sum(:)+Cv(cc1122+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd11122,:)=in(1,:)
      Dv(N+dd11222,:)=in(2,:)
      Dv(N+dd11223,:)=in(3,:)

      in(1,:)=f1*Dv(N+dd1133,:)-c22sum(:)-4*Dv(N+dd00133,:)
      in(2,:)=f2*Dv(N+dd1133,:)-c22sum(:)+Cv(cc1122+C124,:)
      in(3,:)=f3*Dv(N+dd1133,:)-c22sum(:)-4*Dv(N+dd00113,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd11133,:)=in(1,:)
      Dv(N+dd11233,:)=in(2,:)
      Dv(N+dd11333,:)=in(3,:)

      in(1,:)=f1*Dv(N+dd2233,:)-Cv(cc1122+C234,:)+Cv(cc1122+C134,:)
      in(2,:)=f2*Dv(N+dd2233,:)-Cv(cc1122+C234,:)-4*Dv(N+dd00233,:)
      in(3,:)=f3*Dv(N+dd2233,:)-Cv(cc1122+C234,:)-4*Dv(N+dd00223,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd12233,:)=in(1,:)
      Dv(N+dd22233,:)=in(2,:)
      Dv(N+dd22333,:)=in(3,:)

      in(1,:)=f1*Dv(N+dd1123,:)-c12sum(:)-4*Dv(N+dd00123,:)
      in(2,:)=f2*Dv(N+dd1123,:)-c12sum(:)-2*Dv(N+dd00113,:)
      in(3,:)=f3*Dv(N+dd1123,:)-c12sum(:)-2*Dv(N+dd00112,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd11123,:)=in(1,:)
      Dv(N+dd11223,:)=half*(Dv(N+dd11223,:)+in(2,:))
      Dv(N+dd11233,:)=half*(Dv(N+dd11233,:)+in(3,:))

      in(1,:)=f1*Dv(N+dd2223,:)-Cv(cc1112+C234,:)+Cv(cc1112+C134,:)
      in(2,:)=f2*Dv(N+dd2223,:)-Cv(cc1112+C234,:)-6*Dv(N+dd00223,:)
      in(3,:)=f3*Dv(N+dd2223,:)-Cv(cc1112+C234,:)-2*Dv(N+dd00222,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd12223,:)=in(1,:)
      Dv(N+dd22223,:)=half*(Dv(N+dd22223,:)+in(2,:))
      Dv(N+dd22233,:)=half*(Dv(N+dd22233,:)+in(3,:))

      in(1,:)=f1*Dv(N+dd2333,:)-Cv(cc1222+C234,:)+Cv(cc1222+C134,:)
      in(2,:)=f2*Dv(N+dd2333,:)-Cv(cc1222+C234,:)-2*Dv(N+dd00333,:)
      in(3,:)=f3*Dv(N+dd2333,:)-Cv(cc1222+C234,:)-6*Dv(N+dd00233,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd12333,:)=in(1,:)
      Dv(N+dd22333,:)=half*(Dv(N+dd22333,:)+in(2,:))
      Dv(N+dd23333,:)=half*(Dv(N+dd23333,:)+in(3,:))

C------end of five index tensors

      if (maxdindex .eq. 5) return

C------six index tensors
      do ep=-2,0
      do j=dd000000,dd003333
      Dv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 60
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd000000,ep)=Dv(N+dd000000,ep)+idp1(j)*half
     &  *(-f1*Dv(N+dd00001,epmj)-f2*Dv(N+dd00002,epmj)
     &    -f3*Dv(N+dd00003,epmj)+Cv(cc0000+C234,epmj)
     &    +two*m1*Dv(N+dd0000,epmj))
      Dv(N+dd001111,ep)=Dv(N+dd001111,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11111,epmj)-f2*Dv(N+dd11112,epmj)-f3*Dv(N+
     &    dd11113,epmj)+Cv(cc0+C234,epmj)+four*Cv(cc1+C234,
     &    epmj)+six*Cv(cc11+C234,epmj)+four*Cv(cc111+C234,
     &    epmj)+Cv(cc1111+C234,epmj)+four*Cv(cc1112+C234,epmj)
     &    +12._dp*Cv(cc112+C234,epmj)+six*Cv(cc1122+C234
     &    ,epmj)+12._dp*Cv(cc12+C234,epmj)+12._dp*Cv(cc122+
     &    C234,epmj)+four*Cv(cc1222+C234,epmj)+four*Cv(cc2+
     &    C234,epmj)+six*Cv(cc22+C234,epmj)+four*Cv(cc222+
     &    C234,epmj)+Cv(cc2222+C234,epmj)
     &    +two*m1*Dv(N+dd1111,epmj))
      Dv(N+dd001112,ep)=Dv(N+dd001112,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11112,epmj)-f2*Dv(N+dd11122,epmj)-f3*Dv(N+
     &    dd11123,epmj)-Cv(cc1+C234,epmj)-three*Cv(cc11+C234
     &    ,epmj)-three*Cv(cc111+C234,epmj)-Cv(cc1111+C234,epmj)
     &    -three*Cv(cc1112+C234,epmj)-six*Cv(cc112+C234,
     &    epmj)-three*Cv(cc1122+C234,epmj)-three*Cv(cc12+C234
     &    ,epmj)-three*Cv(cc122+C234,epmj)-Cv(cc1222+C234,epmj)
     &    +two*m1*Dv(N+dd1112,epmj))
      Dv(N+dd001113,ep)=Dv(N+dd001113,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11113,epmj)-f2*Dv(N+dd11123,epmj)-f3*Dv(N+
     &    dd11133,epmj)-Cv(cc1112+C234,epmj)-three*Cv(cc112+
     &    C234,epmj)-three*Cv(cc1122+C234,epmj)-three*Cv(cc12
     &    +C234,epmj)-six*Cv(cc122+C234,epmj)-three*Cv(
     &    cc1222+C234,epmj)-Cv(cc2+C234,epmj)-three*Cv(cc22
     &    +C234,epmj)-three*Cv(cc222+C234,epmj)-Cv(cc2222+C234,epmj)
     &    +two*m1*Dv(N+dd1113,epmj))
      Dv(N+dd001122,ep)=Dv(N+dd001122,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11122,epmj)-f2*Dv(N+dd11222,epmj)-f3*Dv(N+
     &    dd11223,epmj)+Cv(cc11+C234,epmj)+two*Cv(cc111+
     &    C234,epmj)+Cv(cc1111+C234,epmj)+two*Cv(cc1112+
     &    C234,epmj)+two*Cv(cc112+C234,epmj)
     &    +Cv(cc1122+C234,epmj)+two*m1*Dv(N+dd1122,epmj))
      Dv(N+dd001123,ep)=Dv(N+dd001123,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11123,epmj)-f2*Dv(N+dd11223,epmj)-f3*Dv(N+
     &    dd11233,epmj)+Cv(cc1112+C234,epmj)+two*Cv(cc112+
     &    C234,epmj)+two*Cv(cc1122+C234,epmj)+Cv(cc12+C234
     &    ,epmj)+two*Cv(cc122+C234,epmj)+Cv(cc1222+C234,epmj)
     &    +two*m1*Dv(N+dd1123,epmj))
      Dv(N+dd001133,ep)=Dv(N+dd001133,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11133,epmj)-f2*Dv(N+dd11233,epmj)-f3*Dv(N+
     &    dd11333,epmj)+Cv(cc1122+C234,epmj)+two*Cv(cc122+
     &    C234,epmj)+two*Cv(cc1222+C234,epmj)+Cv(cc22+C234
     &    ,epmj)+two*Cv(cc222+C234,epmj)+Cv(cc2222+C234,epmj)
     &    +two*m1*Dv(N+dd1133,epmj))
      Dv(N+dd001222,ep)=Dv(N+dd001222,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11222,epmj)-f2*Dv(N+dd12222,epmj)-f3*Dv(N+
     &    dd12223,epmj)-Cv(cc111+C234,epmj)-Cv(cc1111+C234,
     &    epmj)-Cv(cc1112+C234,epmj)
     &    +two*m1*Dv(N+dd1222,epmj))
      Dv(N+dd001223,ep)=Dv(N+dd001223,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11223,epmj)-f2*Dv(N+dd12223,epmj)-f3*Dv(N+
     &    dd12233,epmj)-Cv(cc1112+C234,epmj)-Cv(cc112+C234,
     &    epmj)-Cv(cc1122+C234,epmj)
     &    +two*m1*Dv(N+dd1223,epmj))
      Dv(N+dd001233,ep)=Dv(N+dd001233,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11233,epmj)-f2*Dv(N+dd12233,epmj)-f3*Dv(N+
     &    dd12333,epmj)-Cv(cc1122+C234,epmj)-Cv(cc122+C234,
     &    epmj)-Cv(cc1222+C234,epmj)
     &    +two*m1*Dv(N+dd1233,epmj))
      Dv(N+dd001333,ep)=Dv(N+dd001333,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd11333,epmj)-f2*Dv(N+dd12333,epmj)-f3*Dv(N+
     &    dd13333,epmj)-Cv(cc1222+C234,epmj)-Cv(cc222+C234,
     &    epmj)-Cv(cc2222+C234,epmj)
     &    +two*m1*Dv(N+dd1333,epmj))
      Dv(N+dd002222,ep)=Dv(N+dd002222,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd12222,epmj)-f2*Dv(N+dd22222,epmj)-f3*Dv(N+
     &    dd22223,epmj)+Cv(cc1111+C234,epmj)
     &    +two*m1*Dv(N+dd2222,epmj))
      Dv(N+dd002223,ep)=Dv(N+dd002223,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd12223,epmj)-f2*Dv(N+dd22223,epmj)-f3*Dv(N+
     &    dd22233,epmj)+Cv(cc1112+C234,epmj)
     &    +two*m1*Dv(N+dd2223,epmj))
      Dv(N+dd002233,ep)=Dv(N+dd002233,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd12233,epmj)-f2*Dv(N+dd22233,epmj)-f3*Dv(N+
     &    dd22333,epmj)+Cv(cc1122+C234,epmj)
     &    +two*m1*Dv(N+dd2233,epmj))
      Dv(N+dd002333,ep)=Dv(N+dd002333,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd12333,epmj)-f2*Dv(N+dd22333,epmj)-f3*Dv(N+
     &    dd23333,epmj)+Cv(cc1222+C234,epmj)
     &    +two*m1*Dv(N+dd2333,epmj))
      Dv(N+dd003333,ep)=Dv(N+dd003333,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd13333,epmj)-f2*Dv(N+dd23333,epmj)-f3*Dv(N+
     &    dd33333,epmj)+Cv(cc2222+C234,epmj)
     &    +two*m1*Dv(N+dd3333,epmj))
      Dv(N+dd000011,ep)=Dv(N+dd000011,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd00111,epmj)-f2*Dv(N+dd00112,epmj)-f3*Dv(N+
     &    dd00113,epmj)+Cv(cc00+C234,epmj)+two*Cv(cc001+
     &    C234,epmj)+Cv(cc0011+C234,epmj)+two*Cv(cc0012+
     &    C234,epmj)+two*Cv(cc002+C234,epmj)+Cv(cc0022+
     &    C234,epmj)
     &    +two*m1*Dv(N+dd0011,epmj))
      Dv(N+dd000012,ep)=Dv(N+dd000012,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd00112,epmj)-f2*Dv(N+dd00122,epmj)-f3*Dv(N+
     &    dd00123,epmj)-Cv(cc001+C234,epmj)-Cv(cc0011+C234,
     &    epmj)-Cv(cc0012+C234,epmj)
     &    +two*m1*Dv(N+dd0012,epmj))
      Dv(N+dd000013,ep)=Dv(N+dd000013,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd00113,epmj)-f2*Dv(N+dd00123,epmj)-f3*Dv(N+
     &    dd00133,epmj)-Cv(cc0012+C234,epmj)-Cv(cc002+C234,
     &    epmj)-Cv(cc0022+C234,epmj)
     &    +two*m1*Dv(N+dd0013,epmj))
      Dv(N+dd000022,ep)=Dv(N+dd000022,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd00122,epmj)-f2*Dv(N+dd00222,epmj)-f3*Dv(N+
     &    dd00223,epmj)+Cv(cc0011+C234,epmj)
     &    +two*m1*Dv(N+dd0022,epmj))
      Dv(N+dd000023,ep)=Dv(N+dd000023,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd00123,epmj)-f2*Dv(N+dd00223,epmj)-f3*Dv(N+
     &    dd00233,epmj)+Cv(cc0012+C234,epmj)
     &    +two*m1*Dv(N+dd0023,epmj))
      Dv(N+dd000033,ep)=Dv(N+dd000033,ep)+idp1(j)*half*(-f1*Dv(
     &    N+dd00133,epmj)-f2*Dv(N+dd00233,epmj)-f3*Dv(N+
     &    dd00333,epmj)+Cv(cc0022+C234,epmj)
     &    +two*m1*Dv(N+dd0033,epmj))

      enddo
 60   continue
      enddo

C Dv(ppppp)
      in(1,:) = f1*Dv(N+dd11111,:)-10.D0*Dv(N+dd001111,:)+
     & Cv(cc0+C234,:)+5.D0*Cv(cc1+C234,:)+10.D0*Cv(cc11+
     & C234,:)+10.D0*Cv(cc111+C234,:)+5.D0*Cv(cc1111+C234,:)
     & +Cv(cc11111+C234,:)+5.D0*Cv(cc11112+C234,:)+20.D0*
     & Cv(cc1112+C234,:)+10.D0*Cv(cc11122+C234,:)+30.D0*Cv(
     & cc112+C234,:)+30.D0*Cv(cc1122+C234,:)+10.D0*Cv(cc11222
     & +C234,:)+20.D0*Cv(cc12+C234,:)+30.D0*Cv(cc122+C234,
     & :)+20.D0*Cv(cc1222+C234,:)+5.D0*Cv(cc12222+C234,:)+
     & 5.D0*Cv(cc2+C234,:)+10.D0*Cv(cc22+C234,:)+10.D0*Cv(
     & cc222+C234,:)+5.D0*Cv(cc2222+C234,:)+Cv(cc22222+C234,:)

      in(2,:) = f2*Dv(N+dd11111,:)+Cv(cc0+C234,:)+5.D0*Cv(
     & cc1+C234,:)+10.D0*Cv(cc11+C234,:)+10.D0*Cv(cc111+
     & C234,:)+5.D0*Cv(cc1111+C234,:)+Cv(cc11111+C124,:)+
     & Cv(cc11111+C234,:)+5.D0*Cv(cc11112+C234,:)+20.D0*Cv(
     & cc1112+C234,:)+10.D0*Cv(cc11122+C234,:)+30.D0*Cv(cc112
     & +C234,:)+30.D0*Cv(cc1122+C234,:)+10.D0*Cv(cc11222+
     & C234,:)+20.D0*Cv(cc12+C234,:)+30.D0*Cv(cc122+C234,:)
     & +20.D0*Cv(cc1222+C234,:)+5.D0*Cv(cc12222+C234,:)+5.D0
     & *Cv(cc2+C234,:)+10.D0*Cv(cc22+C234,:)+10.D0*Cv(cc222
     & +C234,:)+5.D0*Cv(cc2222+C234,:)+Cv(cc22222+C234,:)

      in(3,:) = f3*Dv(N+dd11111,:)+Cv(cc0+C234,:)+5.D0*Cv(
     & cc1+C234,:)+10.D0*Cv(cc11+C234,:)+10.D0*Cv(cc111+
     & C234,:)+5.D0*Cv(cc1111+C234,:)+Cv(cc11111+C123,:)+
     & Cv(cc11111+C234,:)+5.D0*Cv(cc11112+C234,:)+20.D0*Cv(
     & cc1112+C234,:)+10.D0*Cv(cc11122+C234,:)+30.D0*Cv(cc112
     & +C234,:)+30.D0*Cv(cc1122+C234,:)+10.D0*Cv(cc11222+
     & C234,:)+20.D0*Cv(cc12+C234,:)+30.D0*Cv(cc122+C234,:)
     & +20.D0*Cv(cc1222+C234,:)+5.D0*Cv(cc12222+C234,:)+5.D0
     & *Cv(cc2+C234,:)+10.D0*Cv(cc22+C234,:)+10.D0*Cv(cc222
     & +C234,:)+5.D0*Cv(cc2222+C234,:)+Cv(cc22222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd111111,:)=in(1,:)
      Dv(N+dd111112,:)=in(2,:)
      Dv(N+dd111113,:)=in(3,:)
     
C Dv(ppppk)
      in(1,:) = f1*Dv(N+dd11112,:)-8._dp*Dv(N+dd001112,:)-Cv(
     & cc1+C234,:)-four*Cv(cc11+C234,:)-six*Cv(cc111+C234,
     & :)-four*Cv(cc1111+C234,:)-Cv(cc11111+C234,:)-four*
     & Cv(cc11112+C234,:)-12._dp*Cv(cc1112+C234,:)-six*Cv(
     & cc11122+C234,:)-12._dp*Cv(cc112+C234,:)-12._dp*Cv(cc1122
     & +C234,:)-four*Cv(cc11222+C234,:)-four*Cv(cc12+C234,
     & :)-six*Cv(cc122+C234,:)-four*Cv(cc1222+C234,:)-Cv(
     & cc12222+C234,:)

      in(2,:) = f2*Dv(N+dd11112,:)-two*Dv(N+dd001111,:)-Cv(
     & cc1+C234,:)-four*Cv(cc11+C234,:)-six*Cv(cc111+C234,
     & :)-four*Cv(cc1111+C234,:)-Cv(cc11111+C234,:)-four*
     & Cv(cc11112+C234,:)-12._dp*Cv(cc1112+C234,:)-six*Cv(
     & cc11122+C234,:)-12._dp*Cv(cc112+C234,:)-12._dp*Cv(cc1122
     & +C234,:)-four*Cv(cc11222+C234,:)-four*Cv(cc12+C234,
     & :)-six*Cv(cc122+C234,:)-four*Cv(cc1222+C234,:)-Cv(
     & cc12222+C234,:)

      in(3,:) = f3*Dv(N+dd11112,:)-Cv(cc1+C234,:)-four*Cv(
     & cc11+C234,:)-six*Cv(cc111+C234,:)-four*Cv(cc1111+
     & C234,:)-Cv(cc11111+C234,:)+Cv(cc11112+C123,:)-four*
     & Cv(cc11112+C234,:)-12._dp*Cv(cc1112+C234,:)-six*Cv(
     & cc11122+C234,:)-12._dp*Cv(cc112+C234,:)-12._dp*Cv(cc1122
     & +C234,:)-four*Cv(cc11222+C234,:)-four*Cv(cc12+C234,
     & :)-six*Cv(cc122+C234,:)-four*Cv(cc1222+C234,:)-Cv(
     & cc12222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd111112,:)=in(1,:)
      Dv(N+dd111122,:)=in(2,:)
      Dv(N+dd111123,:)=in(3,:)
     
C Dv(pppkl)
      in(1,:) = f1*Dv(N+dd11123,:)-six*Dv(N+dd001123,:)+Cv(
     & cc11112+C234,:)+three*Cv(cc1112+C234,:)+three*Cv(cc11122
     & +C234,:)+three*Cv(cc112+C234,:)+six*Cv(cc1122+C234,
     & :)+three*Cv(cc11222+C234,:)+Cv(cc12+C234,:)+three*Cv(
     & cc122+C234,:)+three*Cv(cc1222+C234,:)+Cv(cc12222+C234,:)

      in(2,:) = f2*Dv(N+dd11123,:)-two*Dv(N+dd001113,:)+Cv(
     & cc11112+C234,:)+three*Cv(cc1112+C234,:)+three*Cv(cc11122
     & +C234,:)+three*Cv(cc112+C234,:)+six*Cv(cc1122+C234,
     & :)+three*Cv(cc11222+C234,:)+Cv(cc12+C234,:)+three*Cv(
     & cc122+C234,:)+three*Cv(cc1222+C234,:)+Cv(cc12222+C234,:)

      in(3,:) = f3*Dv(N+dd11123,:)-two*Dv(N+dd001112,:)+Cv(
     & cc11112+C234,:)+three*Cv(cc1112+C234,:)+three*Cv(cc11122
     & +C234,:)+three*Cv(cc112+C234,:)+six*Cv(cc1122+C234,
     & :)+three*Cv(cc11222+C234,:)+Cv(cc12+C234,:)+three*Cv(
     & cc122+C234,:)+three*Cv(cc1222+C234,:)+Cv(cc12222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd111123,:)=in(1,:)
      Dv(N+dd111223,:)=in(2,:)
      Dv(N+dd111233,:)=in(3,:)
     
C Dv(pppll)
      in(1,:) = f1*Dv(N+dd11133,:)-six*Dv(N+dd001133,:)+Cv(
     & cc11122+C234,:)+three*Cv(cc1122+C234,:)+three*Cv(cc11222
     & +C234,:)+three*Cv(cc122+C234,:)+six*Cv(cc1222+C234,
     & :)+three*Cv(cc12222+C234,:)+Cv(cc22+C234,:)+three*Cv(
     & cc222+C234,:)+three*Cv(cc2222+C234,:)+Cv(cc22222+C234
     & ,:)

      in(2,:) = f2*Dv(N+dd11133,:)+Cv(cc11122+C124,:)+Cv(
     & cc11122+C234,:)+three*Cv(cc1122+C234,:)+three*Cv(cc11222
     & +C234,:)+three*Cv(cc122+C234,:)+six*Cv(cc1222+C234,
     & :)+three*Cv(cc12222+C234,:)+Cv(cc22+C234,:)+three*Cv(
     & cc222+C234,:)+three*Cv(cc2222+C234,:)+Cv(cc22222+C234
     & ,:)

      in(3,:) = f3*Dv(N+dd11133,:)-four*Dv(N+dd001113,:)+Cv(
     & cc11122+C234,:)+three*Cv(cc1122+C234,:)+three*Cv(cc11222
     & +C234,:)+three*Cv(cc122+C234,:)+six*Cv(cc1222+C234,
     & :)+three*Cv(cc12222+C234,:)+Cv(cc22+C234,:)+three*Cv(
     & cc222+C234,:)+three*Cv(cc2222+C234,:)+Cv(cc22222+C234
     & ,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd111133,:)=in(1,:)
      Dv(N+dd111233,:)=in(2,:)
      Dv(N+dd111333,:)=in(3,:)
     
C Dv(kkkpp)
      in(1,:) = f1*Dv(N+dd11222,:)-four*Dv(N+dd001222,:)-Cv(
     & cc111+C234,:)-two*Cv(cc1111+C234,:)-Cv(cc11111+C234
     & ,:)-two*Cv(cc11112+C234,:)-two*Cv(cc1112+C234,:)-
     & Cv(cc11122+C234,:)

      in(2,:) = f2*Dv(N+dd11222,:)-six*Dv(N+dd001122,:)-Cv(
     & cc111+C234,:)-two*Cv(cc1111+C234,:)-Cv(cc11111+C234
     & ,:)-two*Cv(cc11112+C234,:)-two*Cv(cc1112+C234,:)-
     & Cv(cc11122+C234,:)

      in(3,:) = f3*Dv(N+dd11222,:)-Cv(cc111+C234,:)-two*Cv(
     & cc1111+C234,:)-Cv(cc11111+C234,:)-two*Cv(cc11112+
     & C234,:)-two*Cv(cc1112+C234,:)-Cv(cc11122+C234,:)+
     & Cv(cc11222+C123,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd111222,:)=in(1,:)
      Dv(N+dd112222,:)=in(2,:)
      Dv(N+dd112223,:)=in(3,:)
     
C Dv(lllpp)
      in(1,:) = f1*Dv(N+dd11333,:)-four*Dv(N+dd001333,:)-Cv(
     & cc11222+C234,:)-two*Cv(cc1222+C234,:)-two*Cv(cc12222
     & +C234,:)-Cv(cc222+C234,:)-two*Cv(cc2222+C234,:)-
     & Cv(cc22222+C234,:)

      in(2,:) = f2*Dv(N+dd11333,:)+Cv(cc11222+C124,:)-Cv(
     & cc11222+C234,:)-two*Cv(cc1222+C234,:)-two*Cv(cc12222
     & +C234,:)-Cv(cc222+C234,:)-two*Cv(cc2222+C234,:)-
     & Cv(cc22222+C234,:)

      in(3,:) = f3*Dv(N+dd11333,:)-six*Dv(N+dd001133,:)-Cv(
     & cc11222+C234,:)-two*Cv(cc1222+C234,:)-two*Cv(cc12222
     & +C234,:)-Cv(cc222+C234,:)-two*Cv(cc2222+C234,:)-
     & Cv(cc22222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd111333,:)=in(1,:)
      Dv(N+dd112333,:)=in(2,:)
      Dv(N+dd113333,:)=in(3,:)
     
C Dv(kkllp)
      in(1,:) = f1*Dv(N+dd12233,:)-two*Dv(N+dd002233,:)+Cv(
     & cc11122+C234,:)+Cv(cc1122+C234,:)+Cv(cc11222+C234,:)

      in(2,:) = f2*Dv(N+dd12233,:)-four*Dv(N+dd001233,:)+Cv(
     & cc11122+C234,:)+Cv(cc1122+C234,:)+Cv(cc11222+C234,:)

      in(3,:) = f3*Dv(N+dd12233,:)-four*Dv(N+dd001223,:)+Cv(
     & cc11122+C234,:)+Cv(cc1122+C234,:)+Cv(cc11222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd112233,:)=in(1,:)
      Dv(N+dd122233,:)=in(2,:)
      Dv(N+dd122333,:)=in(3,:)
     
C Dv(kkkkk)
      in(1,:) = f1*Dv(N+dd22222,:)+Cv(cc11111+C134,:)-Cv(
     & cc11111+C234,:)

      in(2,:) = f2*Dv(N+dd22222,:)-10.D0*Dv(N+dd002222,:)-
     & Cv(cc11111+C234,:)

      in(3,:) = f3*Dv(N+dd22222,:)-Cv(cc11111+C234,:)+Cv(
     & cc22222+C123,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd122222,:)=in(1,:)
      Dv(N+dd222222,:)=in(2,:)
      Dv(N+dd222223,:)=in(3,:)
     
C Dv(kkkkl)
      in(1,:) = f1*Dv(N+dd22223,:)+Cv(cc11112+C134,:)-Cv(
     & cc11112+C234,:)
      in(2,:) = f2*Dv(N+dd22223,:)-8._dp*Dv(N+dd002223,:)-Cv(
     & cc11112+C234,:)
      in(3,:) = f3*Dv(N+dd22223,:)-two*Dv(N+dd002222,:)-Cv(
     & cc11112+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd122223,:)=in(1,:)
      Dv(N+dd222223,:)=in(2,:)
      Dv(N+dd222233,:)=in(3,:)
     
C Dv(kkkll)
      in(1,:) = f1*Dv(N+dd22233,:)+Cv(cc11122+C134,:)-Cv(
     & cc11122+C234,:)
      in(2,:) = f2*Dv(N+dd22233,:)-six*Dv(N+dd002233,:)-Cv(
     & cc11122+C234,:)
      in(3,:) = f3*Dv(N+dd22233,:)-four*Dv(N+dd002223,:)-Cv(
     & cc11122+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd122233,:)=in(1,:)
      Dv(N+dd222233,:)=in(2,:)
      Dv(N+dd222333,:)=in(3,:)
     
C Dv(llllk)
      in(1,:) = f1*Dv(N+dd23333,:)+Cv(cc12222+C134,:)-Cv(
     & cc12222+C234,:)
      in(2,:) = f2*Dv(N+dd23333,:)-two*Dv(N+dd003333,:)-Cv(
     & cc12222+C234,:)
      in(3,:) = f3*Dv(N+dd23333,:)-8._dp*Dv(N+dd002333,:)-Cv(
     & cc12222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd123333,:)=in(1,:)
      Dv(N+dd223333,:)=in(2,:)
      Dv(N+dd233333,:)=in(3,:)
     
C Dv(lllll)
      in(1,:) = f1*Dv(N+dd33333,:)+Cv(cc22222+C134,:)-Cv(
     & cc22222+C234,:)
      in(2,:) = f2*Dv(N+dd33333,:)+Cv(cc22222+C124,:)-Cv(
     & cc22222+C234,:)
      in(3,:) = f3*Dv(N+dd33333,:)-10.D0*Dv(N+dd003333,:)-
     & Cv(cc22222+C234,:)
      call pvBackSubst(G, 3, perm, in)
      Dv(N+dd133333,:)=in(1,:)
      Dv(N+dd233333,:)=in(2,:)
      Dv(N+dd333333,:)=in(3,:)

C------end of six index tensors
      if (maxdindex .eq. 6) return

C----seven index
      do ep=-2,0
      do j=dd0000001,dd0033333
      Dv(N+j,ep) = czip
      enddo
      if (ep .eq. -2) goto 70
      do j=0,ep+2
      epmj=ep-j

      Dv(N+dd0000001,ep)=Dv(N+dd0000001,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd000011,epmj)
     &    + f2*Dv(N+dd000012,epmj)
     &    + f3*Dv(N+dd000013,epmj)
     &    - 2*Dv(N+dd00001,epmj)*m1
     &    + Cv(cc00001+C234,epmj)
     &    + Cv(cc00002+C234,epmj)
     &    + Cv(cc0000+C234,epmj))

      Dv(N+dd0000002,ep)=Dv(N+dd0000002,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd000012,epmj)
     &    + f2*Dv(N+dd000022,epmj)
     &    + f3*Dv(N+dd000023,epmj)
     &    - 2*Dv(N+dd00002,epmj)*m1
     &    - Cv(cc00001+C234,epmj))

      Dv(N+dd0000003,ep)=Dv(N+dd0000003,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd000013,epmj)
     &    + f2*Dv(N+dd000023,epmj)
     &    + f3*Dv(N+dd000033,epmj)
     &    - 2*Dv(N+dd00003,epmj)*m1
     &    - Cv(cc00002+C234,epmj))

      Dv(N+dd0000111,ep)=Dv(N+dd0000111,ep)-half*idp2(j)*(
     &    + Cv(cc00+C234,epmj)
     &    + 3*Cv(cc001+C234,epmj)
     &    + 3*Cv(cc002+C234,epmj)
     &    + f1*Dv(N+dd001111,epmj)
     &    + f2*Dv(N+dd001112,epmj)
     &    + f3*Dv(N+dd001113,epmj)
     &    - 2*Dv(N+dd00111,epmj)*m1
     &    + Cv(cc00111+C234,epmj)
     &    + 3*Cv(cc00112+C234,epmj)
     &    + 3*Cv(cc00122+C234,epmj)
     &    + Cv(cc00222+C234,epmj)
     &    + 3*Cv(cc0011+C234,epmj)
     &    + 3*Cv(cc0022+C234,epmj)
     &    + 6*Cv(cc0012+C234,epmj))

      Dv(N+dd0000112,ep)=Dv(N+dd0000112,ep)-half*idp2(j)*(
     &    - Cv(cc001+C234,epmj)
     &    + f1*Dv(N+dd001112,epmj)
     &    + f2*Dv(N+dd001122,epmj)
     &    + f3*Dv(N+dd001123,epmj)
     &    - 2*Dv(N+dd00112,epmj)*m1
     &    - Cv(cc00111+C234,epmj)
     &    - 2*Cv(cc00112+C234,epmj)
     &    - Cv(cc00122+C234,epmj)
     &    - 2*Cv(cc0011+C234,epmj)
     &    - 2*Cv(cc0012+C234,epmj))


      Dv(N+dd0000113,ep)=Dv(N+dd0000113,ep)-half*idp2(j)*(
     &    - Cv(cc002+C234,epmj)
     &    + f1*Dv(N+dd001113,epmj)
     &    + f2*Dv(N+dd001123,epmj)
     &    + f3*Dv(N+dd001133,epmj)
     &    - 2*Dv(N+dd00113,epmj)*m1
     &    - Cv(cc00112+C234,epmj)
     &    - 2*Cv(cc00122+C234,epmj)
     &    - Cv(cc00222+C234,epmj)
     &    - 2*Cv(cc0022+C234,epmj)
     &    - 2*Cv(cc0012+C234,epmj))

      Dv(N+dd0000122,ep)=Dv(N+dd0000122,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001122,epmj)
     &    + f2*Dv(N+dd001222,epmj)
     &    + f3*Dv(N+dd001223,epmj)
     &    - 2*Dv(N+dd00122,epmj)*m1
     &    + Cv(cc00111+C234,epmj)
     &    + Cv(cc00112+C234,epmj)
     &    + Cv(cc0011+C234,epmj))

      Dv(N+dd0000123,ep)=Dv(N+dd0000123,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001123,epmj)
     &    + f2*Dv(N+dd001223,epmj)
     &    + f3*Dv(N+dd001233,epmj)
     &    - 2*Dv(N+dd00123,epmj)*m1
     &    + Cv(cc00112+C234,epmj)
     &    + Cv(cc00122+C234,epmj)
     &    + Cv(cc0012+C234,epmj))

      Dv(N+dd0000133,ep)=Dv(N+dd0000133,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001133,epmj)
     &    + f2*Dv(N+dd001233,epmj)
     &    + f3*Dv(N+dd001333,epmj)
     &    - 2*Dv(N+dd00133,epmj)*m1
     &    + Cv(cc00122+C234,epmj)
     &    + Cv(cc00222+C234,epmj)
     &    + Cv(cc0022+C234,epmj))

      Dv(N+dd0000222,ep)=Dv(N+dd0000222,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001222,epmj)
     &    + f2*Dv(N+dd002222,epmj)
     &    + f3*Dv(N+dd002223,epmj)
     &    - 2*Dv(N+dd00222,epmj)*m1
     &    - Cv(cc00111+C234,epmj))

      Dv(N+dd0000223,ep)=Dv(N+dd0000223,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001223,epmj)
     &    + f2*Dv(N+dd002223,epmj)
     &    + f3*Dv(N+dd002233,epmj)
     &    - 2*Dv(N+dd00223,epmj)*m1
     &    - Cv(cc00112+C234,epmj))

      Dv(N+dd0000233,ep)=Dv(N+dd0000233,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001233,epmj)
     &    + f2*Dv(N+dd002233,epmj)
     &    + f3*Dv(N+dd002333,epmj)
     &    - 2*Dv(N+dd00233,epmj)*m1
     &    - Cv(cc00122+C234,epmj))

      Dv(N+dd0000333,ep)=Dv(N+dd0000333,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd001333,epmj)
     &    + f2*Dv(N+dd002333,epmj)
     &    + f3*Dv(N+dd003333,epmj)
     &    - 2*Dv(N+dd00333,epmj)*m1
     &    - Cv(cc00222+C234,epmj))

      Dv(N+dd0011111,ep)=Dv(N+dd0011111,ep)-half*idp2(j)*(
     &    + Cv(cc0+C234,epmj)
     &    + 5*Cv(cc1+C234,epmj)
     &    + 5*Cv(cc2+C234,epmj)
     &    + 10*Cv(cc11+C234,epmj)
     &    + 20*Cv(cc12+C234,epmj)
     &    + 10*Cv(cc22+C234,epmj)
     &    + 10*Cv(cc111+C234,epmj)
     &    + 30*Cv(cc112+C234,epmj)
     &    + 30*Cv(cc122+C234,epmj)
     &    + 10*Cv(cc222+C234,epmj)
     &    + 5*Cv(cc1111+C234,epmj)
     &    + 20*Cv(cc1112+C234,epmj)
     &    + 30*Cv(cc1122+C234,epmj)
     &    + 20*Cv(cc1222+C234,epmj)
     &    + 5*Cv(cc2222+C234,epmj)
     &    + f1*Dv(N+dd111111,epmj)
     &    + f2*Dv(N+dd111112,epmj)
     &    + f3*Dv(N+dd111113,epmj)
     &    - 2*Dv(N+dd11111,epmj)*m1
     &    + Cv(cc11111+C234,epmj)
     &    + 5*Cv(cc11112+C234,epmj)
     &    + 10*Cv(cc11122+C234,epmj)
     &    + 10*Cv(cc11222+C234,epmj)
     &    + 5*Cv(cc12222+C234,epmj)
     &    + Cv(cc22222+C234,epmj))

      Dv(N+dd0011112,ep)=Dv(N+dd0011112,ep)-half*idp2(j)*(
     &    - Cv(cc1+C234,epmj)
     &    - 4*Cv(cc11+C234,epmj)
     &    - 4*Cv(cc12+C234,epmj)
     &    - 6*Cv(cc111+C234,epmj)
     &    - 12*Cv(cc112+C234,epmj)
     &    - 6*Cv(cc122+C234,epmj)
     &    - 4*Cv(cc1111+C234,epmj)
     &    - 12*Cv(cc1112+C234,epmj)
     &    - 12*Cv(cc1122+C234,epmj)
     &    - 4*Cv(cc1222+C234,epmj)
     &    + f1*Dv(N+dd111112,epmj)
     &    + f2*Dv(N+dd111122,epmj)
     &    + f3*Dv(N+dd111123,epmj)
     &    - 2*Dv(N+dd11112,epmj)*m1
     &    - Cv(cc11111+C234,epmj)
     &    - 4*Cv(cc11112+C234,epmj)
     &    - 6*Cv(cc11122+C234,epmj)
     &    - 4*Cv(cc11222+C234,epmj)
     &    - Cv(cc12222+C234,epmj))

      Dv(N+dd0011113,ep)=Dv(N+dd0011113,ep)-half*idp2(j)*(
     &    - Cv(cc2+C234,epmj)
     &    - 4*Cv(cc12+C234,epmj)
     &    - 4*Cv(cc22+C234,epmj)
     &    - 6*Cv(cc112+C234,epmj)
     &    - 12*Cv(cc122+C234,epmj)
     &    - 6*Cv(cc222+C234,epmj)
     &    - 4*Cv(cc1112+C234,epmj)
     &    - 12*Cv(cc1122+C234,epmj)
     &    - 12*Cv(cc1222+C234,epmj)
     &    - 4*Cv(cc2222+C234,epmj)
     &    + f1*Dv(N+dd111113,epmj)
     &    + f2*Dv(N+dd111123,epmj)
     &    + f3*Dv(N+dd111133,epmj)
     &    - 2*Dv(N+dd11113,epmj)*m1
     &    - Cv(cc11112+C234,epmj)
     &    - 4*Cv(cc11122+C234,epmj)
     &    - 6*Cv(cc11222+C234,epmj)
     &    - 4*Cv(cc12222+C234,epmj)
     &    - Cv(cc22222+C234,epmj))

      Dv(N+dd0011122,ep)=Dv(N+dd0011122,ep)-half*idp2(j)*(
     &   + Cv(cc11+C234,epmj)
     &    + 3*Cv(cc111+C234,epmj)
     &    + 3*Cv(cc112+C234,epmj)
     &    + 3*Cv(cc1111+C234,epmj)
     &    + 6*Cv(cc1112+C234,epmj)
     &    + 3*Cv(cc1122+C234,epmj)
     &    + f1*Dv(N+dd111122,epmj)
     &    + f2*Dv(N+dd111222,epmj)
     &    + f3*Dv(N+dd111223,epmj)
     &    - 2*Dv(N+dd11122,epmj)*m1
     &    + Cv(cc11111+C234,epmj)
     &    + 3*Cv(cc11112+C234,epmj)
     &    + 3*Cv(cc11122+C234,epmj)
     &    + Cv(cc11222+C234,epmj))

      Dv(N+dd0011123,ep)=Dv(N+dd0011123,ep)-half*idp2(j)*(
     &    + Cv(cc12+C234,epmj)
     &    + 3*Cv(cc112+C234,epmj)
     &    + 3*Cv(cc122+C234,epmj)
     &    + 3*Cv(cc1112+C234,epmj)
     &    + 6*Cv(cc1122+C234,epmj)
     &    + 3*Cv(cc1222+C234,epmj)
     &    + f1*Dv(N+dd111123,epmj)
     &    + f2*Dv(N+dd111223,epmj)
     &    + f3*Dv(N+dd111233,epmj)
     &    - 2*Dv(N+dd11123,epmj)*m1
     &    + Cv(cc11112+C234,epmj)
     &    + 3*Cv(cc11122+C234,epmj)
     &    + 3*Cv(cc11222+C234,epmj)
     &    + Cv(cc12222+C234,epmj))

      Dv(N+dd0011133,ep)=Dv(N+dd0011133,ep)-half*idp2(j)*(
     &    + Cv(cc22+C234,epmj)
     &    + 3*Cv(cc122+C234,epmj)
     &    + 3*Cv(cc222+C234,epmj)
     &    + 3*Cv(cc1122+C234,epmj)
     &    + 6*Cv(cc1222+C234,epmj)
     &    + 3*Cv(cc2222+C234,epmj)
     &    + f1*Dv(N+dd111133,epmj)
     &    + f2*Dv(N+dd111233,epmj)
     &    + f3*Dv(N+dd111333,epmj)
     &    - 2*Dv(N+dd11133,epmj)*m1
     &    + Cv(cc11122+C234,epmj)
     &    + 3*Cv(cc11222+C234,epmj)
     &    + 3*Cv(cc12222+C234,epmj)
     &    + Cv(cc22222+C234,epmj))

      Dv(N+dd0011222,ep)=Dv(N+dd0011222,ep)-half*idp2(j)*(
     &    - Cv(cc111+C234,epmj)
     &    - 2*Cv(cc1111+C234,epmj)
     &    - 2*Cv(cc1112+C234,epmj)
     &    + f1*Dv(N+dd111222,epmj)
     &    + f2*Dv(N+dd112222,epmj)
     &    + f3*Dv(N+dd112223,epmj)
     &    - 2*Dv(N+dd11222,epmj)*m1
     &    - Cv(cc11111+C234,epmj)
     &    - 2*Cv(cc11112+C234,epmj)
     &    - Cv(cc11122+C234,epmj))

      Dv(N+dd0011223,ep)=Dv(N+dd0011223,ep)-half*idp2(j)*(
     &    - Cv(cc112+C234,epmj)
     &    - 2*Cv(cc1112+C234,epmj)
     &    - 2*Cv(cc1122+C234,epmj)
     &    + f1*Dv(N+dd111223,epmj)
     &    + f2*Dv(N+dd112223,epmj)
     &    + f3*Dv(N+dd112233,epmj)
     &    - 2*Dv(N+dd11223,epmj)*m1
     &    - Cv(cc11112+C234,epmj)
     &    - 2*Cv(cc11122+C234,epmj)
     &    - Cv(cc11222+C234,epmj))

      Dv(N+dd0011233,ep)=Dv(N+dd0011233,ep)-half*idp2(j)*(
     &    - Cv(cc122+C234,epmj)
     &    - 2*Cv(cc1122+C234,epmj)
     &    - 2*Cv(cc1222+C234,epmj)
     &    + f1*Dv(N+dd111233,epmj)
     &    + f2*Dv(N+dd112233,epmj)
     &    + f3*Dv(N+dd112333,epmj)
     &    - 2*Dv(N+dd11233,epmj)*m1
     &    - Cv(cc11122+C234,epmj)
     &    - 2*Cv(cc11222+C234,epmj)
     &    - Cv(cc12222+C234,epmj))

      Dv(N+dd0011333,ep)=Dv(N+dd0011333,ep)-half*idp2(j)*(
     &    - Cv(cc222+C234,epmj)
     &    - 2*Cv(cc1222+C234,epmj)
     &    - 2*Cv(cc2222+C234,epmj)
     &    + f1*Dv(N+dd111333,epmj)
     &    + f2*Dv(N+dd112333,epmj)
     &    + f3*Dv(N+dd113333,epmj)
     &    - 2*Dv(N+dd11333,epmj)*m1
     &    - Cv(cc11222+C234,epmj)
     &    - 2*Cv(cc12222+C234,epmj)
     &    - Cv(cc22222+C234,epmj))

      Dv(N+dd0012222,ep)=Dv(N+dd0012222,ep)-half*idp2(j)*(
     &    + Cv(cc1111+C234,epmj)
     &    + f1*Dv(N+dd112222,epmj)
     &    + f2*Dv(N+dd122222,epmj)
     &    + f3*Dv(N+dd122223,epmj)
     &    - 2*Dv(N+dd12222,epmj)*m1
     &    + Cv(cc11111+C234,epmj)
     &    + Cv(cc11112+C234,epmj))

      Dv(N+dd0012223,ep)=Dv(N+dd0012223,ep)-half*idp2(j)*(
     &    + Cv(cc1112+C234,epmj)
     &    + f1*Dv(N+dd112223,epmj)
     &    + f2*Dv(N+dd122223,epmj)
     &    + f3*Dv(N+dd122233,epmj)
     &    - 2*Dv(N+dd12223,epmj)*m1
     &    + Cv(cc11112+C234,epmj)
     &    + Cv(cc11122+C234,epmj))

      Dv(N+dd0012233,ep)=Dv(N+dd0012233,ep)-half*idp2(j)*(
     &   + Cv(cc1122+C234,epmj)
     &    + f1*Dv(N+dd112233,epmj)
     &    + f2*Dv(N+dd122233,epmj)
     &    + f3*Dv(N+dd122333,epmj)
     &    - 2*Dv(N+dd12233,epmj)*m1
     &    + Cv(cc11122+C234,epmj)
     &    + Cv(cc11222+C234,epmj))

      Dv(N+dd0012333,ep)=Dv(N+dd0012333,ep)-half*idp2(j)*(
     &    + Cv(cc1222+C234,epmj)
     &    + f1*Dv(N+dd112333,epmj)
     &    + f2*Dv(N+dd122333,epmj)
     &    + f3*Dv(N+dd123333,epmj)
     &    - 2*Dv(N+dd12333,epmj)*m1
     &    + Cv(cc11222+C234,epmj)
     &    + Cv(cc12222+C234,epmj))

      Dv(N+dd0013333,ep)=Dv(N+dd0013333,ep)-half*idp2(j)*(
     &    + Cv(cc2222+C234,epmj)
     &    + f1*Dv(N+dd113333,epmj)
     &    + f2*Dv(N+dd123333,epmj)
     &    + f3*Dv(N+dd133333,epmj)
     &    - 2*Dv(N+dd13333,epmj)*m1
     &    + Cv(cc12222+C234,epmj)
     &    + Cv(cc22222+C234,epmj))

      Dv(N+dd0022222,ep)=Dv(N+dd0022222,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd122222,epmj)
     &    + f2*Dv(N+dd222222,epmj)
     &    + f3*Dv(N+dd222223,epmj)
     &    - 2*Dv(N+dd22222,epmj)*m1
     &    - Cv(cc11111+C234,epmj))

      Dv(N+dd0022223,ep)=Dv(N+dd0022223,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd122223,epmj)
     &    + f2*Dv(N+dd222223,epmj)
     &    + f3*Dv(N+dd222233,epmj)
     &    - 2*Dv(N+dd22223,epmj)*m1
     &    - Cv(cc11112+C234,epmj))

      Dv(N+dd0022233,ep)=Dv(N+dd0022233,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd122233,epmj)
     &    + f2*Dv(N+dd222233,epmj)
     &    + f3*Dv(N+dd222333,epmj)
     &    - 2*Dv(N+dd22233,epmj)*m1
     &    - Cv(cc11122+C234,epmj))

      Dv(N+dd0022333,ep)=Dv(N+dd0022333,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd122333,epmj)
     &    + f2*Dv(N+dd222333,epmj)
     &    + f3*Dv(N+dd223333,epmj)
     &    - 2*Dv(N+dd22333,epmj)*m1
     &    - Cv(cc11222+C234,epmj))

      Dv(N+dd0023333,ep)=Dv(N+dd0023333,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd123333,epmj)
     &    + f2*Dv(N+dd223333,epmj)
     &    + f3*Dv(N+dd233333,epmj)
     &    - 2*Dv(N+dd23333,epmj)*m1
     &    - Cv(cc12222+C234,epmj))

      Dv(N+dd0033333,ep)=Dv(N+dd0033333,ep)-half*idp2(j)*(
     &    + f1*Dv(N+dd133333,epmj)
     &    + f2*Dv(N+dd233333,epmj)
     &    + f3*Dv(N+dd333333,epmj)
     &    - 2*Dv(N+dd33333,epmj)*m1
     &    - Cv(cc22222+C234,epmj))
      enddo
 70   continue
      enddo

      do j=1,Ndd
c      write(66,*) 'D:zx',j,Dv(N+j,-2),Dv(N+j,-1),Dv(N+j,0)
      enddo
      include 'pvD7.f'

c--- to check recursion identities      
c      call Dfill_alt(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N)

      end

