      subroutine initqqgg(za,zb,nwz)      
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     for a given set of momenta (all outgoing)
c     parton1(-p)+parton2(-q) --> 
c     lepton(l)+antilepton(a)+parton3(k1)+ parton4(k2)
C     calculates the matrix element squared 
C     averaged over initial spins and colours
c     for the various
c     sub-processes and stores them in the common block msq
c-----------------------------------------------------------------

      double precision CoupleC,CoupleV
      double precision xmatpc,xmatpv,tempc,tempv
      double precision PC(-nf:nf,-nf:nf),PV(-nf:nf,-nf:nf)
      double precision msqpc(9),msqpv(9)
      common/msq2/msqpc,msqpv     
      common/qqggcouple/PC,PV     
      include 'ckm.f'
      integer iq_qb_l_a_g_g(6)
      integer iqb_q_l_a_g_g(6)
      integer ig_qb_l_a_qb_g(6)
      integer ig_q_l_a_q_g(6)
      integer iqb_g_l_a_qb_g(6)
      integer iq_g_l_a_q_g(6)
      integer ig_g_l_a_qb_q(6)
      integer j,k,jj,kk,nwz

      data iq_qb_l_a_g_g /1,2,3,4,5,6/
      data iqb_q_l_a_g_g /2,1,3,4,5,6/
      data ig_qb_l_a_qb_g/5,2,3,4,1,6/
      data ig_q_l_a_q_g  /5,1,3,4,2,6/
      data iqb_g_l_a_qb_g/2,5,3,4,1,6/
      data iq_g_l_a_q_g  /1,5,3,4,2,6/
      data ig_g_l_a_qb_q /5,6,3,4,1,2/

case 1 ---  q q
case 2 ---  qbar qbar
case 3 ---  q qb
case 4 ---  qb q
case 5 ---  q g   
case 6 ---  qbar g   
case 7 ---  g qbar
case 8 ---  g q
case 9 ---  g g

case 1 ---- qj qk
      msqpc(1)=0d0
      msqpv(1)=0d0
case 2 --- qj qj or qbarj qbarj
      msqpc(2)=0d0
      msqpv(2)=0d0

case 3 ---  qj qbk
      call xmqqgg(iq_qb_l_a_g_g,za,zb,xmatpc,xmatpv)
      msqpc(3)=xmatpc/(four*xnsq)
      msqpv(3)=xmatpv/(four*xnsq)

case 4 ---   qbj qk
      call xmqqgg(iqb_q_l_a_g_g,za,zb,xmatpc,xmatpv)
      msqpc(4)=xmatpc/(four*xnsq)
      msqpv(4)=xmatpv/(four*xnsq)

case 5 ---   qj g   
      call xmqqgg(iq_g_l_a_q_g,za,zb,xmatpc,xmatpv)
      msqpc(5)=xmatpc/(four*XN*V)
      msqpv(5)=xmatpv/(four*XN*V)

case 6 ---   qbarj g   
      call xmqqgg(iqb_g_l_a_qb_g,za,zb,xmatpc,xmatpv)
      msqpc(6)=xmatpc/(four*XN*V)
      msqpv(6)=xmatpv/(four*XN*V)

case 7 ---   g qbarj
      call xmqqgg(ig_qb_l_a_qb_g,za,zb,xmatpc,xmatpv)
      msqpc(7)=xmatpc/(four*XN*V)
      msqpv(7)=xmatpv/(four*XN*V)

case 8 ---  g qj
      call xmqqgg(ig_q_l_a_q_g,za,zb,xmatpc,xmatpv)
      msqpc(8)=xmatpc/(four*XN*V)
      msqpv(8)=xmatpv/(four*XN*V)

case 9 ---   g g

      call xmqqgg(ig_g_l_a_qb_q,za,zb,xmatpc,xmatpv)
      msqpc(9)=xmatpc/(four*V**2)
      msqpv(9)=xmatpv/(four*V**2)

      do j=-nf,nf
      do k=-nf,nf
C ----case q--qb
      if ((j .gt. 0) .and. (k .lt. 0)) then
      PC(j,k)=CoupleC(j,k,nwz)
      PV(j,k)=CoupleV(j,k,nwz)

C ----case qb--q
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      PC(j,k)=CoupleC(j,k,nwz)
      PV(j,k)=CoupleV(j,k,nwz)
c      write(6,*) 'in initqggg j,k',j,k
c      write(6,*) 'in initqggg PC(j,k)',PC(j,k)
c      write(6,*) 
C ----case g--g
      elseif ((j.eq.0) .and. (k.eq. 0)) then
      tempC=0d0
      tempV=0d0
      do jj=-nf,nf
      do kk=-nf,nf
      tempC=tempC+CoupleC(jj,kk,nwz)
      tempV=tempV+CoupleV(jj,kk,nwz)
      enddo
      enddo
      PC(j,k)=tempC
      PV(j,k)=tempV


C ----case q--g
      elseif (k .eq. 0) then
      tempC=0d0
      tempV=0d0
      do kk=-nf,nf
      tempC=tempC+CoupleC(j,kk,nwz)
      tempV=tempV+CoupleV(j,kk,nwz)
      enddo
      PC(j,k)=tempC
      PV(j,k)=tempV

C ----case g--q
      elseif (j .eq. 0) then
      tempC=0d0
      tempV=0d0
      do jj=-nf,nf
      tempC=tempC+CoupleC(jj,k,nwz)
      tempV=tempV+CoupleV(jj,k,nwz)
      enddo
      PC(j,k)=tempC
      PV(j,k)=tempV
          
      endif      

      enddo
      enddo

      return        
      end






