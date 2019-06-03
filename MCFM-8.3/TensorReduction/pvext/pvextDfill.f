      subroutine pvextDfill(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N)
      implicit none
!     N is the offset in the common block
!     p1,p2,p3,p4 are the invariant masses squared of external lines
!     m1,m2,m3,m4 are the masses squared of internal lines
!     Formula based on Ellis,Kunszt,Melnikov,Zanderighi
!     Phys. Report 518 (2012) 141-250
      include 'types.f'
      include 'TRconstants.f'
      include 'TRscale.f'
      include 'pvCnames.f'
      include 'pvDnames.f'
      include 'pvextCv.f'
      include 'pvextDv.f'
      include 'TRmaxindex.f'
      integer::C234,C134,C124,C123,ep,epmj,N,j,perm(3),pvextCcache
      integer,parameter::np=3
      real(dp):: p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,f1,f2,f3
      complex(dp):: G(np,np),in(3,-2:0),trI4
      logical,save:: first=.true.
!$omp threadprivate(first)
      real(dp),save::idm1(0:2),idm2(0:2),idm3(0:2)
!$omp threadprivate(idm1,idm2,idm3)
      include 'cplx.h'
     
      if (first) then
      first=.false.
!--idm1=1/[D-1]
      idm1(0)=one/three
      idm1(1)=idm1(0)*two/three
      idm1(2)=idm1(1)*two/three
!--idm2=1/[D-2]
      idm2(0)=half
      idm2(1)=idm2(0)
      idm2(2)=idm2(1)
!--idm3=1/[D-3]
      idm3(0)=one
      idm3(1)=two*idm3(0)
      idm3(2)=four*idm3(1)
      endif

      f1=m2-m1-p1
      f2=m3-m2-p1p2+p1
      f3=m4-m3-p4+p1p2

!---double up the Gram matrix to remove factors of 1/2 in Eqs.
      G(1,1)=cplx1(two*p1)
      G(2,2)=cplx1(two*p2)
      G(3,3)=cplx1(two*p3)
      G(1,2)=cplx1(p1p2-p1-p2)
      G(2,1)=G(1,2)
      G(1,3)=cplx1(p4+p2-p1p2-p2p3)
      G(3,1)=G(1,3)
      G(2,3)=cplx1(p2p3-p2-p3)
      G(3,2)=G(2,3)
 
!--- initialize integrals
      do ep=-2,0
      do j=1,Ndd
      Dv(N+j,ep)=cplx2(1d5,-1d5)
      enddo
      enddo

!--- Set up relevant triangle pinchings
      C234=pvextCcache(p2,p3,p2p3,m2,m3,m4)
      C134=pvextCcache(p1p2,p3,p4,m1,m3,m4)
      C124=pvextCcache(p1,p2p3,p4,m1,m2,m4)
      C123=pvextCcache(p1,p2,p1p2,m1,m2,m3)

      call XLUDecomp(G, 3, perm)
      
!--- Initialize form-factors
      do j=1,Ndd
      Dv(N+j,:)=cplx1(10000._dp)
      enddo

!--- Initialize box integral
      do ep=-2,0
      Dv(N+dd0,ep) =trI4(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,musq,ep)
      enddo

! d   Eq.(A.29)
      in(1,:)=f1*Dv(N+dd0,:)+Cv(cc0+C134,:)-Cv(cc0+C234,:)
      in(2,:)=f2*Dv(N+dd0,:)+Cv(cc0+C124,:)-Cv(cc0+C134,:)
      in(3,:)=f3*Dv(N+dd0,:)+Cv(cc0+C123,:)-Cv(cc0+C124,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1,:)=in(1,:)
      Dv(N+dd2,:)=in(2,:)
      Dv(N+dd3,:)=in(3,:)

      if (maxdindex .eq. 1) return

C--- two index tensors
!     Eq.(A.61)
      do ep=-2,0
      Dv(N+dd00,ep)=czip
      if (ep .eq. -2) goto 20
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd00,ep)=Dv(N+dd00,ep)
     & +half*idm3(j)*(two*m1*Dv(N+dd0,epmj)+Cv(cc0+C234,epmj)
     & -f1*Dv(N+dd1,epmj)-f2*Dv(N+dd2,epmj)-f3*Dv(N+dd3,epmj))
      enddo
 20   continue
      enddo 

! d1   Eq.(A.30)
      in(1,:)=f1*Dv(N+dd1,:)+Cv(cc0+C234,:)+Cv(cc1+C134,:)
     & -two*Dv(N+dd00,:)
      in(2,:)=f2*Dv(N+dd1,:)-Cv(cc1+C134,:)+Cv(cc1+C124,:)
      in(3,:)=f3*Dv(N+dd1,:)+Cv(cc1+C123,:)-Cv(cc1+C124,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd11,:)=in(1,:)
      Dv(N+dd12,:)=in(2,:)
      Dv(N+dd13,:)=in(3,:)

! d2   Eq.(A.30)
      in(1,:)=f1*Dv(N+dd2,:)-Cv(cc1+C234,:)+Cv(cc1+C134,:)
      in(2,:)=f2*Dv(N+dd2,:)-Cv(cc1+C134,:)+Cv(cc2+C124,:)
     & -two*Dv(N+dd00,:)
      in(3,:)=f3*Dv(N+dd2,:)+Cv(cc2+C123,:)-Cv(cc2+C124,:)
      call pvBackSubst(G,3,perm,in)
!      Dv(N+dd12,:)=in(1,:)
      Dv(N+dd22,:)=in(2,:)
      Dv(N+dd23,:)=in(3,:)

! d3   Eq.(A.30)
      in(1,:)=f1*Dv(N+dd3,:)-Cv(cc2+C234,:)+Cv(cc2+C134,:)
      in(2,:)=f2*Dv(N+dd3,:)-Cv(cc2+C134,:)+Cv(cc2+C124,:)
      in(3,:)=f3*Dv(N+dd3,:)-Cv(cc2+C124,:)-two*Dv(N+dd00,:)
      call pvBackSubst(G,3,perm,in)
!      Dv(N+dd13,:)=in(1,:)
!      Dv(N+dd23,:)=in(2,:)
      Dv(N+dd33,:)=in(3,:)

C --- end of two index tensors

      if (maxdindex .eq. 2) return

C--- three index tensors
!     Eqs.(A.58-A.60)
      do ep=-2,0
      do j=dd001,dd003
      Dv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 30
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd001,ep)=Dv(N+dd001,ep)
     & +half*idm2(j)*(two*m1*Dv(N+dd1,epmj)-Cv(cc0+C234,epmj)
     & -f1*Dv(N+dd11,epmj)-f2*Dv(N+dd12,epmj)-f3*Dv(N+dd13,epmj))
      Dv(N+dd002,ep)=Dv(N+dd002,ep)
     & +half*idm2(j)*(two*m1*Dv(N+dd2,epmj)+Cv(cc1+C234,epmj)
     & -f1*Dv(N+dd12,epmj)-f2*Dv(N+dd22,epmj)-f3*Dv(N+dd23,epmj))
      Dv(N+dd003,ep)=Dv(N+dd003,ep)
     & +half*idm2(j)*(two*m1*Dv(N+dd3,epmj)+Cv(cc2+C234,epmj)
     & -f1*Dv(N+dd13,epmj)-f2*Dv(N+dd23,epmj)-f3*Dv(N+dd33,epmj))
      enddo
 30   continue
      enddo

! d11   Eq.(A.34)
      in(1,:)=f1*Dv(N+dd11,:)-Cv(cc0+C234,:)+Cv(cc11+C134,:)
     & -four*Dv(N+dd001,:)
      in(2,:)=f2*Dv(N+dd11,:)-Cv(cc11+C134,:)+Cv(cc11+C124,:)
      in(3,:)=f3*Dv(N+dd11,:)+Cv(cc11+C123,:)-Cv(cc11+C124,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd111,:)=in(1,:)
      Dv(N+dd112,:)=in(2,:)
      Dv(N+dd113,:)=in(3,:)

! d22   Eq.(A.35)
      in(1,:)=f1*Dv(N+dd22,:)-Cv(cc11+C234,:)+Cv(cc11+C134,:)
      in(2,:)=f2*Dv(N+dd22,:)-Cv(cc11+C134,:)+Cv(cc22+C124,:)
     &-four*Dv(N+dd002,:)
      in(3,:)=f3*Dv(N+dd22,:)+Cv(cc22+C123,:)-Cv(cc22+C124,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd122,:)=in(1,:)
      Dv(N+dd222,:)=in(2,:)
      Dv(N+dd223,:)=in(3,:)

! d33   Eq.(A.36)
      in(1,:)=f1*Dv(N+dd33,:)-Cv(cc22+C234,:)+Cv(cc22+C134,:)
      in(2,:)=f2*Dv(N+dd33,:)-Cv(cc22+C134,:)+Cv(cc22+C124,:)
      in(3,:)=f3*Dv(N+dd33,:)-Cv(cc22+C124,:)-four*Dv(N+dd003,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd133,:)=in(1,:)
      Dv(N+dd233,:)=in(2,:)
      Dv(N+dd333,:)=in(3,:)

! d13   Eq.(A.32)
      in(1,:)=f1*Dv(N+dd13,:)+Cv(cc2+C234,:)+Cv(cc12+C134,:)
     & -two*Dv(N+dd003,:)
      in(2,:)=f2*Dv(N+dd13,:)-Cv(cc12+C134,:)+Cv(cc12+C124,:)
      in(3,:)=f3*Dv(N+dd13,:)-Cv(cc12+C124,:)-two*Dv(N+dd001,:)
      call pvBackSubst(G,3,perm,in)
!      Dv(N+dd113,:)=in(1,:)
      Dv(N+dd123,:)=in(2,:)
!      Dv(N+dd133,:)=in(3,:)

      if (maxdindex .eq. 3) return
      
C--- four index tensors
!     Eqs.(A.51-A.57)
      do ep=-2,0
      do j=dd0000,dd0033
      Dv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 40
      do j=0,ep+2
      epmj=ep-j
      Dv(N+dd0000,ep)=Dv(N+dd0000,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd00,epmj)+Cv(cc00+C234,epmj)
     & -f1*Dv(N+dd001,epmj)-f2*Dv(N+dd002,epmj)-f3*Dv(N+dd003,epmj))
      Dv(N+dd0011,ep)=Dv(N+dd0011,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd11,epmj)+Cv(cc0+C234,epmj)
     & -f1*Dv(N+dd111,epmj)-f2*Dv(N+dd112,epmj)-f3*Dv(N+dd113,epmj))
      Dv(N+dd0012,ep)=Dv(N+dd0012,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd12,epmj)-Cv(cc1+C234,epmj)
     & -f1*Dv(N+dd112,epmj)-f2*Dv(N+dd122,epmj)-f3*Dv(N+dd123,epmj))
      Dv(N+dd0013,ep)=Dv(N+dd0013,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd13,epmj)-Cv(cc2+C234,epmj)
     & -f1*Dv(N+dd113,epmj)-f2*Dv(N+dd123,epmj)-f3*Dv(N+dd133,epmj))
      Dv(N+dd0022,ep)=Dv(N+dd0022,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd22,epmj)+Cv(cc11+C234,epmj)
     & -f1*Dv(N+dd122,epmj)-f2*Dv(N+dd222,epmj)-f3*Dv(N+dd223,epmj))
      Dv(N+dd0023,ep)=Dv(N+dd0023,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd23,epmj)+Cv(cc12+C234,epmj)
     & -f1*Dv(N+dd123,epmj)-f2*Dv(N+dd223,epmj)-f3*Dv(N+dd233,epmj))
      Dv(N+dd0033,ep)=Dv(N+dd0033,ep)
     & +half*idm1(j)*(two*m1*Dv(N+dd33,epmj)+Cv(cc22+C234,epmj)
     & -f1*Dv(N+dd133,epmj)-f2*Dv(N+dd233,epmj)-f3*Dv(N+dd333,epmj))
      enddo
 40   continue
      enddo

! d111   Eq.(A.38)
      in(1,:)=f1*Dv(N+dd111,:)+Cv(cc111+C134,:)+Cv(cc0+C234,:)
     & -six*Dv(N+dd0011,:)
      in(2,:)=f2*Dv(N+dd111,:)-Cv(cc111+C134,:)+Cv(cc111+C124,:)
      in(3,:)=f3*Dv(N+dd111,:)-Cv(cc111+C124,:)+Cv(cc111+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1111,:)=in(1,:)
      Dv(N+dd1112,:)=in(2,:)
      Dv(N+dd1113,:)=in(3,:)

! d113   Eq.(A.42)
      in(1,:)=f1*Dv(N+dd113,:)+Cv(cc112+C134,:)-Cv(cc2+C234,:)
     & -four*Dv(N+dd0013,:)
      in(2,:)=f2*Dv(N+dd113,:)-Cv(cc112+C134,:)+Cv(cc112+C124,:)
      in(3,:)=f3*Dv(N+dd113,:)-Cv(cc112+C124,:)-two*Dv(N+dd0011,:)
      call pvBackSubst(G,3,perm,in)
!      Dv(N+dd1113,:)=in(1,:)
      Dv(N+dd1123,:)=in(2,:)
      Dv(N+dd1133,:)=in(3,:)

! d122   Eq.(A.43)
      in(1,:)=f1*Dv(N+dd122,:)+Cv(cc111+C134,:)+Cv(cc11+C234,:)
     & -two*Dv(N+dd0022,:)
      in(2,:)=f2*Dv(N+dd122,:)-Cv(cc111+C134,:)+Cv(cc122+C124,:)
     & -four*Dv(N+dd0012,:)
      in(3,:)=f3*Dv(N+dd122,:)-Cv(cc122+C124,:)+Cv(cc122+C123,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1122,:)=in(1,:)
      Dv(N+dd1222,:)=in(2,:)
      Dv(N+dd1223,:)=in(3,:)

!d222   Eq.(A.39)
      in(1,:)=f1*Dv(N+dd222,:)+Cv(cc111+C134,:)-Cv(cc111+C234,:)
      in(2,:)=f2*Dv(N+dd222,:)-Cv(cc111+C134,:)+Cv(cc222+C124,:)
     & -six*Dv(N+dd0022,:)
      in(3,:)=f3*Dv(N+dd222,:)-Cv(cc222+C124,:)+Cv(cc222+C123,:)
      call pvBackSubst(G,3,perm,in)
!      Dv(N+dd1222,:)=in(1,:)
      Dv(N+dd2222,:)=in(2,:)
      Dv(N+dd2223,:)=in(3,:)

!d233   Eq.(A.46)
      in(1,:)=f1*Dv(N+dd233,:)+Cv(cc122+C134,:)-Cv(cc122+C234,:)
      in(2,:)=f2*Dv(N+dd233,:)-Cv(cc122+C134,:)+Cv(cc222+C124,:)
     & -two*Dv(N+dd0033,:)
      in(3,:)=f3*Dv(N+dd233,:)-Cv(cc222+C124,:)-four*Dv(N+dd0023,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1233,:)=in(1,:)
      Dv(N+dd2233,:)=in(2,:)
      Dv(N+dd2333,:)=in(3,:)

!d333   Eq.(A.40)
      in(1,:)=f1*Dv(N+dd333,:)+Cv(cc222+C134,:)-Cv(cc222+C234,:)
      in(2,:)=f2*Dv(N+dd333,:)-Cv(cc222+C134,:)+Cv(cc222+C124,:)
      in(3,:)=f3*Dv(N+dd333,:)-Cv(cc222+C124,:)-six*Dv(N+dd0033,:)
      call pvBackSubst(G,3,perm,in)
      Dv(N+dd1333,:)=in(1,:)
!      Dv(N+dd2333,:)=in(2,:)
      Dv(N+dd3333,:)=in(3,:)

      return
      end

