      subroutine pvextCfill(p1,p2,p1p2,m1s,m2s,m3s,N)
      implicit none
C     Calculate the form factors for massless triangle diagrams
C     p1=p1sq,p2=p2sq,p1p2=(p1+p2)^2
C     N is the offset in the common block
C     Formula based on Ellis,Kunszt,Melnikov,Zanderighi
C     Phys. Report 518 (2012) 141-250
      include 'types.f'
      include 'TRconstants.f'
      include 'TRmaxindex.f'
      include 'TRscale.f'
      include 'pvCnames.f'
      include 'pvBnames.f'
      include 'pvextBv.f'
      include 'pvextCv.f'
      integer::B12,B23,B13,ep,epmj,N,j,perm(2),pvextBcache
      integer,parameter::np=2
      complex(dp)::G(np,np),in(2,-2:0),trI3
      real(dp)::p1,p2,p1p2,m1s,m2s,m3s,f1,f2
      real(dp),save::idm1(0:2),idm2(0:2)

      logical,save:: first=.true.
!$omp threadprivate(first,idm1,idm2)
      include 'cplx.h'
      if (first) then
      first=.false.
C--idm1=1/[D-1]
      idm1(0)=one/three
      idm1(1)=idm1(0)*two/three
      idm1(2)=idm1(1)*two/three
C--idm2=1/[D-2]
      idm2(0)=half
      idm2(1)=idm2(0)
      idm2(2)=idm2(1)
      endif

      B12=pvextBcache(p1,m1s,m2s)
      B23=pvextBcache(p2,m2s,m3s)
      B13=pvextBcache(p1p2,m1s,m3s)

      f1=m2s-m1s-p1
      f2=m3s-m2s+p1-p1p2

!---double up the Gram matrix to remove factors of 1/2 in Eqs.
      G(1,1)=cplx1(two*p1)
      G(2,2)=cplx1(two*p2)
      G(1,2)=cplx1(p1p2-p1-p2)
      G(2,1)=G(1,2)

!--- initialize integrals
      do ep=-2,0
      do j=1,Ncc
      Cv(N+j,ep)=cplx2(1d5,-1d5)
      enddo
      enddo

      call XLUDecomp(G, 2, perm)

!Initialize scalar integral
      do ep=-2,0
      Cv(N+cc0,ep)=trI3(p1,p2,p1p2,m1s,m2s,m3s,musq,ep) !trI3 is a switching routine
      enddo

! cz, Eq.(A.11)
      in(1,:)=f1*Cv(N+cc0,:)+Bv(bb0+B13,:)-Bv(bb0+B23,:)
      in(2,:)=f2*Cv(N+cc0,:)+Bv(bb0+B12,:)-Bv(bb0+B13,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc1,:)=in(1,:)
      Cv(N+cc2,:)=in(2,:)
      
      if (maxcindex .eq. 1) return

C---two index form factors
!     Eq.(A.20)
      do ep=-2,0
      Cv(N+cc00,ep)=czip
      if (ep .eq. -2) goto 20
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc00,ep)=Cv(N+cc00,ep)
     & +half*idm2(j)*(two*m1s*Cv(N+cc0,epmj)+Bv(bb0+B23,epmj)
     & -f1*Cv(N+cc1,epmj)-f2*Cv(N+cc2,epmj))
      enddo
 20   continue
      enddo

!c1   Eq.(A.12)
      in(1,:)=f1*Cv(N+cc1,:)+Bv(bb1+B13,:)+Bv(bb0+B23,:)
     & -two*Cv(N+cc00,:)
      in(2,:)=f2*Cv(N+cc1,:)+Bv(bb1+B12,:)-Bv(bb1+B13,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc11,:)=in(1,:)
      Cv(N+cc12,:)=in(2,:)

!c2   Eq.(A.13)
      in(1,:)=f1*Cv(N+cc2,:)+Bv(bb1+B13,:)-Bv(bb1+B23,:)
      in(2,:)=f2*Cv(N+cc2,:)-Bv(bb1+B13,:)-two*Cv(N+cc00,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc12,:)=in(1,:)
      Cv(N+cc22,:)=in(2,:)
       
      if (maxcindex .eq. 2) return

!---three index form factors

!c11   Eq.(A.21,A.22)
      do ep=-2,0
      Cv(N+cc001,ep)=czip
      Cv(N+cc002,ep)=czip
      if (ep .eq. -2) goto 30
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc001,ep)=Cv(N+cc001,ep)
     & +half*idm1(j)*(two*m1s*Cv(N+cc1,epmj)-Bv(bb0+B23,epmj)
     & -f1*Cv(N+cc11,epmj)-f2*Cv(N+cc12,epmj))
      Cv(N+cc002,ep)=Cv(N+cc002,ep)
     & +half*idm1(j)*(two*m1s*Cv(N+cc2,epmj)+Bv(bb1+B23,epmj)
     & -f1*Cv(N+cc12,epmj)-f2*Cv(N+cc22,epmj))
      enddo
 30   continue
      enddo

!c11   Eq.(A.14)
      in(1,:)=f1*Cv(N+cc11,:)+Bv(bb11+B13,:)-Bv(bb0+B23,:)
     & -four*Cv(N+cc001,:)
      in(2,:)=f2*Cv(N+cc11,:)+Bv(bb11+B12,:)-Bv(bb11+B13,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc111,:)=in(1,:)
      Cv(N+cc112,:)=in(2,:)

!c22   Eq.(A.15)
      in(1,:)=f1*Cv(N+cc22,:)+Bv(bb11+B13,:)-Bv(bb11+B23,:)
      in(2,:)=f2*Cv(N+cc22,:)-Bv(bb11+B13,:)-four*Cv(N+cc002,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc122,:)=in(1,:)
      Cv(N+cc222,:)=in(2,:)

!c12   Eq.(A.16)
      in(1,:)=f1*Cv(N+cc12,:)+Bv(bb11+B13,:)+Bv(bb1+B23,:)
     & -two*Cv(N+cc002,:)
      in(2,:)=f2*Cv(N+cc12,:)-Bv(bb11+B13,:)-two*Cv(N+cc001,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc112,:)=in(1,:)
!      Cv(N+cc122,:)=in(2,:)
       
      return
      end
