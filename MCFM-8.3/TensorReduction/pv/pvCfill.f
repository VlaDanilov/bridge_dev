      subroutine pvCfill(p1,p2,p1p2,m1s,m2s,m3s,N)
      implicit none
C     Calculate the form factors for massless triangle diagrams
C     p1=p1sq,p2=p2sq,p1p2=(p1+p2)^2
C     N is the offset in the common block
      include 'types.f'
      include 'TRconstants.f'
      include 'TRscale.f'
      include 'pvCnames.f'
      include 'pvBnames.f'
      include 'pvBv.f'
      include 'pvCv.f'
      include 'pvverbose.f'
      include 'pvrecurflags.f'
      integer:: B12,B23,B13,ep,epmj,N,j,perm(2),pvBcache
      integer,parameter:: np=2
      complex(dp):: G(np,np),in(2,-2:0),trI3,
     & bsum(-2:0),b0sum(-2:0),b1sum(-2:0),b11sum(-2:0),b111sum(-2:0),
     & b1111sum(-2:0),b11111sum(-2:0),
     & b00sum(-2:0),b001sum(-2:0),b0011sum(-2:0),b0000sum(-2:0)
      real(dp)::p1,p2,p1p2,m1s,m2s,m3s,f1,f2
      logical::exceptional
      integer,save:: icall,irecur,irecur2,irecur3,irecur4
      real(dp),save::idp3(0:2),idp2(0:2),idp1(0:2),id(0:2),
     & idm1(0:2),idm2(0:2)
      logical,save:: first=.true.
      logical,save:: scaleset=.false.
!$omp threadprivate(first,idp3,idp2,idp1,id,idm1,idm2)
!$omp threadprivate(icall,irecur,irecur2,irecur3,irecur4)
      include 'cplx.h'
      if (first) then
      first=.false.
C--idp3=1/[D+3]
      idp3(0)=one/7._dp
      idp3(1)=idp3(0)*two/7._dp
      idp3(2)=idp3(1)*two/7._dp
C--idp2=1/[D+2]
      idp2(0)=one/6._dp
      idp2(1)=idp2(0)/three
      idp2(2)=idp2(1)/three
C--idp1=1/[D+1]
      idp1(0)=0.2_dp
      idp1(1)=idp1(0)*0.4_dp
      idp1(2)=idp1(1)*0.4_dp
C--id=1/D
      id(0)=0.25_dp
      id(1)=id(0)*0.5_dp
      id(2)=id(1)*0.5_dp
C--idm1=1/[D-1]
      idm1(0)=one/three
      idm1(1)=idm1(0)*two/three
      idm1(2)=idm1(1)*two/three
C--idm2=1/[D-2]
      idm2(0)=half
      idm2(1)=idm2(0)
      idm2(2)=idm2(1)
c--- variables for statistics reporting
      irecur=0
      irecur2=0
      irecur3=0
      irecur4=0
      icall=0      
c--- print out flags for recursion
c      write(6,*) 'pvCfill recursion flags:'
c      write(6,*) '  doGsing  ',doGsing
c      write(6,*) '  doGYsing ',doGYsing
c      write(6,*) '  doPsing  ',doPsing
c      write(6,*) '  doPFsing ',doPFsing
      endif

c--- statistics accounting and reporting      
      icall=icall+1
      if (pvverbose) then
      if (mod(icall,100000) .eq. 0) then
        write(6,77) icall,
     & 1d2*dfloat(icall-irecur-irecur2-irecur3-irecur4)/dfloat(icall),
     & 1d2*dfloat(irecur)/dfloat(icall),
     & 1d2*dfloat(irecur2)/dfloat(icall),
     & 1d2*dfloat(irecur3)/dfloat(icall),
     & 1d2*dfloat(irecur4)/dfloat(icall)
      endif
      endif
   77 format(' Cfill ',i9,': ',5(f6.2,'% : '))
      
      B12=pvBcache(p1,m1s,m2s)
      B23=pvBcache(p2,m2s,m3s)
      B13=pvBcache(p1p2,m1s,m3s)

      f1=m2s-m1s-p1
      f2=m3s-m1s-p1p2

!---double up the Gram matrix to remove factors of 1/2 in Eqs.
      G(1,1)=cplx1(two*p1)
      G(2,2)=cplx1(two*p1p2)
      G(1,2)=cplx1(p1+p1p2-p2)
      G(2,1)=G(1,2)

C     Y(i,j)=mi^2+mj^2-(q_i-q_j)^2
C     where q_1=0,  q_2=p1,  q_3=p_1+p_2;

c      Y(1,1) = cplx1(two*m1s)
c      Y(1,2) = cplx1(m1s + m2s - p1)
c      Y(2,1) = Y(1,2)
c      Y(1,3) = cplx1(m1s + m3s - p1p2)
c      Y(3,1) = Y(1,3)
c      Y(2,2) = cplx1(two*m2s)
c      Y(2,3) = cplx1(m2s + m3s - p2)
c      Y(3,2) = Y(2,3)
c      Y(3,3) = cplx1(two*m3s)
      
c      if (pvverbose) write(6,*) 'Check triangle Ysing'
c      Ysing=pvGramsing(Y,3)

c-- find maximum entry in Gram matrix
c      Gmax=zip
c      do j=1,2
c      do k=j,2
c      if (abs(G(j,k)) .gt. Gmax) Gmax=abs(G(j,k))
c      enddo
c      enddo
      
c      if (pvverbose) write(6,*) 'Gmax=',Gmax
c      Psing=.false.
c--- criterion for small momenta recursion
c      if (Gmax .lt. weenumber) Psing=.true.
      
c-- find maximum of f1 and f2
c      fmax=max(abs(f1),abs(f2))
      
c      if (pvverbose) write(6,*) 'fmax=',fmax
c      Fsing=.false.
c--- criterion for small momenta and small f(k) recursion
c      if (fmax .lt. weenumber) Fsing=.true.
      
c      write(6,*) 'Gsing,Ysing,Psing,Fsing',Gsing,Ysing,Psing,Fsing
      
      exceptional=.false.
      
      if     (doPFsing) then
c--- for small momenta and small f(k)
        if (pvverbose) then
          write(6,*) 'USING TRIANGLE SMALL MOMENTA AND f(k) RECURSION'
      endif
        call Cfill_recur4(p1,p2,p1p2,m1s,m2s,m3s,N)
      irecur4=irecur4+1
      return
      elseif (doPsing) then
        if (pvverbose) then
        write(6,*) 'USING TRIANGLE SMALL MOMENTA RECURSION'
      endif
        call Cfill_recur3(p1,p2,p1p2,m1s,m2s,m3s,N)
      irecur3=irecur3+1
      return
      elseif (doGYsing) then
c--- for small Gram and small Y
        if (pvverbose) then
          write(6,*) 'USING TRIANGLE SMALL Y AND SMALL G RECURSION'
      endif
        call Cfill_recur2(p1,p2,p1p2,m1s,m2s,m3s,N,exceptional)
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
          write(6,*) 'USING TRIANGLE SMALL G RECURSION'
      endif
        call Cfill_recur (p1,p2,p1p2,m1s,m2s,m3s,N)
      irecur=irecur+1
      return
      endif
c--- otherwise, usual PV is fine


c--- initialize integrals      
      do ep=-2,0
      do j=1,Ncc
      Cv(N+j,ep)=cplx2(1d5,-1d5)
      enddo
      enddo

      call XLUDecomp(G, 2, perm)

C---one index form factors
c'B'Id,Bsum111(P?,K?,m1?,m2?)=MM(B111,P,K,m1,m2)+MM(B1111,P,K,m1,m2);
c'B'Id,Bsum001(P?,K?,m1?,m2?)=MM(B001,P,K,m1,m2)+MM(B0011,P,K,m1,m2);
c'B'Id,Bsum00(P?,K?,m1?,m2?)=MM(B00,P,K,m1,m2)+MM(B001,P,K,m1,m2);
c'B'Id,Bsum11(P?,K?,m1?,m2?)=MM(B11,P,K,m1,m2)+MM(B111,P,K,m1,m2);
c'B'Id,Bsum1(P?,K?,m1?,m2?)=MM(B1,P,K,m1,m2)+MM(B11,P,K,m1,m2);
c'B'Id,Bsum0(P?,K?,m1?,m2?)=MM(B0,P,K,m1,m2)+MM(B1,P,K,m1,m2);

      do ep=-2,0
      Cv(N+cc0,ep)=trI3(p1,p2,p1p2,m1s,m2s,m3s,musq,ep)
      enddo

      bsum(:)=Bv(bb0+B23,:)+Bv(bb1+B23,:)

      b0sum(:)=Bv(bb0+B23,:)+Bv(bb1+B23,:)
      b00sum(:)=Bv(bb00+B23,:)+Bv(bb001+B23,:)
      b001sum(:)=Bv(bb001+B23,:)+Bv(bb0011+B23,:)
      b0011sum(:)=Bv(bb0011+B23,:)+Bv(bb00111+B23,:)

      b0000sum(:)=Bv(bb0000+B23,:)+Bv(bb00001+B23,:)

      b1sum(:)=Bv(bb1+B23,:)+Bv(bb11+B23,:)
      b11sum(:)=Bv(bb11+B23,:)+Bv(bb111+B23,:)
      b111sum(:)=Bv(bb111+B23,:)+Bv(bb1111+B23,:)
      b1111sum(:)=Bv(bb1111+B23,:)+Bv(bb11111+B23,:)
      b11111sum(:)=Bv(bb11111+B23,:)+Bv(bb111111+B23,:)

      in(1,:)=f1*Cv(N+cc0,:)-Bv(bb0+B23,:)+Bv(bb0+B13,:)
      in(2,:)=f2*Cv(N+cc0,:)-Bv(bb0+B23,:)+Bv(bb0+B12,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc1,:)=in(1,:)
      Cv(N+cc2,:)=in(2,:)
      
C---two index form factors
      do ep=-2,0
      Cv(N+cc00,ep)=czip
      if (ep .eq. -2) goto 20
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc00,ep)=Cv(N+cc00,ep)+idm2(j)*(m1s*Cv(N+cc0,epmj)
     & -half*(f1*Cv(N+cc1,epmj)+f2*Cv(N+cc2,epmj)-Bv(bb0+B23,epmj)))
      enddo
 20   continue
      enddo

      in(1,:)=f1*Cv(N+cc1,:)+bsum(:)-two*Cv(N+cc00,:)
      in(2,:)=f2*Cv(N+cc1,:)+bsum(:)+Bv(bb1+B12,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc11,:)=in(1,:)
      Cv(N+cc12,:)=in(2,:)
      
      in(1,:)=f1*Cv(N+cc2,:)-Bv(bb1+B23,:)+Bv(bb1+B13,:)
      in(2,:)=f2*Cv(N+cc2,:)-Bv(bb1+B23,:)-two*Cv(N+cc00,:)
      call pvBackSubst(G,2,perm,in)
!      Cv(N+cc12,:)=half*(Cv(N+cc12,:)+in(1,:))
      Cv(N+cc22,:)=in(2,:)
       
c      if ((maxcindex .eq. 2) .and. (pvRespectmaxcindex)) return

C---three index form factors
      do ep=-2,0
      Cv(N+cc001,ep)=czip
      Cv(N+cc002,ep)=czip
      if (ep .eq. -2) goto 30
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc001,ep)=Cv(N+cc001,ep)+idm1(j)*(m1s*Cv(N+cc1,epmj)
     & -half*(f1*Cv(N+cc11,epmj)+f2*Cv(N+cc12,epmj)+bsum(epmj)))
      Cv(N+cc002,ep)=Cv(N+cc002,ep)+idm1(j)*(m1s*Cv(N+cc2,epmj)
     & -half*(f1*Cv(N+cc12,epmj)+f2*Cv(N+cc22,epmj)-Bv(bb1+B23,epmj)))
      enddo
 30   continue
      enddo

      bsum(:)=bsum(:)+b1sum(:)
C--- bsum is now equal to 
c--- Bv(bb1+B23,ep)+2*Bv(bb1+B23,ep)+Bv(bb11+B23,ep)
      in(1,:)=f1*Cv(N+cc11,:)-bsum(:)-four*Cv(N+cc001,:)
      in(2,:)=f2*Cv(N+cc11,:)-bsum(:)+Bv(bb11+B12,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc111,:)=in(1,:)
      Cv(N+cc112,:)=in(2,:)

      in(1,:)=f1*Cv(N+cc22,:)-Bv(bb11+B23,:)+Bv(bb11+B13,:)
      in(2,:)=f2*Cv(N+cc22,:)-Bv(bb11+B23,:)-four*Cv(N+cc002,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc122,:)=in(1,:)
      Cv(N+cc222,:)=in(2,:)

c      b1sum(ep)=Bv(bb1+B23,ep)+Bv(bb11+B23,ep)
      in(1,:)=f1*Cv(N+cc12,:)+b1sum(:)-two*Cv(N+cc002,:)
      in(2,:)=f2*Cv(N+cc12,:)+b1sum(:)-two*Cv(N+cc001,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc112,:)=in(1,:)
      Cv(N+cc122,:)=in(2,:)


c      if ((maxcindex .eq. 3) .and. (pvRespectmaxcindex)) return

C---four index form factors
      do ep=-2,0
      do j=cc0000,cc0022
      Cv(N+j,ep)=czip
      enddo
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc0000,ep)=Cv(N+cc0000,ep)+id(j)*(
     & +m1s*Cv(N+cc00,epmj)
     & -half*(f1*Cv(N+cc001,epmj)+f2*Cv(N+cc002,epmj)
     & -Bv(bb00+B23,epmj)))
      Cv(N+cc0011,ep)=Cv(N+cc0011,ep)+id(j)*(
     & +m1s*Cv(N+cc11,epmj)
     & -half*(f1*Cv(N+cc111,epmj)+f2*Cv(N+cc112,epmj)
     & -b0sum(epmj)-b1sum(epmj)))
      Cv(N+cc0012,ep)=Cv(N+cc0012,ep)+id(j)*(
     & +m1s*Cv(N+cc12,epmj)
     & -half*(f1*Cv(N+cc112,epmj)+f2*Cv(N+cc122,epmj)
     & +b1sum(epmj)))
      Cv(N+cc0022,ep)=Cv(N+cc0022,ep)+id(j)*(
     & +m1s*Cv(N+cc22,epmj)
     & -half*(f1*Cv(N+cc122,epmj)+f2*Cv(N+cc222,epmj)
     & -Bv(bb11+B23,epmj)))

      enddo
c      write(66,*) 'ox',cc0000,ep,Cv(N+cc0000,ep)
c      write(66,*) 'ox',cc0011,ep,Cv(N+cc0011,ep)
c      write(66,*) 'ox',cc0012,ep,Cv(N+cc0012,ep)
c      write(66,*) 'ox',cc0022,ep,Cv(N+cc0022,ep)
      enddo

      bsum(:)=bsum(:)+b1sum(:)+b11sum(:)
C--- bsum is now equal to 
c--- Bv(bb1+B23,ep)+3*Bv(bb1+B23,ep)+3*Bv(bb11+B23,ep)+Bv(bb111+B23,ep)

      in(1,:)=f1*Cv(N+cc111,:)+bsum(:)-six*Cv(N+cc0011,:)
      in(2,:)=f2*Cv(N+cc111,:)+bsum(:)+Bv(bb111+B12,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc1111,:)=in(1,:)
      Cv(N+cc1112,:)=in(2,:)

      in(1,:)=f1*Cv(N+cc112,:)-b1sum(:)-b11sum(:)
     & -four*Cv(N+cc0012,:)
      in(2,:)=f2*Cv(N+cc112,:)-b1sum(:)-b11sum(:)
     & -two*Cv(N+cc0011,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc1112,:)=in(1,:)
      Cv(N+cc1122,:)=in(2,:)

      in(1,:)=f1*Cv(N+cc222,:)-Bv(bb111+B23,:)+Bv(bb111+B13,:)
      in(2,:)=f2*Cv(N+cc222,:)-Bv(bb111+B23,:)-six*Cv(N+cc0022,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc1222,:)=in(1,:)
      Cv(N+cc2222,:)=in(2,:)

      in(1,:)=f1*Cv(N+cc122,:)+b11sum(:)-two*Cv(N+cc0022,:)
      in(2,:)=f2*Cv(N+cc122,:)+b11sum(:)-four*Cv(N+cc0012,:)
      call pvBackSubst(G,2,perm,in)
      Cv(N+cc1122,:)=in(1,:)
      Cv(N+cc1222,:)=half*(Cv(N+cc1222,:)+in(2,:))

c      if ((maxcindex .eq. 4) .and. (pvRespectmaxcindex)) return

C---five index form factors
      do ep=-2,0
      do j=cc00001,cc00222
      Cv(N+j,ep)=czip
      enddo
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc00001,ep)=Cv(N+cc00001,ep)+idp1(j)*(
     & +m1s*Cv(N+cc001,epmj)
     & -half*(f1*Cv(N+cc0011,epmj)+f2*Cv(N+cc0012,epmj)
     & +b00sum(epmj)))
      Cv(N+cc00002,ep)=Cv(N+cc00002,ep)+idp1(j)*(
     & +m1s*Cv(N+cc002,epmj)
     & -half*(f1*Cv(N+cc0012,epmj)+f2*Cv(N+cc0022,epmj)
     & -Bv(bb001+B23,epmj)))
      Cv(N+cc00111,ep)=Cv(N+cc00111,ep)+idp1(j)*(
     & +m1s*Cv(N+cc111,epmj)
     & -half*(f1*Cv(N+cc1111,epmj)+f2*Cv(N+cc1112,epmj)
     & +b0sum(epmj)+two*b1sum(epmj)+b11sum(epmj)))
      Cv(N+cc00112,ep)=Cv(N+cc00112,ep)+idp1(j)*(
     & +m1s*Cv(N+cc112,epmj)
     & -half*(f1*Cv(N+cc1112,epmj)+f2*Cv(N+cc1122,epmj)
     & -b1sum(epmj)-b11sum(epmj)))
      Cv(N+cc00122,ep)=Cv(N+cc00122,ep)+idp1(j)*(
     & +m1s*Cv(N+cc122,epmj)
     & -half*(f1*Cv(N+cc1122,epmj)+f2*Cv(N+cc1222,epmj)
     & +b11sum(epmj)))
      Cv(N+cc00222,ep)=Cv(N+cc00222,ep)+idp1(j)*(
     & +m1s*Cv(N+cc222,epmj)
     & -half*(f1*Cv(N+cc1222,epmj)+f2*Cv(N+cc2222,epmj)
     & -Bv(bb111+B23,epmj)))
      enddo
      enddo

C Cv(pppp
      in(1,:) = f1*Cv(N+cc1111,:)-8.D0*Cv(N+cc00111,:)- Bv(bb0+B23,
     & :)-4.D0*Bv(bb1+B23,:)-6.D0*Bv(bb11+B23,:)-4.D0*Bv(
     & bb111+B23,:)-Bv(bb1111+B23,:)
      in(2,:) = f2*Cv(N+cc1111,:)-Bv(bb0+B23,:)-4.D0*Bv(bb1+
     & B23,:)-6.D0*Bv(bb11+B23,:)-4.D0*Bv(bb111+B23,:)+Bv(
     & bb1111+B12,:)-Bv(bb1111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc11111,:)=in(1,:)
      Cv(N+cc11112,:)=in(2,:)
     
C Cv(pppk
      in(1,:) = f1*Cv(N+cc1112,:)-6.D0*Cv(N+cc00112,:)+ Bv(bb1+B23,
     & :)+3.D0*Bv(bb11+B23,:)+3.D0*Bv(bb111+B23,:)+Bv(
     & bb1111+B23,:)
      in(2,:) = f2*Cv(N+cc1112,:)-2.D0*Cv(N+cc00111,:) +Bv(bb1+B23,
     & :)+3.D0*Bv(bb11+B23,:)+3.D0*Bv(bb111+B23,:)+Bv(
     & bb1111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc11112,:)=in(1,:)
      Cv(N+cc11122,:)=in(2,:)
     
C Cv(pkkk
      in(1,:) = f1*Cv(N+cc1222,:)-2.D0*Cv(N+cc00222,:)+Bv(bb111+
     & B23,:)+Bv(bb1111+B23,:)
      in(2,:) = f2*Cv(N+cc1222,:)-6.D0*Cv(N+cc00122,:)+Bv(bb111+
     & B23,:)+Bv(bb1111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc11222,:)=in(1,:)
      Cv(N+cc12222,:)=in(2,:)
     
C Cv(kkkk
      in(1,:) = f1*Cv(N+cc2222,:)+Bv(bb1111+B13,:)-Bv(bb1111+B23,:)
      in(2,:) = f2*Cv(N+cc2222,:)-8.D0*Cv(N+cc00222,:)
     & -Bv(bb1111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc12222,:)=in(1,:)
      Cv(N+cc22222,:)=in(2,:)

c      if ((maxcindex .eq. 5) .and. (pvRespectmaxcindex)) return

C---six index form factors
      do ep=-2,0
      do j=cc000000,cc002222
      Cv(N+j,:)=czip
      enddo
      if (ep .eq. -2) goto 60
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc000000,ep)=Cv(N+cc000000,ep)+idp2(j)*(
     & +m1s*Cv(N+cc0000,epmj)
     & -half*(f1*Cv(N+cc00001,epmj)+f2*Cv(N+cc00002,epmj)
     & -Bv(bb0000+B23,epmj)))

      Cv(N+cc000011,ep)=Cv(N+cc000011,ep)+idp2(j)*(
     & +m1s*Cv(N+cc0011,epmj)
     & -half*(f1*Cv(N+cc00111,epmj)+f2*Cv(N+cc00112,epmj)
     & -B00sum(epmj)-B001sum(epmj)))

      Cv(N+cc000012,ep)=Cv(N+cc000012,ep)+idp2(j)*(
     & +m1s*Cv(N+cc0012,epmj)
     & -half*(f1*Cv(N+cc00112,epmj)+f2*Cv(N+cc00122,epmj)
     & +B001sum(epmj)))

      Cv(N+cc000022,ep)=Cv(N+cc000022,ep)+idp2(j)*(
     & +m1s*Cv(N+cc0022,epmj)
     & -half*(f1*Cv(N+cc00122,epmj)+ f2*Cv(N+cc00222,epmj)
     & -Bv(bb0011+B23,epmj)))

      Cv(N+cc001111,ep)=Cv(N+cc001111,ep)+idp2(j)*(
     & +m1s*Cv(N+cc1111,epmj)
     & -half*(f1*Cv(N+cc11111,epmj)+f2*Cv(N+cc11112,epmj)
     & -B111sum(epmj)-B0sum(epmj)-three*B1sum(epmj)-three*B11sum(epmj)))

      Cv(N+cc001112,ep)=Cv(N+cc001112,ep)+idp2(j)*(
     & +m1s*Cv(N+cc1112,epmj)
     & -half*(f1*Cv(N+cc11112,epmj)+f2*Cv(N+cc11122,epmj)
     & +B1sum(epmj)+two*B11sum(epmj)+B111sum(epmj)))
      Cv(N+cc001122,ep)=Cv(N+cc001122,ep)+idp2(j)*(
     & +m1s*Cv(N+cc1122,epmj)
     & -half*(f1*Cv(N+cc11122,epmj)+f2*Cv(N+cc11222,epmj)
     & -B11sum(epmj)-B111sum(epmj))) 
      Cv(N+cc001222,ep)=Cv(N+cc001222,ep)+idp2(j)*(
     & +m1s*Cv(N+cc1222,epmj)
     & -half*(f1*Cv(N+cc11222,epmj)+f2*Cv(N+cc12222,epmj)
     & +B111sum(epmj)))

      Cv(N+cc002222,ep)=Cv(N+cc002222,ep)+idp2(j)*(
     & +m1s*Cv(N+cc2222,epmj)
     & -half*(f1*Cv(N+cc12222,epmj)+f2*Cv(N+cc22222,epmj)
     & -Bv(bb1111+B23,epmj)))
      enddo
 60   continue
      enddo

      
C    (Cv(N+cc000011,Cv(N+cc000012,zzzzp)
      in(1,:)=
     &  + f1*Cv(N+cc00001,:)
     &  -2*Cv(N+cc000000,:)
     &  + B0000sum(:)
      in(2,:)=
     &  + f2*Cv(N+cc00001,:)
     &  + Bv(bb00001+B12,:)
     &  + B0000sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc000011,:)=in(1,:)
      Cv(N+cc000012,:)=in(2,:)
     
C    (Cv(N+cc000012,Cv(N+cc000022,zzzzk)
      in(1,:)=
     &  + f1*Cv(N+cc00002,:)
     &  + Bv(bb00001+B13,:)
     &  - Bv(bb00001+B23,:)
      in(2,:)=
     &  + f2*Cv(N+cc00002,:)
     &  -2*Cv(N+cc000000,:)
     &  - Bv(bb00001+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc000012,:)=in(1,:)
      Cv(N+cc000022,:)=in(2,:)

C    (Cv(N+cc001111,Cv(N+cc001112,zzppp)
      in(1,:)=
     &  + f1*Cv(N+cc00111,:)-6*Cv(N+cc000011,:)
     &  + B0011sum(:)+ two*B001sum(:)+ B00sum(:)
      in(2,:)=
     & + f2*Cv(N+cc00111,:)+ Bv(bb00111+B12,:)
     &  + B0011sum(:)+ two*B001sum(:)+ B00sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc001111,:)=in(1,:)
      Cv(N+cc001112,:)=in(2,:)
     
C    (Cv(N+cc001112,Cv(N+cc001122,zzppk)
      in(1,:)=
     &  + f1*Cv(N+cc00112,:)-4*Cv(N+cc000012,:)
     &   -B001sum(:)-B0011sum(:)
      in(2,:)=
     &  + f2*Cv(N+cc00112,:)
     &  -2*Cv(N+cc000011,:)
     &   -B001sum(:)-B0011sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc001112,:)=in(1,:)
      Cv(N+cc001122,:)=in(2,:)
     
C    (Cv(N+cc001122,Cv(N+cc001222,zzpkk)
      in(1,:)=
     &  + f1*Cv(N+cc00122,:)-2*Cv(N+cc000022,:)
     &  + B0011sum(:)
      in(2,:)=
     &  + f2*Cv(N+cc00122,:)-4*Cv(N+cc000012,:)
     &  + B0011sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc001122,:)=in(1,:)
      Cv(N+cc001222,:)=in(2,:)
     
C    (Cv(N+cc001222,Cv(N+cc002222,zzkkk)
      in(1,:)=
     &  + f1*Cv(N+cc00222,:)
     &  + Bv(bb00111+B13,:)- Bv(bb00111+B23,:)
      in(2,:)=
     &  + f2*Cv(N+cc00222,:)-6*Cv(N+cc000022,:)
     &  - Bv(bb00111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc001222,:)=in(1,:)
      Cv(N+cc002222,:)=in(2,:)
     


C    (Cv(N+cc111111,Cv(N+cc111112,ppppp)
      in(1,:)=
     &  +f1*Cv(N+cc11111,:)-10*Cv(N+cc001111,:)
     &  +B0sum(:)+4*B1sum(:)+4*B111sum(:)
     &  +6*B11sum(:)+B1111sum(:)
      in(2,:)=
     &  +f2*Cv(N+cc11111,:)+Bv(bb11111+B12,:)
     &  +B0sum(:)+4*B1sum(:)
     &  +4*B111sum(:)+6*B11sum(:)+B1111sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc111111,:)=in(1,:)
      Cv(N+cc111112,:)=in(2,:)
     
C    (Cv(N+cc111112,Cv(N+cc111122,ppppk)
      in(1,:)=
     &   +f1*Cv(N+cc11112,:)-8*Cv(N+cc001112,:)
     &  -B1sum(:)-3*B11sum(:)-3*B111sum(:)-B1111sum(:)
      in(2,:)=
     &  + f2*Cv(N+cc11112,:)-2*Cv(N+cc001111,:)
     &  -B1sum(:)-3*B11sum(:)-3*B111sum(:)-B1111sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc111112,:)=in(1,:)
      Cv(N+cc111122,:)=in(2,:)
     
C    (Cv(N+cc111122,Cv(N+cc111222,pppkk)
      in(1,:)=
     &  + f1*Cv(N+cc11122,:)-6*Cv(N+cc001122,:)
     &  +B11sum(:)+2*B111sum(:)+B1111sum(:)
      in(2,:)=
     &  + f2*Cv(N+cc11122,:)-4*Cv(N+cc001112,:)
     &  +B11sum(:)+2*B111sum(:)+B1111sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc111122,:)=in(1,:)
      Cv(N+cc111222,:)=in(2,:)
     
C    (Cv(N+cc111222,Cv(N+cc112222,ppkkk)
      in(1,:)=
     &  + f1*Cv(N+cc11222,:)-4*Cv(N+cc001222,:)
     &  -B111sum(:)-B1111sum(:)
      in(2,:)=
     &  + f2*Cv(N+cc11222,:)-6*Cv(N+cc001122,:)
     &  -B111sum(:)-B1111sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc111222,:)=in(1,:)
      Cv(N+cc112222,:)=in(2,:)

C    (Cv(N+cc112222,Cv(N+cc122222,pkkkk)
      in(1,:)=f1*Cv(N+cc12222,:)-2*Cv(N+cc002222,:)
     &  + B1111sum(:)

      in(2,:)=f2*Cv(N+cc12222,:)-8*Cv(N+cc001222,:)
     &  + B1111sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc112222,:)=in(1,:)
      Cv(N+cc122222,:)=in(2,:)

     
C    (Cv(N+cc122222,Cv(N+cc222222,kkkkk)
      in(1,:)=
     &  + f1*Cv(N+cc22222,:)
     &  + Bv(bb11111+B13,:)-Bv(bb11111+B23,:)
      in(2,:)=
     &  + f2*Cv(N+cc22222,:)-10* Cv(N+cc002222,:)
     &  - Bv(bb11111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc122222,:)=in(1,:)
      Cv(N+cc222222,:)=in(2,:)

c      if ((maxcindex .eq. 6) .and. (pvRespectmaxcindex)) return

C---seven index form factors
      do ep=-2,0
      do j=cc0000001,cc0022222
      Cv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 61
      do j=0,ep+2
      epmj=ep-j

      Cv(N+cc0011111,ep)=Cv(N+cc0011111,ep)-idp3(j)*(
     & +half*b1111sum(epmj)
     & +2*b111sum(epmj)
     & +3*b11sum(epmj)
     & +2*b1sum(epmj)
     & +half*b0sum(epmj)
     & +half*f1*Cv(N+cc111111,epmj)
     & +half*f2*Cv(N+cc111112,epmj)
     & -Cv(N+cc11111,epmj)*m1s)


      Cv(N+cc0011112,ep)=Cv(N+cc0011112,ep)-idp3(j)*(
     & -half*b1sum(epmj)
     & -three*half*b11sum(epmj)
     & -three*half*b111sum(epmj)
     & -half*b1111sum(epmj)
     & +half*f1*Cv(N+cc111112,epmj)
     & +half*f2*Cv(N+cc111122,epmj)
     & -Cv(N+cc11112,epmj)*m1s)


      Cv(N+cc0012222,ep)=Cv(N+cc0012222,ep)-idp3(j)*(
     & + half*Bv(bb1111+B23,epmj)
     & + half*Bv(bb11111+B23,epmj)
     & + half*f1*Cv(N+cc112222,epmj)
     & + half*f2*Cv(N+cc122222,epmj)
     & - Cv(N+cc12222,epmj)*m1s)


      Cv(N+cc0022222,ep)=Cv(N+cc0022222,ep)-idp3(j)*(
     & - half*Bv(bb11111+B23,epmj)
     & + half*f1*Cv(N+cc122222,epmj)
     & + half*f2*Cv(N+cc222222,epmj)
     & - Cv(N+cc22222,epmj)*m1s)


      Cv(N+cc0011222,ep)=Cv(N+cc0011222,ep)-idp3(j)*(
     & - half*Bv(bb111+B23,epmj)
     & - Bv(bb1111+B23,epmj)
     & - half*Bv(bb11111+B23,epmj)
     & + half*f1*Cv(N+cc111222,epmj)
     & + half*f2*Cv(N+cc112222,epmj)
     & - Cv(N+cc11222,epmj)*m1s)


      Cv(N+cc0011122,ep)=Cv(N+cc0011122,ep)-idp3(j)*(
     & + half*Bv(bb11+B23,epmj)
     & + three*half*Bv(bb111+B23,epmj)
     & + three*half*Bv(bb1111+B23,epmj)
     & + half*Bv(bb11111+B23,epmj)
     & + half*f1*Cv(N+cc111122,epmj)
     & + half*f2*Cv(N+cc111222,epmj)
     & - Cv(N+cc11122,epmj)*m1s)


      Cv(N+cc0000001,ep)=Cv(N+cc0000001,ep)-idp3(j)*(
     & + half*Bv(bb0000+B23,epmj)
     & + half*Bv(bb00001+B23,epmj)
     & + half*f1*Cv(N+cc000011,epmj)
     & + half*f2*Cv(N+cc000012,epmj)
     & - Cv(N+cc00001,epmj)*m1s)

      Cv(N+cc0000002,ep)=Cv(N+cc0000002,ep)-idp3(j)*(
     & - half*Bv(bb00001+B23,epmj)
     & + half*f1*Cv(N+cc000012,epmj)
     & + half*f2*Cv(N+cc000022,epmj)
     & - Cv(N+cc00002,epmj)*m1s)

      Cv(N+cc0000111,ep)=Cv(N+cc0000111,ep)-idp3(j)*(
     & + half*Bv(bb00+B23,epmj)
     & + three*half*Bv(bb001+B23,epmj)
     & + three*half*Bv(bb0011+B23,epmj)
     & + half*Bv(bb00111+B23,epmj)
     & + half*f1*Cv(N+cc001111,epmj)
     & + half*f2*Cv(N+cc001112,epmj)
     & - Cv(N+cc00111,epmj)*m1s)

      Cv(N+cc0000112,ep)=Cv(N+cc0000112,ep)-idp3(j)*(
     & - half*Bv(bb001+B23,epmj)
     & - Bv(bb0011+B23,epmj)
     & - half*Bv(bb00111+B23,epmj)
     & + half*f1*Cv(N+cc001112,epmj)
     & + half*f2*Cv(N+cc001122,epmj)
     & - Cv(N+cc00112,epmj)*m1s)

      Cv(N+cc0000122,ep)=Cv(N+cc0000122,ep)-idp3(j)*(
     & + half*Bv(bb0011+B23,epmj)
     & + half*Bv(bb00111+B23,epmj)
     & + half*f1*Cv(N+cc001122,epmj)
     & + half*f2*Cv(N+cc001222,epmj)
     & - Cv(N+cc00122,epmj)*m1s)

      Cv(N+cc0000222,ep)=Cv(N+cc0000222,ep)-idp3(j)*(
     & - half*Bv(bb00111+B23,epmj)
     & + half*f1*Cv(N+cc001222,epmj)
     & + half*f2*Cv(N+cc002222,epmj)
     & - Cv(N+cc00222,epmj)*m1s)
      enddo
 61   continue
      enddo
      

C   Cv(N+cc0000001,Cv(N+cc0000002,zzzzzz)
      in(1,:)=
     & + f1*Cv(N+cc000000,:)
     & + Bv(bb000000+B13,:)
     & - Bv(bb000000+B23,:)
      in(2,:)=
     & + f2*Cv(N+cc000000,:)
     & + Bv(bb000000+B12,:)
     & - Bv(bb000000+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0000001,:)=in(1,:)
      Cv(N+cc0000002,:)=in(2,:)


C   Cv(N+cc0000111,Cv(N+cc0000112,zzzzpp)
      in(1,:)=
     & - Bv(bb0000+B23,:)
     & - 2*Bv(bb00001+B23,:)
     & + f1*Cv(N+cc000011,:)
     & - 4*Cv(N+cc0000001,:)
     & - Bv(bb000011+B23,:)
      in(2,:)=
     & - Bv(bb0000+B23,:)
     & - 2*Bv(bb00001+B23,:)
     & + f2*Cv(N+cc000011,:)
     & + Bv(bb000011+B12,:)
     & - Bv(bb000011+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0000111,:)=in(1,:)
      Cv(N+cc0000112,:)=in(2,:)

C   Cv(N+cc0000112,Cv(N+cc0000122,zzzzpk)
      in(1,:)=
     & + Bv(bb00001+B23,:)
     & + f1*Cv(N+cc000012,:)
     & - 2*Cv(N+cc0000002,:)
     & + Bv(bb000011+B23,:)
      in(2,:)=
     & + Bv(bb00001+B23,:)
     & + f2*Cv(N+cc000012,:)
     & - 2*Cv(N+cc0000001,:)
     & + Bv(bb000011+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0000112,:)=in(1,:)
      Cv(N+cc0000122,:)=in(2,:)

C   Cv(N+cc0000122,Cv(N+cc0000222,zzzzkk)
      in(1,:)=
     & + f1*Cv(N+cc000022,:)
     & + Bv(bb000011+B13,:)
     & - Bv(bb000011+B23,:)
      in(2,:)=
     & + f2*Cv(N+cc000022,:)
     & - 4*Cv(N+cc0000002,:)
     & - Bv(bb000011+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0000122,:)=in(1,:)
      Cv(N+cc0000222,:)=in(2,:)

C   Cv(N+cc0011111,Cv(N+cc0011112,zzpppp)
      in(1,:)=
     & - Bv(bb00+B23,:)
     & - 4*Bv(bb001+B23,:)
     & - 6*Bv(bb0011+B23,:)
     & - 4*Bv(bb00111+B23,:)
     & + f1*Cv(N+cc001111,:)
     & - 8*Cv(N+cc0000111,:)
     & - Bv(bb001111+B23,:)
      in(2,:)=
     & - Bv(bb00+B23,:)
     & - 4*Bv(bb001+B23,:)
     & - 6*Bv(bb0011+B23,:)
     & - 4*Bv(bb00111+B23,:)
     & + f2*Cv(N+cc001111,:)
     & + Bv(bb001111+B12,:)
     & - Bv(bb001111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0011111,:)=in(1,:)
      Cv(N+cc0011112,:)=in(2,:)

C   Cv(N+cc0011112,Cv(N+cc0011122,zzpppk)
      in(1,:)=
     & + Bv(bb001+B23,:)
     & + 3*Bv(bb0011+B23,:)
     & + 3*Bv(bb00111+B23,:)
     & + f1*Cv(N+cc001112,:)
     & - 6*Cv(N+cc0000112,:)
     & + Bv(bb001111+B23,:)
      in(2,:)=
     & + Bv(bb001+B23,:)
     & + 3*Bv(bb0011+B23,:)
     & + 3*Bv(bb00111+B23,:)
     & + f2*Cv(N+cc001112,:)
     & - 2*Cv(N+cc0000111,:)
     & + Bv(bb001111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0011112,:)=in(1,:)
      Cv(N+cc0011122,:)=in(2,:)

C   Cv(N+cc0011122,Cv(N+cc0011222,zzppkk)
      in(1,:)=
     & - Bv(bb0011+B23,:)
     & - 2*Bv(bb00111+B23,:)
     & + f1*Cv(N+cc001122,:)
     & - 4*Cv(N+cc0000122,:)
     & - Bv(bb001111+B23,:)
      in(2,:)=
     & - Bv(bb0011+B23,:)
     & - 2*Bv(bb00111+B23,:)
     & + f2*Cv(N+cc001122,:)
     & - 4*Cv(N+cc0000112,:)
     & - Bv(bb001111+B23,:)


      call pvBackSubst(G,2,perm, in)

      Cv(N+cc0011122,:)=in(1,:)
      Cv(N+cc0011222,:)=in(2,:)

C   Cv(N+cc0011222,Cv(N+cc0012222,zzpkkk)
      in(1,:)=
     & + Bv(bb00111+B23,:)
     & + f1*Cv(N+cc001222,:)
     & - 2*Cv(N+cc0000222,:)
     & + Bv(bb001111+B23,:)
      in(2,:)=
     & + Bv(bb00111+B23,:)
     & + f2*Cv(N+cc001222,:)
     & - 6*Cv(N+cc0000122,:)
     & + Bv(bb001111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0011222,:)=in(1,:)
      Cv(N+cc0012222,:)=in(2,:)

C   Cv(N+cc0012222,Cv(N+cc0022222,zzkkkk)
      in(1,:)=
     & + f1*Cv(N+cc002222,:)
     & + Bv(bb001111+B13,:)
     & - Bv(bb001111+B23,:)
      in(2,:)=
     & + f2*Cv(N+cc002222,:)
     & - 8*Cv(N+cc0000222,:)
     & - Bv(bb001111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc0012222,:)=in(1,:)
      Cv(N+cc0022222,:)=in(2,:)


C   Cv(N+cc1111111,Cv(N+cc1111112,pppppp)
      in(1,:)=
     & + f1*Cv(N+cc111111,:)
     & - 12*Cv(N+cc0011111,:)
     & - b11111sum(:)
     & - 5*b1111sum(:)
     & - 10*b111sum(:)
     & - 10*b11sum(:)
     & - 5*b1sum(:)
     & - b0sum(:)
      in(2,:)=
     & + f2*Cv(N+cc111111,:)
     & - b11111sum(:)
     & - 5*b1111sum(:)
     & - 10*b111sum(:)
     & - 10*b11sum(:)
     & - 5*b1sum(:)
     & - b0sum(:)
     & + Bv(bb111111+B12,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1111111,:)=in(1,:)
      Cv(N+cc1111112,:)=in(2,:)

      
C   Cv(N+cc1111112,Cv(N+cc1111122,pppppk)
      in(1,:)=
     & + b1sum(:)
     & + 4*b11sum(:)
     & + 6*b111sum(:)
     & + 4*b1111sum(:)
     & + b11111sum(:)
     & + f1*Cv(N+cc111112,:)
     & - 10*Cv(N+cc0011112,:)
      in(2,:)=
     & + b1sum(:)
     & + 4*b11sum(:)
     & + 6*b111sum(:)
     & + 4*b1111sum(:)
     & + b11111sum(:)
     & + f2*Cv(N+cc111112,:)
     & - 2*Cv(N+cc0011111,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1111112,:)=in(1,:)
      Cv(N+cc1111122,:)=in(2,:)

C   Cv(N+cc1111122,Cv(N+cc1111222,ppppkk)
      in(1,:)=
     & - b11sum(:)
     & - 3*b111sum(:)
     & - 3*b1111sum(:)
     & - b11111sum(:)
     & + f1*Cv(N+cc111122,:)
     & - 8*Cv(N+cc0011122,:)
      in(2,:)=
     & - b11sum(:)
     & - 3*b111sum(:)
     & - 3*b1111sum(:)
     & - b11111sum(:)
     & + f2*Cv(N+cc111122,:)
     & - 4*Cv(N+cc0011112,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1111122,:)=in(1,:)
      Cv(N+cc1111222,:)=in(2,:)

C   Cv(N+cc1111222,Cv(N+cc1112222,pppkkk)
      in(1,:)=
     & + b111sum(:)
     & + 2*b1111sum(:)
     & + b11111sum(:)
     & + f1*Cv(N+cc111222,:)
     & - 6*Cv(N+cc0011222,:)
      in(2,:)=
     & + b111sum(:)
     & + 2*b1111sum(:)
     & + b11111sum(:)
     & + f2*Cv(N+cc111222,:)
     & - 6*Cv(N+cc0011122,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1111222,:)=in(1,:)
      Cv(N+cc1112222,:)=in(2,:)

C   Cv(N+cc1112222,Cv(N+cc1122222,ppkkkk)
      in(1,:)=
     & - b1111sum(:)
     & - b11111sum(:)
     & + f1*Cv(N+cc112222,:)
     & - 4*Cv(N+cc0012222,:)
      in(2,:)=
     & - b1111sum(:)
     & - b11111sum(:)
     & + f2*Cv(N+cc112222,:)
     & - 8*Cv(N+cc0011222,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1112222,:)=in(1,:)
      Cv(N+cc1122222,:)=in(2,:)


C   Cv(N+cc1122222,Cv(N+cc1222222,pkkkkk)
      in(1,:)=
     & + f1*Cv(N+cc122222,:)
     & - 2*Cv(N+cc0022222,:)
     & + b11111sum(:)
      in(2,:)=
     & + f2*Cv(N+cc122222,:)
     & - 10*Cv(N+cc0012222,:)
     & + b11111sum(:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1122222,:)=in(1,:)
      Cv(N+cc1222222,:)=in(2,:)

C   Cv(N+cc1222222,Cv(N+cc2222222,kkkkkk)
      in(1,:)=
     & + f1*Cv(N+cc222222,:)
     & + Bv(bb111111+B13,:)
     & - Bv(bb111111+B23,:)
      in(2,:)=
     & + f2*Cv(N+cc222222,:)
     & - 12*Cv(N+cc0022222,:)
     & - Bv(bb111111+B23,:)
      call pvBackSubst(G,2,perm, in)
      Cv(N+cc1222222,:)=in(1,:)
      Cv(N+cc2222222,:)=in(2,:)



c--- to check recursion identities      
c      call Cfill_alt(p1,p2,p1p2,m1s,m2s,m3s,N)
       
      return
      end
