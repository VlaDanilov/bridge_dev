      subroutine pvBfill(p1sq,m0sq,m1sq,N)
!    N is the offset in the storage
!    Implementation of the formula of Denner and Dittmaier
!@article{Denner:2005nn,
!      author         = "Denner, Ansgar and Dittmaier, S.",
!      title          = "{Reduction schemes for one-loop tensor integrals}",
!      journal        = "Nucl. Phys.",
!      volume         = "B734",
!      year           = "2006",
!      pages          = "62-115",
!      doi            = "10.1016/j.nuclphysb.2005.11.007",
!      eprint         = "hep-ph/0509141",
!      archivePrefix  = "arXiv",
!      primaryClass   = "hep-ph",
!      reportNumber   = "MPP-2005-84, PSI-PR-05-08",
!      SLACcitation   = "%%CITATION = HEP-PH/0509141;%%"
! }
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'pvAnames.f'
      include 'pvBnames.f'
      include 'TRconstants.f'
      include 'TRonshellcutoff.f'
      include 'pvAv.f'
      include 'pvBv.f'
      include 'TRscale.f'
      include 'pvverbose.f'
      integer:: N,Np,ep,A2,pvAcache,nl,in,iP
      real(dp)::p1sq,m0sq,m1sq,f1,iep
      complex(dp):: xp,xm,rt,arg,arg1,pvfndd,xpvfndd
      real(dp)::fac,facnp
      real(dp),save:: idp1(0:2),id(0:2),idm1(0:2)
      logical:: p1sqnonzero
      logical,save:: first=.true.
      logical,save:: scaleset=.false.
!$omp threadprivate(first,scaleset,idp1,id,idm1)
C-----statement functions
      include 'cplx.h'
      fac(n)=(-one)**n/dfloat(n+1)
      facnp(in,iP)=(-one)**(iP-2*in-1)
      
      if (first) then
      first=.false.
C--idp1=1/[D+1]
      idp1(0)=0.2_dp
      idp1(1)=idp1(0)*0.4_dp
      idp1(2)=idp1(1)*0.4_dp
C--id=1/D
      id(0)=0.25_dp
      id(1)=id(0)*half
      id(2)=id(1)*half
C--idm1=1/[D-1]
      idm1(0)=one/three
      idm1(1)=idm1(0)*two/three
      idm1(2)=idm1(1)*two/three
      endif

      if (scaleset .neqv. .true.) then
      scaleset=.true.
      if ((scale .eq. -1d12) .and. (musq .eq. -1d12)) then
      write(6,*) 'Did you forget to call setmudim?'
      write(6,*) 'Setting scale to scale=one'
      scale=one
      musq=one
      endif
      endif

      A2=pvAcache(m1sq)

C----self energies never contain double poles -- set to zero
      do Np=N+1,N+Nbb
      Bv(Np,-2)=czip
      enddo 

C----calculate B0
c      do ep=-1,0
c      Bv(bb0+N,ep)=trI2(p1sq,m0sq,m1sq,musq,ep)
c      enddo

      if ((abs(p1sq/musq) .lt. onshellcutoff) 
     . .and. (abs(m0sq/musq) .lt. onshellcutoff) 
     . .and. (abs(m1sq/musq) .lt. onshellcutoff)) then
      
      if (pvverbose) then
        write(6,*) 'setting zero psq, zero mass self-energy to zero'
        write(6,*) 'p1sq=',p1sq,m0sq,m1sq
      endif
      
      do Np=N+1,N+Nbb
      do ep=-1,0
      Bv(Np,ep)=czip
      enddo 
      enddo 
      return

      elseif (abs(m0sq/musq) .lt. onshellcutoff) then

C---deal with special cases for m0sq=0, p1sq=m1sq, DD(4.13)
      if (abs((p1sq-m1sq)/musq) .lt. onshellcutoff) then
      nl=0
      Bv(bb0+N,-1)=cplx1(fac(nl))
      Bv(bb0+N,0)=Bv(bb0+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      nl=1
      Bv(bb1+N,-1)=cplx1(fac(nl))
      Bv(bb1+N,0)=Bv(bb1+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      nl=2
      Bv(bb11+N,-1)=cplx1(fac(nl))
      Bv(bb11+N,0)=Bv(bb11+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      nl=3
      Bv(bb111+N,-1)=cplx1(fac(nl))
      Bv(bb111+N,0)=Bv(bb111+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      nl=4
      Bv(bb1111+N,-1)=cplx1(fac(nl))
      Bv(bb1111+N,0)=Bv(bb1111+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      nl=5
      Bv(bb11111+N,-1)=cplx1(fac(nl))
      Bv(bb11111+N,0)=Bv(bb11111+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      nl=6
      Bv(bb111111+N,-1)=cplx1(fac(nl))
      Bv(bb111111+N,0)=Bv(bb111111+N,-1)
     . *cplx1(log(musq/m1sq)+two/dfloat(nl+1))

      else
C---deal with special cases for m0sq=0, DD(4.12)
      arg=cplx1(one-m1sq/p1sq)
      arg1=cplx1(m1sq-p1sq)
      iep=sign(one,p1sq)

c--- if p1sq=0 too, root is formally infinite and pvfndd should be zero
      if (abs(p1sq/musq) .gt. onshellcutoff) then
        p1sqnonzero=.true.
      else
        p1sqnonzero=.false.
        xpvfndd=czip
      endif

      nl=0
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb0+N,-1)=cplx1(fac(nl))
      Bv(bb0+N,0)=Bv(bb0+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)
      nl=1
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb1+N,-1)=cplx1(fac(nl))
      Bv(bb1+N,0)=Bv(bb1+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)
      nl=2
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb11+N,-1)=cplx1(fac(nl))
      Bv(bb11+N,0)=Bv(bb11+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)
      nl=3
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb111+N,-1)=cplx1(fac(nl))
      Bv(bb111+N,0)=Bv(bb111+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)
      nl=4
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb1111+N,-1)=cplx1(fac(nl))
      Bv(bb1111+N,0)=Bv(bb1111+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)
      nl=5
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb11111+N,-1)=cplx1(fac(nl))
      Bv(bb11111+N,0)=Bv(bb11111+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)
      nl=6
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      Bv(bb111111+N,-1)=cplx1(fac(nl))
      Bv(bb111111+N,0)=Bv(bb111111+N,-1)
     . *(cplx1(log(musq)+one/dfloat(nl+1))-cln(arg1,-one)-xpvfndd)

      endif

      elseif (abs(p1sq/musq) .lt. onshellcutoff) then
C---deal with special case, p1sq=0
      xp=cplx1(m0sq/(m0sq-m1sq))  ! other root is formally infinite

      nl=0
      Bv(bb0+N,-1)=cplx1(fac(nl))
      Bv(bb0+N,0)=Bv(bb0+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))
      nl=1
      Bv(bb1+N,-1)=cplx1(fac(nl))
      Bv(bb1+N,0)=Bv(bb1+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))
      nl=2
      Bv(bb11+N,-1)=cplx1(fac(nl))
      Bv(bb11+N,0)=Bv(bb11+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))
      nl=3
      Bv(bb111+N,-1)=cplx1(fac(nl))
      Bv(bb111+N,0)=Bv(bb111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))
      nl=4
      Bv(bb1111+N,-1)=cplx1(fac(nl))
      Bv(bb1111+N,0)=Bv(bb1111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))
      nl=5
      Bv(bb11111+N,-1)=cplx1(fac(nl))
      Bv(bb11111+N,0)=Bv(bb11111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))
      nl=6
      Bv(bb111111+N,-1)=cplx1(fac(nl))
      Bv(bb111111+N,0)=Bv(bb111111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one))

      else
C----general case, DD (4.8)
      rt=sqrt(cplx1((m1sq-m0sq-p1sq)**2-four*p1sq*m0sq))
      xp=half*(-cplx1(m1sq-m0sq-p1sq)+rt)/p1sq
      xm=half*(-cplx1(m1sq-m0sq-p1sq)-rt)/p1sq

      nl=0
      Bv(bb0+N,-1)=cplx1(fac(nl))
      Bv(bb0+N,0)=Bv(bb0+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))
      nl=1
      Bv(bb1+N,-1)=cplx1(fac(nl))
      Bv(bb1+N,0)=Bv(bb1+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))
      nl=2
      Bv(bb11+N,-1)=cplx1(fac(nl))
      Bv(bb11+N,0)=Bv(bb11+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))
      nl=3
      Bv(bb111+N,-1)=cplx1(fac(nl))
      Bv(bb111+N,0)=Bv(bb111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))
      nl=4
      Bv(bb1111+N,-1)=cplx1(fac(nl))
      Bv(bb1111+N,0)=Bv(bb1111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))
      nl=5
      Bv(bb11111+N,-1)=cplx1(fac(nl))
      Bv(bb11111+N,0)=Bv(bb11111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))
      nl=6
      Bv(bb111111+N,-1)=cplx1(fac(nl))
      Bv(bb111111+N,0)=Bv(bb111111+N,-1)
     . *(cplx1(log(musq/m0sq))-pvfndd(nl,xp,one)-pvfndd(nl,xm,-one))

      endif

c--- Construct tensors involving the metric (these are general)
c---  following the formulae of DD (4.4)
      f1=m1sq-m0sq-p1sq
      
c--- one factor of the metric
      in=0
      
      iP=2
      do ep=-1,0
      Bv(bb00+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa0+A2,ep)-f1*Bv(bb1+N,ep)
     . +two*p1sq*Bv(bb11+N,ep))
      enddo
      
      iP=3
      do ep=-1,0
      Bv(bb001+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa0+A2,ep)-f1*Bv(bb11+N,ep)
     . +two*p1sq*Bv(bb111+N,ep))
      enddo
      
      iP=4
      do ep=-1,0
      Bv(bb0011+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa0+A2,ep)-f1*Bv(bb111+N,ep)
     . +two*p1sq*Bv(bb1111+N,ep))
      enddo
      
      iP=5
      do ep=-1,0
      Bv(bb00111+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa0+A2,ep)-f1*Bv(bb1111+N,ep)
     . +two*p1sq*Bv(bb11111+N,ep))
      enddo
      
      iP=6
      do ep=-1,0
      Bv(bb001111+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa0+A2,ep)-f1*Bv(bb11111+N,ep)
     . +two*p1sq*Bv(bb111111+N,ep))
      enddo
      
c--- two factors of the metric
      in=1
      
      iP=4
      do ep=-1,0
      Bv(bb0000+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa00+A2,ep)-f1*Bv(bb001+N,ep)
     . +two*p1sq*Bv(bb0011+N,ep))
      enddo
      
      iP=5
      do ep=-1,0
      Bv(bb00001+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa00+A2,ep)-f1*Bv(bb0011+N,ep)
     . +two*p1sq*Bv(bb00111+N,ep))
      enddo
      
      iP=6
      do ep=-1,0
      Bv(bb000011+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa00+A2,ep)-f1*Bv(bb00111+N,ep)
     . +two*p1sq*Bv(bb001111+N,ep))
      enddo
      
c--- three factors of the metric

      in=2
      iP=6
      do ep=-1,0
      Bv(bb000000+N,ep)=-half/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*Av(aa0000+A2,ep)-f1*Bv(bb00001+N,ep)
     . +two*p1sq*Bv(bb000011+N,ep))
      enddo
      
      return 
      end

