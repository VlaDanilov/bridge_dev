!---- J. Campbell, April 2017
!---- Power corrections for q-qbar initiated singlet processes
!----
!---- Formulae taken from:
!----   I.~Moult, L.~Rothen, I.~W.~Stewart, F.~J.~Tackmann and H.~X.~Zhu,
!----   ``Subleading Power Corrections for N-Jettiness Subtractions,''
!----   arXiv:1612.00450 [hep-ph].
!
!---- Eqs. (42) and (43) (MCFM-8.0 definition of tau)
!----  which simplify to Eqs. (33),(34),(36) and (37) in "tauboost" case
!----
!---- Alternative results taken from:
!----  R.~Boughezal, X.~Liu and F.~Petriello,
!----   ``Power Corrections in the N-jettiness Subtraction Scheme,''
!----   [arXiv:1612.02911 [hep-ph]].
!
!---- Eqs. (2.41), (2.42), (3.8) and (3.9) (MCFM-8.0 definition of tau)
!----   after manipulation to expand logs and integrate by expandpc.frm (appended)

      subroutine tau0_powcorr_qa(order,x1,x2,Qs,beama,beamb,msqpow)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'taucut.f'
      include 'scale.f'
      include 'facscale.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'kpart.f'
      real(dp):: x1,x2,msqpow(-nf:nf,-nf:nf),beama(-5:5),beamb(-5:5)
!      real(dp):: delpdfsa(-5:5),delpdfsb(-5:5)
      real(dp):: squaba(-5:5),squabb(-5:5),Qs,expY,L,Ltauc,tauc
      real(dp):: dxpdfa(-5:5),dxpdfb(-5:5),err,dx,faca,facb,Lmu,Y
      real(dp):: logs,logsa1,logsa2,logsa3,logsb1,logsb2,logsb3,fac,incsl
      real(dp), parameter:: xfrac=0.05_dp
      integer j,k,ih1,ih2
      integer order
      common/density/ih1,ih2
! Set this flag to .true. to use pure leading-power results of Moult et al.
! Set this flag to .false. to use expressions of Boughezal et al. that contain subleading logs
      logical, parameter:: useMoultetal=.true.
   
      msqpow(:,:)=zip
      
      dx=min(xfrac*x1,one-x1)      ! keep variation of x1 inside [0,1]
      call dxpdf_dfridr(ih1,1,x1,dx,err,dxpdfa)
      dx=min(xfrac*x2,one-x2)      ! keep variation of x2 inside [0,1]
      call dxpdf_dfridr(ih2,2,x2,dx,err,dxpdfb)

! apply additional factors of x1,x2, c.f. Eq. (35)
! and note that delta^\prime * f = - delta * f^\prime
      dxpdfa=-dxpdfa*x1
      dxpdfb=-dxpdfb*x2

      tauc=taucut/Qs
      L=log(tauc)
      Lmu=log(scale/Qs)

! can easily check that Y = log(x1/x2)/two = yraptwo(3,4,p)
      if (tauboost) then
        expY=one
      else
        expY=sqrt(x1/x2)
      endif
      Y=log(expY)

! Flag to turn off sub-leading terms generated when integrating over tau (default=0)
      incsl=zip
!      incsl=one

!---- fill PDFs
      do j=-nf,nf
      do k=-nf,nf
      
! do not fill qq and aa entries
      if (j*k > 0) cycle
      
      if (useMoultetal) then
!-------------------------------
! expressions from Moult et al |
!-------------------------------
        if( (order == 1) .or.
     &     ((order == 2) .and. (coeffonly .eqv. .false.)) ) then
          if((j .ne. 0) .and. (k .ne. 0)) then  
            msqpow(j,k)=ason4pi*four*CF
     &       *(beama(j)*(beamb(k)+dxpdfb(k))*expY*tauc*(L-one*incsl)
     &        +beamb(k)*(beama(j)+dxpdfa(j))/expY*tauc*(L-one*incsl))
          elseif((j .ne. 0) .and. (k == 0)) then 
            msqpow(j,k)=-ason4pi*two*TR*expY*beama(j)*beamb(k)*tauc*(L-one*incsl)
          elseif((j == 0) .and. (k .ne. 0)) then
            msqpow(j,k)=-ason4pi*two*TR/expY*beama(j)*beamb(k)*tauc*(L-one*incsl)
          endif
        endif

        if (order >= 2) then
          if((j .ne. 0) .and. (k .ne. 0)) then  
            msqpow(j,k)=msqpow(j,k)-ason4pi**2*sixteen*CF**2
     &       *(beama(j)*(beamb(k)+dxpdfb(k))*expY
     &         *tauc*(L**3-three*(L**2-two*L+two)*incsl)
     &        +beamb(k)*(beama(j)+dxpdfa(j))/expY
     &         *tauc*(L**3-three*(L**2-two*L+two)*incsl))
          elseif((j .ne. 0) .and. (k == 0)) then 
            msqpow(j,k)=msqpow(j,k)+ason4pi**2*four*TR*(CF+CA)*expY*beama(j)*beamb(k)
     &                 *tauc*(L**3-three*(L**2-two*L+two)*incsl)
          elseif((j == 0) .and. (k .ne. 0)) then
            msqpow(j,k)=msqpow(j,k)+ason4pi**2*four*TR*(CF+CA)/expY*beama(j)*beamb(k)
     &                 *tauc*(L**3-three*(L**2-two*L+two)*incsl)
          endif
        endif
      
      else
!-----------------------------------
! expressions from Boughezal et al |
!-----------------------------------
        if( (order == 1) .or.
     &     ((order == 2) .and. (coeffonly .eqv. .false.)) ) then
          if((j .ne. 0) .and. (k .ne. 0)) then  
            msqpow(j,k)=ason4pi*four*CF
     &       *(beama(j)*(beamb(k)+dxpdfb(k))*expY*tauc*(L-(one+Y)*incsl)
     &        +beamb(k)*(beama(j)+dxpdfa(j))/expY*tauc*(L-(one-Y)*incsl))
          elseif((j .ne. 0) .and. (k == 0)) then 
            msqpow(j,k)=-ason4pi*two*TR*expY*beama(j)*beamb(k)*tauc*(L-(one+Y)*incsl)
          elseif((j == 0) .and. (k .ne. 0)) then
            msqpow(j,k)=-ason4pi*two*TR/expY*beama(j)*beamb(k)*tauc*(L-(one-Y)*incsl)
          endif
        endif

        if (order >= 2) then
          if((j .ne. 0) .and. (k .ne. 0)) then  
            msqpow(j,k)=msqpow(j,k)-ason4pi**2*sixteen*CF**2
     &       *(beama(j)*(beamb(k)+dxpdfb(k))*expY
     &         *tauc*( - 6.D0*incsl + 2.D0*Lmu**2*incsl + 6.D0*L*incsl - 2.D0*L*Lmu**2 - 3.
     &    D0*L**2*incsl + L**3 - 2.D0*Y*incsl + 2.D0*Y*Lmu**2*incsl + 2.
     &    D0*Y*L*incsl - Y*L**2 + Y**2*incsl - Y**2*L + Y**3*incsl
     &      )
     &        +beamb(k)*(beama(j)+dxpdfa(j))/expY
     &         *tauc*( - 6.D0*incsl + 2.D0*Lmu**2*incsl + 6.D0*L*incsl - 2.D0
     &    *L*Lmu**2 - 3.D0*L**2*incsl + L**3 + 2.D0*Y*incsl - 2.D0*Y*
     &    Lmu**2*incsl - 2.D0*Y*L*incsl + Y*L**2 + Y**2*incsl - Y**2*L - Y**3*incsl
     &      ))
          elseif((j .ne. 0) .and. (k == 0)) then 
            msqpow(j,k)=msqpow(j,k)+ason4pi**2*four*TR*(CF+CA)*expY*beama(j)*beamb(k)
     &                 *tauc*( - 6.D0*incsl + 6.D0*L
     &    *incsl - 3.D0*L**2*incsl + L**3 - 2.D0*Y*incsl + 2.D0*Y*L*
     &    incsl - Y*L**2 + Y**2*incsl - Y**2*L + Y**3*incsl
     &      )
          elseif((j == 0) .and. (k .ne. 0)) then
            msqpow(j,k)=msqpow(j,k)+ason4pi**2*four*TR*(CF+CA)/expY*beama(j)*beamb(k)
     &                 *tauc*( - 6.D0*incsl + 6.D0*
     &    L*incsl - 3.D0*L**2*incsl + L**3 + 2.D0*Y*incsl - 2.D0*Y*L*
     &    incsl + Y*L**2 + Y**2*incsl - Y**2*L - Y**3*incsl
     &      )
          endif
        endif
      
      endif
      
      enddo
      enddo

      return
      end



















!--- s ason2pi,ason4pi,t0,m,mu,musq,Y,x1,x2,dx1,dx2,dt0,L,Lmu;
!--- s [pdf2+dxpdf2],[pdf1+dxpdf1],[-sixteen],two,four,CF,TR,[CF+CA],emY,eY,tauc,incsl;
!--- cf log;
!--- 
!--- g blpqq=1/2*(ason2pi*CF)^2
!---  *(8*log(t0/mu)^2-2*log(t0*m*eY/musq)^2-2*log(t0*m*emY/musq)^2)
!---  *(eY/m*log(t0/m*emY)*(2*x1*dx1)
!---   +emY/m*log(t0/m*eY)*(2*x2*dx2));
!--- 
!--- g blpqg=(ason2pi)^2*[CF+CA]*TR
!---  *(eY/m*log(m*eY/t0)^2*log(t0*eY/m));
!--- 
!--- g blpgq=(ason2pi)^2*[CF+CA]*TR
!---  *(emY/m*log(m*emY/t0)^2*log(t0*emY/m));
!--- 
!--- g blpqqnlo=ason2pi*CF
!---  *(eY/m*log(m*eY/t0)*(2*x1*dx1)
!---   +emY/m*log(m*emY/t0)*(2*x2*dx2));
!--- 
!--- g blpqgnlo=ason2pi*TR*(eY/m*log(m*eY/t0));
!--- 
!--- g blpgqnlo=ason2pi*TR*(emY/m*log(m*emY/t0));
!--- 
!--- id,log(m*eY/t0)=-log(t0/m*emY);
!--- id,log(m*emY/t0)=-log(t0/m*eY);
!--- 
!--- id,log(t0*m*eY/musq)=log(t0/m)+log(m^2*eY/musq);
!--- id,log(t0*m*emY/musq)=log(t0/m)+log(m^2*emY/musq);
!--- id,log(t0/mu)=log(t0/m)+1/2*log(m^2/musq);
!--- 
!--- id,log(m^2*eY/musq)=log(m^2/musq)+Y;
!--- id,log(m^2*emY/musq)=log(m^2/musq)-Y;
!--- 
!--- id,log(m^2/musq)=2*log(m/mu);
!--- *id,Y^2=log(eY)^2;
!--- 
!--- id,ason2pi=ason4pi*2;
!--- id,x1*dx1=-[pdf1+dxpdf1];
!--- id,x2*dx2=-[pdf2+dxpdf2];
!--- 
!--- InExpression,blpqq;
!--- multiply,-1/16*[-sixteen];
!--- EndInExpression;
!--- InExpression,blpqg,blpgq,blpqqnlo;
!--- multiply,1/4*four;
!--- EndInExpression;
!--- InExpression,blpqgnlo,blpgqnlo;
!--- multiply,1/2*two;
!--- EndInExpression;
!--- 
!--- multiply,dt0;
!--- 
!--- id,log(t0/m*emY)=log(t0/m)-Y;
!--- id,log(t0/m*eY)=log(t0/m)+Y;
!--- 
!--- *print +s;
!--- *.end
!--- 
!--- *** Integral log(b*t)^3*dt
!--- *       + t * (
!--- *          - 6
!--- *          + 6*log(t*b)
!--- *          - 3*log(t*b)^2
!--- *          + log(t*b)^3
!--- *          )
!--- 
!--- id,dt0*log(t0/m)^3=t0*(
!---           (- 6
!---           + 6*log(t0/m)
!---           - 3*log(t0/m)^2)*incsl
!---           + log(t0/m)^3);
!--- 
!--- *** Integral log(b*t)^2*dt
!--- *       + t * (
!--- *          + 2
!--- *          - 2*log(t*b)
!--- *          + log(t*b)^2
!--- *          )
!--- 
!--- id,dt0*log(t0/m)^2=t0*(
!---           (+ 2
!---           - 2*log(t0/m))*incsl
!---           + log(t0/m)^2);
!--- 
!--- *** Integral log(b*t)*dt
!--- *   bit2 =
!--- *       + t * (
!--- *          - 1
!--- *          + log(t*b)
!--- *          )
!--- 
!--- id,dt0*log(t0/m)=t0*(log(t0/m)-1*incsl);
!--- 
!--- *** Constants
!--- id,dt0=t0*incsl;
!--- 
!--- id,log(t0/m)=L;
!--- id,log(m/mu)=-Lmu;
!--- id,t0/m=tauc;
!--- 
!--- format doublefortran;
!--- b eY,emY,ason4pi,CF,x1,x2,dx1,dx2,m,[pdf1+dxpdf1],[pdf2+dxpdf2],[-sixteen],t0,dt0,
!---  [CF+CA],TR,two,four,tauc;
!--- print;
!--- .end

