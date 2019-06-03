!---- CW Dec 16
!---- Power corrections for gg initiated singlet processes
!-----formula taken from
!-----1612.02911 eq 40 (NLO order =1) and
!-----eq 50
      subroutine tau0_powcorr_gg(order,x1,x2,Qs,beama,beamb,msqpow)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'taucut.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'kpart.f'
      real(dp):: x1,x2,msqpow(-nf:nf,-nf:nf),beama(-5:5),beamb(-5:5)
!      real(dp):: delpdfsa(-5:5),delpdfsb(-5:5)
      real(dp):: squaba(-5:5),squabb(-5:5),Qs,expY,La,Lb,La2,Lb2,Ltau
      real(dp):: dxpdfa(-5:5),dxpdfb(-5:5),err,dx,faca,facb
      real(dp):: logs,logsa1,logsa2,logsa3,logsb1,logsb2,logsb3
      real(dp), parameter:: xfrac=0.05_dp
      integer j,k,ih1,ih2
      integer order
      common/density/ih1,ih2

      msqpow(:,:)=zip

!      call midpoint_devpdf(ih1,1,x1,delpdfsa)
!      call midpoint_devpdf(ih2,2,x2,delpdfsb)
!      dxpdfa=delpdfsa
!      dxpdfb=delpdfsb

      dx=min(xfrac*x1,one-x1)      ! keep variation of x1 inside [0,1]
      call dxpdf_dfridr(ih1,1,x1,dx,err,dxpdfa)
      dx=min(xfrac*x2,one-x2)      ! keep variation of x2 inside [0,1]
      call dxpdf_dfridr(ih2,2,x2,dx,err,dxpdfb)
      
!      do j=0,0
!      write(6,*) x1,j,delpdfsa(j),dxpdfa(j)
!      enddo
!      do j=0,0
!      write(6,*) x2,j,delpdfsb(j),dxpdfb(j)
!      enddo
!      pause

! can easily check that Y = log(x1/x2)/two = yraptwo(3,4,p)
      if (tauboost) then
        expY=one
      else
        expY=sqrt(x1/x2)
      endif
      faca=expY/Qs
      La=log(Qs*expY/taucut)
      La2=log(taucut*Qs*expY/scale**2)
      facb=one/expY/Qs
      Lb=log(Qs/expY/taucut)
      Lb2=log(taucut*Qs/expY/scale**2)
      Ltau=log(taucut/scale)

!---- fill PDFs
      do j=-nf,nf
      do k=-nf,nf

        squaba(j) = two*(x1*dxpdfa(j)+beama(j))
        squabb(k) = two*(x2*dxpdfb(k)+beamb(k))

        if( (order == 1) .or.
     &     ((order == 2) .and. (coeffonly .eqv. .false.)) ) then
          if((j.eq.0).and.(k.eq.0)) then  
            msqpow(j,k)=ason2pi*CA*taucut*
     &        (faca*(zip*one+La)*squaba(j)*beamb(k)
     &        +facb*(zip*one+Lb)*beama(j)*squabb(k))
          elseif((j==0).and.(k.ne.0)) then
            msqpow(j,k)=ason2pi*CF*taucut*faca*(zip*one+La)*beama(j)*beamb(k)
          elseif((k==0).and.(j.ne.0)) then 
            msqpow(j,k)=ason2pi*CF*taucut*facb*(zip*one+Lb)*beama(j)*beamb(k)
          endif
        endif

        if (order >= 2) then
          if((j==0) .and. (k==0)) then  
            logsa1=La*Ltau**2+(Ltau**2-two*La*Ltau-four*Ltau+two*La+six)!*zip
            logsa2=La*La2**2+(La2**2-two*La*La2-four*La2+two*La+six)!*zip
            logsa3=La*Lb2**2+(Lb2**2-two*La*Lb2-four*Lb2+two*La+six)!*zip
            logsb1=Lb*Ltau**2+(Ltau**2-two*Lb*Ltau-four*Ltau+two*Lb+six)!*zip
            logsb2=Lb*La2**2+(La2**2-two*Lb*La2-four*La2+two*Lb+six)!*zip
            logsb3=Lb*Lb2**2+(Lb2**2-two*Lb*Lb2-four*Lb2+two*Lb+six)!*zip
            msqpow(j,k)=msqpow(j,k)+half*(ason2pi*CA)**2*taucut*(
     &            faca*squaba(j)*beamb(k)*(-eight*logsa1+two*logsa2+two*logsa3)
     &           +facb*beama(j)*squabb(k)*(-eight*logsb1+two*logsb2+two*logsb3))
          elseif((j==0) .and. (k .ne. 0)) then 
            logs=La**2*log(taucut*expY/Qs)
     &          +(two*La*log(taucut*expY/Qs)
     &           -La**2+two*log(taucut*expY/Qs)-four*La-six)!*zip
            msqpow(j,k)=msqpow(j,k)+ason2pi**2*(CF+CA)*CF*taucut
     &           *faca*logs*beama(j)*beamb(k)
          elseif((k == 0) .and. (j .ne. 0)) then
            logs=Lb**2*log(taucut/expY/Qs)
     &          +(two*Lb*log(taucut/expY/Qs)
     &           -Lb**2+two*log(taucut/expY/Qs)-four*Lb-six)!*zip
            msqpow(j,k)=msqpow(j,k)+ason2pi**2*(CF+CA)*CF*taucut
     &           *facb*logs*beama(j)*beamb(k)
          endif
        endif
         
      enddo
      enddo

      return
      end

