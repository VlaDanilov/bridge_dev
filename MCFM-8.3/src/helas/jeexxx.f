C
C ----------------------------------------------------------------------
C
      SUBROUTINE JEEXXX(EB,EF,SHLF,CHLF,PHI,NHB,NHF,NSF , JEE)
C
C This subroutine computes an off-shell photon wavefunction emitted from
C the electron or positron beam, with a special care for the small angle
C region.  The momenta are measured in the laboratory frame, where the  
C e- (e+) beam is along the positive (negative) z axis.                 
C                                                                       
C INPUT:                                                                
C       real    EB             : energy (GeV)    of beam  e-/e+         
C       real    EF             : energy (GeV)    of final e-/e+         
C       real    SHLF           : sin(theta/2)    of final e-/e+         
C       real    CHLF           : cos(theta/2)    of final e-/e+         
C       real    PHI            : azimuthal angle of final e-/e+         
C       integer NHB  = -1 or 1 : helicity        of beam  e-/e+         
C       integer NHF  = -1 or 1 : helicity        of final e-/e+         
C       integer NSF  = -1 or 1 : +1 for electron, -1 for positron       
C                                                                       
C OUTPUT:                                                               
C       complex JEE(6)         : off-shell photon          J^mu(<e|A|e>)
C
      implicit none
      COMPLEX*16 JEE(6),COEFF
      REAL*8  CS(2),EB,EF,SHLF,CHLF,PHI,ME,ALPHA,GAL,HI,SF,SFH,X,ME2,Q2,
     &        RFP,RFM,SNP,CSP,RXC,C,S
      INTEGER NHB,NHF,NSF
C
      ME   =0.51099906D-3
      ALPHA=1./128.
      GAL  =SQRT(ALPHA*4.*3.14159265D0)
C
      HI =NHB
      SF =NSF
      SFH=NHB*NSF
      CS((3+NSF)/2)=SHLF
      CS((3-NSF)/2)=CHLF
C CS(1)=CHLF and CS(2)=SHLF for electron
C CS(1)=SHLF and CS(2)=CHLF for positron
      X=EF/EB
      ME2=ME**2
      Q2=-4.*CS(2)**2*(EF*EB-ME2)
     &   +SF*(1.-X)**2/X*(SHLF+CHLF)*(SHLF-CHLF)*ME2
      RFP=(1+NSF)
      RFM=(1-NSF)
      SNP=SIN(PHI)
      CSP=COS(PHI)
C
      IF (NHB.EQ.NHF) THEN
         RXC=2.*X/(1.-X)*CS(1)**2
         COEFF= GAL*2.*EB*SQRT(X)*CS(2)/Q2
     &         *(dCMPLX( RFP )-RFM*dCMPLX( CSP ,-SNP*HI ))*.5
         JEE(1) =  dCMPLX( 0.d0 )
         JEE(2) =  COEFF*dCMPLX( (1.+RXC)*CSP ,-SFH*SNP )
         JEE(3) =  COEFF*dCMPLX( (1.+RXC)*SNP , SFH*CSP )
         JEE(4) =  COEFF*(-SF*RXC/CS(1)*CS(2))
      ELSE
         COEFF= GAL*ME/Q2/SQRT(X)
     &         *(dCMPLX( RFP )+RFM*dCMPLX( CSP , SNP*HI ))*.5*HI
         JEE(1) = -COEFF*(1.+X)*CS(2)*dCMPLX( CSP , SFH*SNP )
         JEE(2) =  COEFF*(1.-X)*CS(1)
         JEE(3) =  JEE(2)*dCMPLX( 0.d0 , SFH )
         JEE(4) =  JEE(1)*SF*(1.-X)/(1.+X)
      ENDIF
C
      C=(CHLF+SHLF)*(CHLF-SHLF)
      S=2.*CHLF*SHLF
C
      JEE(5) = -EB*dCMPLX( 1.-X , SF-X*C )
      JEE(6) =  EB*X*S*dCMPLX( CSP , SNP )
C
      RETURN          
      END
C
