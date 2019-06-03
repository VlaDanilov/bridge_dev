
      double precision sqrts,tauh,rmur,rmuf,lfh,lft,lfr,lrt,lth,elwfac
      double precision apimuR,mbpole,ataut2,nf
      double precision mbMSbarmuR,mbMSbarmuRbbh,mbMSbarmt,mbMSbarmuB,
     &mbMSbarmuD
      double precision c1sm0,c1sm1,c1sm2,c2sm1
      double precision c1eff0,c1eff1,c1eff2,c2eff1
      double precision gth,gth11,gth12,gth21,gth22
      double precision gbh,gbh11,gbh12,gbh21,gbh22
      
      double precision rho0,apimz,myapimz,gfermi,mt,mbmb
      double precision gt,gb,dgb

      double precision alphamuR,muRfactor,muBfactor,AlfasMZ,betanull
      double precision sigmanull,cme,Mbin,Msb1in,Msb2in,Mtin
      double precision Mst1in,Mst2in,Mw,Mh,muR,muB,muD,Mz

      integer norderggh,norderbbh,renscheme,tanbresum
      integer ew

      common/params/ sqrts,tauh,ataut2,rmur,rmuf,lfh,lft,lfr,lrt,lth,
     &elwfac,apimuR,mbMSbarmuR,mbMSbarmuRbbh,mbMSbarmt,mbMSbarmuB,
     &mbMSbarmuD,mbpole,nf,c1sm0,c1sm1,c1sm2,c2sm1,c1eff0,
     &c1eff1,c1eff2,c2eff1,gth,gth11,gth12,gth21,gth22,gbh,gbh11,
     &gbh12,gbh21,gbh22

      common/vars/rho0,apimz,myapimz,gfermi,mt,mbmb,gt,gb,dgb

      common /inp/ alphamuR,muRfactor,muBfactor,AlfasMZ,cme,
     &norderggh,norderbbh,renscheme,ew,tanbresum

      common /constants/ betanull,sigmanull,Mw,Mz,Mh,muR,muB,muD
