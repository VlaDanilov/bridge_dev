      ! See the module zajj_treeamps_m for the amplitudes itself.
      
      module zajj_anomcoup_m
        use zajj_treeamps_m
        implicit none

      contains

      subroutine amp_zqqQQa_anomZZ(j1,j2,j3,j4,j5,j6,j7,za,zb,ha,hb,ai)
        use zajj_treeamps_m
        implicit none
        
        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: ha,hb
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: ai(2,2,2,2)

        ! imo it's much cleaner to just enumerate all 16 possibilities
        ! than trying to write loops with several convoluted conditions
        ! and fixing/adjusting signs at various places

        ! hq,hQ,ha,hl

        ai(2,2,2,2) = zajj_tree_qqQQ_anomzz_plus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb)
        ai(2,2,1,2) = zajj_tree_qqQQ_anomzz_minus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb)

        ! lepton exchange
        ai(2,2,2,1) = zajj_tree_qqQQ_anomzz_plus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb)
        ai(2,2,1,1) = zajj_tree_qqQQ_anomzz_minus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb)

        ! Q exchange
        ai(2,1,2,2) = zajj_tree_qqQQ_anomzz_plus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb)
        ai(2,1,1,2) = zajj_tree_qqQQ_anomzz_minus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb)

        ! lepton +  Q exchange
        ai(2,1,2,1) = zajj_tree_qqQQ_anomzz_plus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb)
        ai(2,1,1,1) = zajj_tree_qqQQ_anomzz_minus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb)

        ai(1,1,1,1) = conjg(ai(2,2,2,2))
        ai(1,1,2,1) = conjg(ai(2,2,1,2))
        
        ai(1,1,1,2) = conjg(ai(2,2,2,1))
        ai(1,1,2,2) = conjg(ai(2,2,1,1))

        ai(1,2,1,1) = conjg(ai(2,1,2,2))
        ai(1,2,2,1) = conjg(ai(2,1,1,2))

        ai(1,2,1,2) = conjg(ai(2,1,2,1))
        ai(1,2,2,2) = conjg(ai(2,1,1,1))

c       ai(1,1,1,1) = conjg(zajj_tree_qqQQ_anomzz_plus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb))
c       ai(1,1,2,1) = conjg(zajj_tree_qqQQ_anomzz_minus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb))

c       ai(1,1,1,2) = conjg(zajj_tree_qqQQ_anomzz_plus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb))
c       ai(1,1,2,2) = conjg(zajj_tree_qqQQ_anomzz_minus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb))

c       ai(1,2,1,1) = conjg(zajj_tree_qqQQ_anomzz_plus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb))
c       ai(1,2,2,1) = conjg(zajj_tree_qqQQ_anomzz_minus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb))

c       ai(1,2,1,2) = conjg(zajj_tree_qqQQ_anomzz_plus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb))
c       ai(1,2,2,2) = conjg(zajj_tree_qqQQ_anomzz_minus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb))

      end

      subroutine amp_zqqQQa_anomZa(j1,j2,j3,j4,j5,j6,j7,za,zb,ha,hb,ai)
        use zajj_treeamps_m
        implicit none
        
        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: ha,hb
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: ai(2,2,2,2)

        ! hq,hQ,ha,hl

        ai(2,2,2,2) = zajj_tree_qqQQ_anomZa_plus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb)
        ai(2,2,1,2) = zajj_tree_qqQQ_anomZa_minus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb)

        ! lepton exchange
        ai(2,2,2,1) = zajj_tree_qqQQ_anomZa_plus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb)
        ai(2,2,1,1) = zajj_tree_qqQQ_anomZa_minus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb)

        ! Q exchange
        ai(2,1,2,2) = zajj_tree_qqQQ_anomZa_plus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb)
        ai(2,1,1,2) = zajj_tree_qqQQ_anomZa_minus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb)

        ! lepton +  Q exchange
        ai(2,1,2,1) = zajj_tree_qqQQ_anomZa_plus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb)
        ai(2,1,1,1) = zajj_tree_qqQQ_anomZa_minus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb)

        ai(1,1,1,1) = conjg(ai(2,2,2,2))
        ai(1,1,2,1) = conjg(ai(2,2,1,2))
        
        ai(1,1,1,2) = conjg(ai(2,2,2,1))
        ai(1,1,2,2) = conjg(ai(2,2,1,1))

        ai(1,2,1,1) = conjg(ai(2,1,2,2))
        ai(1,2,2,1) = conjg(ai(2,1,1,2))

        ai(1,2,1,2) = conjg(ai(2,1,2,1))
        ai(1,2,2,2) = conjg(ai(2,1,1,1))

c       ai(1,1,1,1) = conjg(zajj_tree_qqQQ_anomZa_plus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb))
c       ai(1,1,2,1) = conjg(zajj_tree_qqQQ_anomZa_minus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb))

c       ai(1,1,1,2) = conjg(zajj_tree_qqQQ_anomZa_plus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb))
c       ai(1,1,2,2) = conjg(zajj_tree_qqQQ_anomZa_minus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb))

c       ai(1,2,1,1) = conjg(zajj_tree_qqQQ_anomZa_plus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb))
c       ai(1,2,2,1) = conjg(zajj_tree_qqQQ_anomZa_minus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb))

c       ai(1,2,1,2) = conjg(zajj_tree_qqQQ_anomZa_plus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb))
c       ai(1,2,2,2) = conjg(zajj_tree_qqQQ_anomZa_minus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb))

      end

      subroutine amp_zqqQQa_anomaZ(j1,j2,j3,j4,j5,j6,j7,za,zb,ha,hb,ai)
        use zajj_treeamps_m
        implicit none

        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb
        complex(dp), intent(out) :: ai(2,2,2,2)

        ! hq,hQ,ha,hl

        ai(2,2,2,2) = zajj_tree_qqQQ_anomaZ_plus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb)
        ai(2,2,1,2) = zajj_tree_qqQQ_anomaZ_minus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb)

        ! lepton exchange
        ai(2,2,2,1) = zajj_tree_qqQQ_anomaZ_plus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb)
        ai(2,2,1,1) = zajj_tree_qqQQ_anomaZ_minus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb)

        ! Q exchange
        ai(2,1,2,2) = zajj_tree_qqQQ_anomaZ_plus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb)
        ai(2,1,1,2) = zajj_tree_qqQQ_anomaZ_minus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb)

        ! lepton +  Q exchange
        ai(2,1,2,1) = zajj_tree_qqQQ_anomaZ_plus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb)
        ai(2,1,1,1) = zajj_tree_qqQQ_anomaZ_minus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb)

        ai(1,1,1,1) = conjg(ai(2,2,2,2))
        ai(1,1,2,1) = conjg(ai(2,2,1,2))
        
        ai(1,1,1,2) = conjg(ai(2,2,2,1))
        ai(1,1,2,2) = conjg(ai(2,2,1,1))

        ai(1,2,1,1) = conjg(ai(2,1,2,2))
        ai(1,2,2,1) = conjg(ai(2,1,1,2))

        ai(1,2,1,2) = conjg(ai(2,1,2,1))
        ai(1,2,2,2) = conjg(ai(2,1,1,1))


c       ai(1,1,1,1) = conjg(zajj_tree_qqQQ_anomaZ_plus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb))
c       ai(1,1,2,1) = conjg(zajj_tree_qqQQ_anomaZ_minus(j2,j1,j6,j7,j5,j3,j4,za,zb,ha,hb))

c       ai(1,1,1,2) = conjg(zajj_tree_qqQQ_anomaZ_plus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb))
c       ai(1,1,2,2) = conjg(zajj_tree_qqQQ_anomaZ_minus(j2,j1,j7,j6,j5,j3,j4,za,zb,ha,hb))

c       ai(1,2,1,1) = conjg(zajj_tree_qqQQ_anomaZ_plus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb))
c       ai(1,2,2,1) = conjg(zajj_tree_qqQQ_anomaZ_minus(j2,j1,j6,j7,j5,j4,j3,za,zb,ha,hb))

c       ai(1,2,1,2) = conjg(zajj_tree_qqQQ_anomaZ_plus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb))
c       ai(1,2,2,2) = conjg(zajj_tree_qqQQ_anomaZ_minus(j2,j1,j7,j6,j5,j4,j3,za,zb,ha,hb))

      end

      ! this is copy and pasted from amp_qqagg_ql, with added arguments for
      ! anomalous couplings
      function amp_qqagg_anomZZ(i1,qh,i2,h2,i3,h3,i4,h4,
     & i5,lh,j6,j7,za,zb,ha,hb)
        use zajj_treeamps_m
        implicit none
        !include 'types.f'
        complex(dp):: amp_qqagg_anomZZ
        
        include 'constants.f'
        include 'nf.f'
        !include 'mxpart.f'
        include 'cplx.h'

        integer, intent(in) :: i1,i2,i3,i4,i5,j6,j7
        integer, intent(in) :: qh,h2,h3,h4,lh
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        integer:: i6,i7,j,k,g1,g2,g3,g4,lg
        complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)
        complex(dp):: A7h1,A7h2,A7h3,A7h4,A7h5,A7h6,A7h7,A7h8
c-----
C---A(qh,h2,h3,h4,lh)
C---h=1 LH
C---h=2 RH
        g1=qh
        g2=h2
        g3=h3
        g4=h4
        lg=lh
        i6=j6
        i7=j7
        if (g1*lg==2)then
          i7=j6
          i6=j7
          lg=3-lg
        endif

        if (g1==2) then
          xa(:,:) = za(:,:)
          xb(:,:) = zb(:,:)
        elseif (g1==1) then
          g1=2
          g2=3-g2
          g3=3-g3
          g4=3-g4
          lg=3-lg
          
          xa(:,:) = zb(:,:)
          xb(:,:) = za(:,:)
        endif

c-----A7h1
c-----A(2,2,2,2,2)
      if (
     & (g1==2).and.(g2==2).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,2,2,2)
      A7h1=czip
      A7h1 = zajj_tree_qqgg_anomZZ_ppp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h1
c-----A7h2
c-----A(2,2,2,1,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,2,1,2)
      A7h2=czip
      A7h2 = zajj_tree_qqgg_anomZZ_ppm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h2
c-----A7h3
c-----A(2,2,1,2,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,1,2,2)
      A7h3=czip
      A7h3 = zajj_tree_qqgg_anomZZ_pmp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h3
c-----A7h4
c-----A(2,1,2,2,2)
      elseif
     .((g1==2).and.(g2==1).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,2,2,2)
      A7h4=czip
      A7h4 = zajj_tree_qqgg_anomZZ_mpp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h4
c-----A7h5
c-----A(2,2,1,1,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,1,1,2)
      A7h5=czip
      A7h5 = zajj_tree_qqgg_anomZZ_pmm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h5
c-----A7h6
c-----A(2,1,2,1,2) =
      elseif
     .((g1==2).and.(g2==1).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,2,1,2)
      A7h6=czip
      A7h6 = zajj_tree_qqgg_anomZZ_mpm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h6
c-----A7h7
c-----A(2,1,1,2,2) =
      elseif
     .((g1==2).and.(g2==1).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,1,2,2)
      A7h7=czip
      A7h7 = zajj_tree_qqgg_anomZZ_mmp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h7
c-----A7h8
c-----A(2,1,1,1,2)
      elseif
     .((g1==2).and.(g2==1).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,1,1,2)
      A7h8=czip
      A7h8 = zajj_tree_qqgg_anomZZ_mmm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZZ = A7h8
c----
      else
      amp_qqagg_anomZZ = czip
      call abort
      endif

      end

      function amp_qqagg_anomZa(i1,qh,i2,h2,i3,h3,i4,h4,
     & i5,lh,j6,j7,za,zb,ha,hb)
        use zajj_treeamps_m
        implicit none
        !include 'types.f'
        complex(dp):: amp_qqagg_anomZa
        
        include 'constants.f'
        include 'nf.f'
        !include 'mxpart.f'
        include 'cplx.h'

        integer, intent(in) :: i1,i2,i3,i4,i5,j6,j7
        integer, intent(in) :: qh,h2,h3,h4,lh
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        integer:: i6,i7,j,k,g1,g2,g3,g4,lg
        complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)
        complex(dp):: A7h1,A7h2,A7h3,A7h4,A7h5,A7h6,A7h7,A7h8
        complex(dp) :: test
c-----
C---A(qh,h2,h3,h4,lh)
C---h=1 LH
C---h=2 RH
        g1=qh
        g2=h2
        g3=h3
        g4=h4
        lg=lh
        i6=j6
        i7=j7
        if (g1*lg==2)then
          i7=j6
          i6=j7
          lg=3-lg
        endif

        if (g1==2) then
          xa(:,:)=za(:,:)
          xb(:,:)=zb(:,:)
        elseif (g1==1) then
          g1=2
          g2=3-g2
          g3=3-g3
          g4=3-g4
          lg=3-lg

          xa(:,:)=zb(:,:)
          xb(:,:)=za(:,:)
        endif

c-----A7h1
c-----A(2,2,2,2,2)
      if (
     & (g1==2).and.(g2==2).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,2,2,2)
      A7h1=czip
      A7h1 = zajj_tree_qqgg_anomZa_ppp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h1
c-----A7h2
c-----A(2,2,2,1,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,2,1,2)
      A7h2=czip
      A7h2 = zajj_tree_qqgg_anomZa_ppm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h2
c-----A7h3
c-----A(2,2,1,2,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,1,2,2)
      A7h3=czip
      A7h3 = zajj_tree_qqgg_anomZa_pmp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h3
c-----A7h4
c-----A(2,1,2,2,2)
      elseif
     .((g1==2).and.(g2==1).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,2,2,2)
      A7h4=czip
      A7h4 = zajj_tree_qqgg_anomZa_mpp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h4
c-----A7h5
c-----A(2,2,1,1,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,1,1,2)
      A7h5=czip
      A7h5 = zajj_tree_qqgg_anomZa_pmm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h5
c-----A7h6
c-----A(2,1,2,1,2) =
      elseif
     .((g1==2).and.(g2==1).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,2,1,2)
      A7h6=czip
      A7h6 = zajj_tree_qqgg_anomZa_mpm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h6
c-----A7h7
c-----A(2,1,1,2,2) =
      elseif
     .((g1==2).and.(g2==1).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,1,2,2)
      A7h7=czip
      A7h7 = zajj_tree_qqgg_anomZa_mmp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h7
c-----A7h8
c-----A(2,1,1,1,2)
      elseif
     .((g1==2).and.(g2==1).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,1,1,2)
      A7h8=czip
      A7h8 = zajj_tree_qqgg_anomZa_mmm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomZa = A7h8
c----
      else
      amp_qqagg_anomZa = czip
      call abort
      endif

      end

      function amp_qqagg_anomaZ(i1,qh,i2,h2,i3,h3,i4,h4,
     & i5,lh,j6,j7,za,zb,ha,hb)
        use zajj_treeamps_m
        implicit none
        !include 'types.f'
        complex(dp):: amp_qqagg_anomaZ
        
        include 'constants.f'
        include 'nf.f'
        !include 'mxpart.f'
        include 'cplx.h'

        integer, intent(in) :: i1,i2,i3,i4,i5,j6,j7
        integer, intent(in) :: qh,h2,h3,h4,lh
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        integer:: i6,i7,j,k,g1,g2,g3,g4,lg
        complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)
        complex(dp):: A7h1,A7h2,A7h3,A7h4,A7h5,A7h6,A7h7,A7h8
c-----
C---A(qh,h2,h3,h4,lh)
C---h=1 LH
C---h=2 RH
        g1=qh
        g2=h2
        g3=h3
        g4=h4
        lg=lh
        i6=j6
        i7=j7
        if (g1*lg==2)then
          i7=j6
          i6=j7
          lg=3-lg
        endif

        if (g1==2) then
          xa(:,:) = za(:,:)
          xb(:,:) = zb(:,:)
        elseif (g1==1) then
        g1=2
        g2=3-g2
        g3=3-g3
        g4=3-g4
        lg=3-lg

        xa(:,:)=zb(:,:)
        xb(:,:)=za(:,:)
        endif

c-----A7h1
c-----A(2,2,2,2,2)
      if (
     & (g1==2).and.(g2==2).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,2,2,2)
      A7h1=czip
      A7h1 = zajj_tree_qqgg_anomaZ_ppp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h1
c-----A7h2
c-----A(2,2,2,1,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,2,1,2)
      A7h2=czip
      A7h2 = zajj_tree_qqgg_anomaZ_ppm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h2
c-----A7h3
c-----A(2,2,1,2,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,1,2,2)
      A7h3=czip
      A7h3 = zajj_tree_qqgg_anomaZ_pmp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h3
c-----A7h4
c-----A(2,1,2,2,2)
      elseif
     .((g1==2).and.(g2==1).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,2,2,2)
      A7h4=czip
      A7h4 = zajj_tree_qqgg_anomaZ_mpp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h4
c-----A7h5
c-----A(2,2,1,1,2) =
      elseif
     .((g1==2).and.(g2==2).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,1,1,2)
      A7h5=czip
      A7h5 = zajj_tree_qqgg_anomaZ_pmm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h5
c-----A7h6
c-----A(2,1,2,1,2) =
      elseif
     .((g1==2).and.(g2==1).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,2,1,2)
      A7h6=czip
      A7h6 = zajj_tree_qqgg_anomaZ_mpm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h6
c-----A7h7
c-----A(2,1,1,2,2) =
      elseif
     .((g1==2).and.(g2==1).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,1,2,2)
      A7h7=czip
      A7h7 = zajj_tree_qqgg_anomaZ_mmp(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h7
c-----A7h8
c-----A(2,1,1,1,2)
      elseif
     .((g1==2).and.(g2==1).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,1,1,2)
      A7h8=czip
      A7h8 = zajj_tree_qqgg_anomaZ_mmm(i5,i1,i6,i7,i2,i3,i4,xa,xb,ha,hb)
      amp_qqagg_anomaZ = A7h8
c----
      else
      amp_qqagg_anomaZ = czip
      call abort
      endif

      end


      ! this is the xzqqagg_ql function with added arguments for the
      ! amplitude function and anomalous couplings
      subroutine xzqqagg_anom(j1,j2,j3,j4,j5,j6,j7,za,zb,ampFun,ha,hb,a70h3)
        include 'types.f'
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'cplx.h'

        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb
        complex(dp), intent(out) :: a70h3(2,2,2,2,2,2)

        integer :: lh,h2,h3,h4,hq
        integer :: hqc,lhc,h2c,h3c,h4c
        complex(dp) :: amp_qqagg_ql

        interface 
          function ampFun(i1,qh,i2,h2,i3,h3,i4,h4,i5,lh,j6,j7,za,zb,ha,hb)
            implicit none
            include 'types.f'
            include 'mxpart.f'
            complex(dp) :: ampFun
            integer, intent(in) :: i1,i2,i3,i4,i5,j6,j7
            integer, intent(in) :: qh,h2,h3,h4,lh
            complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
            complex(dp), intent(in) :: ha,hb
          end function
        end interface


        a70h3 = 0._dp

        hq = 2
        hqc = 1

        do lh=1,2; do h2=1,2; do h3=1,2; do h4=1,2
          a70h3(1,hq,h2,h3,h4,lh) =
     &          ampFun(j1,hq,j2,h2,j3,h3,j4,h4,j5,lh,j6,j7,za,zb,ha,hb)
          a70h3(2,hq,h2,h3,h4,lh) =
     &          ampFun(j1,hq,j2,h2,j4,h4,j3,h3,j5,lh,j6,j7,za,zb,ha,hb)

          lhc=mod(lh,2)+1
          h2c=mod(h2,2)+1
          h3c=mod(h3,2)+1
          h4c=mod(h4,2)+1

          a70h3(1,hqc,h2c,h3c,h4c,lhc) = conjg(a70h3(1,hq,h2,h3,h4,lh))
          a70h3(2,hqc,h2c,h3c,h4c,lhc) = conjg(a70h3(2,hq,h2,h3,h4,lh))
! manually fix signs due to crossing that are ruined by conjg relation above
          if ((j4 == 4) .and. (j3 /= 3)) then
            a70h3(1,hqc,h2c,h3c,h4c,lhc)=-a70h3(1,hqc,h2c,h3c,h4c,lhc)
            a70h3(2,hqc,h2c,h3c,h4c,lhc)=-a70h3(2,hqc,h2c,h3c,h4c,lhc)
          endif
        enddo; enddo; enddo; enddo
        
c       write (*,*) "lepton plus"
c       write (*,*) "anomzz ppp", a70h3(1, 2,2,2,2,2)/zajj_tree_qqgg_anomzz_ppp(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mpp", a70h3(1, 2,1,2,2,2)/zajj_tree_qqgg_anomzz_mpp(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "anomzz ppm", a70h3(1, 2,2,2,1,2)/zajj_tree_qqgg_anomzz_ppm(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mpm", a70h3(1, 2,1,2,1,2)/zajj_tree_qqgg_anomzz_mpm(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "anomzz pmm", a70h3(1, 2,2,1,1,2)/zajj_tree_qqgg_anomzz_pmm(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mmm", a70h3(1, 2,1,1,1,2)/zajj_tree_qqgg_anomzz_mmm(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "anomzz pmp", a70h3(1, 2,2,1,2,2)/zajj_tree_qqgg_anomzz_pmp(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mmp", a70h3(1, 2,1,1,2,2)/zajj_tree_qqgg_anomzz_mmp(j5,j1,j6,j7,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "lepton minus"
c       write (*,*) "anomzz ppp", a70h3(1, 2,2,2,2,1)/zajj_tree_qqgg_anomzz_ppp(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mpp", a70h3(1, 2,1,2,2,1)/zajj_tree_qqgg_anomzz_mpp(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "anomzz ppm", a70h3(1, 2,2,2,1,1)/zajj_tree_qqgg_anomzz_ppm(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mpm", a70h3(1, 2,1,2,1,1)/zajj_tree_qqgg_anomzz_mpm(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "anomzz pmm", a70h3(1, 2,2,1,1,1)/zajj_tree_qqgg_anomzz_pmm(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mmm", a70h3(1, 2,1,1,1,1)/zajj_tree_qqgg_anomzz_mmm(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "anomzz pmp", a70h3(1, 2,2,1,2,1)/zajj_tree_qqgg_anomzz_pmp(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)
c       write (*,*) "anomzz mmp", a70h3(1, 2,1,1,2,1)/zajj_tree_qqgg_anomzz_mmp(j5,j1,j7,j6,j2,j3,j4,za,zb,ha,hb)

c       write (*,*) "quark Minus, za <-> zb"

c       write (*,*) "anomZZ MMM", a70h3(1, 1,1,1,1,1)/zajj_tree_qqgg_anomZZ_ppp(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PMM", a70h3(1, 1,2,1,1,1)/zajj_tree_qqgg_anomZZ_mpp(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "anomZZ MMP", a70h3(1, 1,1,1,2,1)/zajj_tree_qqgg_anomZZ_ppm(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PMP", a70h3(1, 1,2,1,2,1)/zajj_tree_qqgg_anomZZ_mpm(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "anomZZ MPP", a70h3(1, 1,1,2,2,1)/zajj_tree_qqgg_anomZZ_pmm(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PPP", a70h3(1, 1,2,2,2,1)/zajj_tree_qqgg_anomZZ_mmm(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "anomZZ MPM", a70h3(1, 1,1,2,1,1)/zajj_tree_qqgg_anomZZ_pmp(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PPM", a70h3(1, 1,2,2,1,1)/zajj_tree_qqgg_anomZZ_mmp(j5,j1,j6,j7,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "quark Minus, lepton plus, za <-> zb"

c       write (*,*) "anomZZ MMM", a70h3(1, 1,1,1,1,2)/zajj_tree_qqgg_anomZZ_ppp(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PMM", a70h3(1, 1,2,1,1,2)/zajj_tree_qqgg_anomZZ_mpp(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "anomZZ MMP", a70h3(1, 1,1,1,2,2)/zajj_tree_qqgg_anomZZ_ppm(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PMP", a70h3(1, 1,2,1,2,2)/zajj_tree_qqgg_anomZZ_mpm(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "anomZZ MPP", a70h3(1, 1,1,2,2,2)/zajj_tree_qqgg_anomZZ_pmm(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PPP", a70h3(1, 1,2,2,2,2)/zajj_tree_qqgg_anomZZ_mmm(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)

c       write (*,*) "anomZZ MPM", a70h3(1, 1,1,2,1,2)/zajj_tree_qqgg_anomZZ_pmp(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)
c       write (*,*) "anomZZ PPM", a70h3(1, 1,2,2,1,2)/zajj_tree_qqgg_anomZZ_mmp(j5,j1,j7,j6,j2,j3,j4,zb,za,ha,hb)


      end subroutine

      end module
