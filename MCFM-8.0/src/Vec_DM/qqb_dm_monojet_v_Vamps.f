
      subroutine qqb_dm_monojet_v_Vamps(p,i1,i2,i3,i4,i5,qgqb) 
      implicit none
      include 'types.f'
       
!----- combined colour and Born intefered amplitudes 
!----- as a function of quark line helicity 
!------ default is q(i1)+g(i2)+qb(i3)+x(i4)+x(i5) 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'qcdcouple.f' 
!      include 'zprods_decl.f'
      real(dp):: qgqb(2) 
      integer:: i1,i2,i3,i4,i5 
      real(dp):: p(mxpart,4) 
      integer:: h1,h2,h3,h4 
      complex(dp):: amp_tree(2,2,2,2),amp_lc(2,2,2,2),
     & amp_slc(2,2,2,2) 

!------ tree
      call qqb_dm_monojet_Vamps(p,i1,i2,i3,i4,i5,amp_tree)  
!------ leading colour      
       call qqb_dm_monojet_lc_Vamps(p,i1,i2,i3,i4,i5,amp_lc) 
!------ subleading colour (swap i2 and i3) 
      call qqb_dm_monojet_slc_Vamps(p,i1,i3,i2,i4,i5,amp_slc)

    
      qgqb(:)=zero 
      do h1=1,2 
         do h2=1,2
            do h3=1,2 
               do h4=1,2 
                  qgqb(h1)=qgqb(h1)
     &                 +ason2pi*Dble(conjg(amp_tree(h1,h2,h3,h4))*
     &                 (amp_lc(h1,h2,h3,h4)
     &                 +amp_slc(h1,h2,h3,h4)/xnsq))
               enddo
            enddo
         enddo
      enddo


      
      return 
      end 

      
