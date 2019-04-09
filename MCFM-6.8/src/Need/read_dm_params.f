

      subroutine read_dm_params
      implicit none 
      include 'dm_params.f' 
      character*200 line 
      integer j
      open(unit=35,file='dm_parameters.DAT',status='old',err=440) 

      read(35,*) xmass  
      read(35,*) dm_lam
      read(35,*) effective_th
      read(35,*) yukawa_scal 
      read(35,*) line
!-----
      do j=1,5 
         read(35,*) dmL(j) 
      enddo
      read(35,*) line 
!-----
      do j=1,5 
         read(35,*) dmR(j) 
      enddo
      read(35,*) line 
      
      read(35,*) medmass
      read(35,*) medwidth
      read(35,*) g_dmx
      read(35,*) g_dmq


      write(6,*) '************** DM PARAMETERS ****************' 
      write(6,43) '* DM Mass ',xmass,'          *' 
!      write(6,43) '* DM coupling ',gx,'           *' 
      write(6,43) '* Lambda  ',dm_lam,'              *'      
      write(6,*) '* Effective theory : ',effective_th,'  * '
      write(6,44)  '* dm L (d,u,s,c,b),',dmL,' *'
      write(6,44)  '* dm R (d,u,s,c,b),',dmR,' *'
      if(effective_th.eqv..false.) then 
         write(6,*) 'Mediator Mass,',medmass
         write(6,*) 'Mediator Width,',medwidth
         write(6,*) 'g_x ',g_dmx
         write(6,*) 'g_q ',g_dmq
      endif
 !     write(6,43) '* Mediator Mass ',medmass,'        *' 
      write(6,*) '***************************************' 

 43   format(1x,a10,f8.3,a20) 

 44   format(1x,a20,5(f8.4,' '),a4)
      return 
 440  write(6,*) 'Couldnt open dm_parameters.DAT' 
      stop 
      
      return 
      end 
