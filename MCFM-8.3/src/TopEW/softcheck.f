      subroutine softcheck(p)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      integer j,k
      real(dp):: p(mxpart,4),msqg(-nf:nf,-nf:nf),
     .     msqs(-nf:nf,-nf:nf)

      call writeout(p)

      write(6,*) 'j, k, msqg, msqs, g/s: '

      call qqb_QQb_mix_g(p,msqg)
      call qqb_QQb_mix_sft(p,msqs)


      do j=-nf,nf
         k = -j
         if (j .ne. 0) then
            write(6,*) j, k, 
     .           msqg(j,k),msqs(j,k),msqg(j,k)/msqs(j,k)
         end if
      end do

      pause

      end subroutine softcheck
