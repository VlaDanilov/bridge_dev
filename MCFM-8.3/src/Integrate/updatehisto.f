      subroutine updatehisto(part,ips,ncall)
      implicit none
      include 'types.f'
      include 'histo.f'
      include 'accumhist.f'
      include 'tiny.f'
      integer(kind(lord)), intent(in) :: part
      integer, intent(in) :: ips
      integer(kind=8), intent(in) :: ncall
      integer :: n,j
      real(dp):: thiswgt,sumprev
      integer nplotmax
      common/nplotmax/nplotmax

c--- compute errors from n, maxhisto+n and store in 2*maxhisto+n
      do n=1,nplotmax
        !call mopera(n,'U',maxhisto+n,2*maxhisto+n,1._dp,1._dp,ncall)
        hist(2*maxhisto+n,:) = sqrt(abs(hist(maxhisto+n,:)/hdel(n) - hist(n,:)**2/ncall))
      enddo


      ! we are accumulating different parts separately
      do n=1,nplotmax
        do j=1,nbin(n)
          ! assuming we don't integrate something that results in a
          ! tiny error, like maybe a constant?
          
          if (abs(hist(2*maxhisto+n,j)) > 100*tiny) then
c            print*, " updatehisto.f:updatehisto: " 
c            print*, " accumhist(part,ips,n,j) 
c     &      =", accumhist(part,ips,n,j)
c            print*, " hist(maxhisto+n,j)", maxhisto+n,j,hist(maxhisto+n,
c     &      j)
c            print*, " hist(maxhisto+n,j)/hdel(n) = ",
c     &      hist(maxhisto+n,j)/hdel(n)
c            print*, " hist(n,j)**2/ncall = ", hist(n,j)**2/ncall
c            print*, " hist(n,j) = ", n, j, 
c     &      hist(n,j)
c            print*, " hist(2*maxhisto+n,j) = ",
c     &      2*maxhisto+n, j, hist(2*maxhisto+n,j) 
c            print*, " hist(n,j)/hist(2*maxhist
c     &      o+n, j)**2 = ", hist(n,j)/hist(2*maxhisto+n, j)**2
c            print*, " part = ", part
c            print*, " ips = ", ips            
            accumhist(part,ips,n,j) =  accumhist(part,ips,n,j)
     &            + hist(n,j)/hist(2*maxhisto+n, j)**2
c           print*, " accumhist(part,ips,n,j)
c    &            + hist(n,j)/hist(2*maxhisto+n, j)**2
c     &      = ", accumhist(part,ips,n,j)
c            print*,"............................................"
c            print*, " accumhist(part,ips,
c     &      2*maxhisto+n,j)=accumhist(part,ips,2*maxhisto+n,j)
c     &      + 1/hist(2*maxhisto+n,j)**2"
c            print*, "accumhist(part,ips,2*maxhisto+n,j) =
c     &      ", accumhist(part,ips,2*maxhisto+n,j)
c            print*," accumhist(part,ips,
c     &       2*maxhisto+n,j) =", accumhist(part,ips,2*maxhisto+n,j)
            accumhist(part,ips,2*maxhisto+n,j) = 
     &                 accumhist(part,ips,2*maxhisto+n,j)
     &                  + 1/hist(2*maxhisto+n, j)**2
c            print*, " 1/hist(2*maxhisto+n,j)**
c     &      2 = ", 1/hist(2*maxhisto+n,j)**2
c            print*," accumhist(part,ips,
c     &       2*maxhisto+n,j)=", accumhist(part,ips,2*maxhisto+n,j)
c            print*, "-----------updatehisto end-----------------"
          endif
        enddo
      enddo

      ! reset for next iteration
      ihis(:,:) = 0
      hist(:,:) = 0._dp
      ient(:) = 0
      iuscore(:) = 0
      ioscore(:) = 0

      end subroutine
      
      subroutine finalizehist()
        implicit none
        include 'types.f'
        include 'histo.f'
        include 'accumhist.f'
        include 'tiny.f'
        integer :: part,n,j,ips
        integer :: nplotmax
        common/nplotmax/nplotmax

        finalhist(:,:) = 0._dp

        do n=1,nplotmax
          do j=1,nbin(n)
            do ips=1,maxIps
              do part=1,maxParts
                if (accumhist(part,ips,2*maxhisto+n,j) > 100*tiny) then
c                  print*, "finalhist(n,j) = finalhist(n,j)+accumhist(
c     &            part,ips,n,j)/accumhist(part,ips,2*maxhisto+n,j)"
c                  print*, "finalhist(n,j) = ", finalhist(n,j)
c                  print*, "updatehisto.f:finalizehist accumhist(part,
c     &            ips,n,j) = ", part, ips, n, j, accumhist(part,ips,n,j)
c                  print*, "updatehisto.f:finalizehist accumhist(part,
c     &            ips,2*maxhisto+n,j) = ", part, ips, n, j, 
c     &            accumhist(part,ips,2*maxhisto+n,j)
                  finalhist(n,j) = finalhist(n,j) + 
     &                  accumhist(part,ips,n,j)
     &                    /accumhist(part,ips,2*maxhisto+n,j)
                  finalhist(2*maxhisto+n,j) = finalhist(2*maxhisto+n,j) +
     &                      1._dp/accumhist(part,ips,2*maxhisto+n,j)
c                 print*, "updatehisto.f:finalizehist FINAL 
c     &            finalhist(n,j) = ", n, j , finalhist(n,j)
c                  print*, "--------------------------------------------"
                endif
              enddo
            enddo
            finalhist(2*maxhisto+n,j) = sqrt(finalhist(2*maxhisto+n,j))
          enddo
        enddo

      end subroutine

      subroutine zerohist()
        implicit none
        include 'types.f'
        include 'histo.f'

        hist(:,:) = 0._dp
        ihis(:,:) = 0
      end subroutine

      subroutine zeroaccumhist()
        implicit none
        include 'types.f'
        include 'histo.f'
        include 'accumhist.f'

        accumhist(:,:,:,:) = 0._dp
      end subroutine

      integer function hist_readsnapshot()
        implicit none
        include 'types.f'
        include 'histo.f'
        include 'accumhist.f'

        character*255 runname
        common/runname/runname

        integer :: ierr
        integer :: i,j,k,l

        ! How do I make sure that really the complete record length has been read?
        open(unit=11, file=trim(runname)//"_histogram_snapshot.dat",
     &        status='OLD', iostat=ierr,
     &        form='unformatted', access='direct', recl=maxParts*maxIps*nplot*maxnbin*8)
          if (ierr == 0) then
            read(11, rec=1, iostat=ierr) accumhist
            close(unit=11)
            if (ierr /= 0) then
              write (*,*) "Broken histogram_snapshot.dat!"
              call abort
            endif
          endif

          if (ierr == 0) then
            hist_readsnapshot = 0
          else
            hist_readsnapshot = ierr
          endif
      end function

      integer function hist_writesnapshot()
        implicit none
        include 'types.f'
        include 'histo.f'
        include 'accumhist.f'

        character*255 runname
        common/runname/runname

        integer :: ierr
        integer :: i,j,k,l

        open(unit=11, file=trim(runname)//"_histogram_snapshot.dat",
     &        status='REPLACE', iostat=ierr,
     &        form='unformatted', access='direct', recl=maxParts*maxIps*nplot*maxnbin*8)
          if (ierr == 0) then
            write(11, rec=1) accumhist
            close(unit=11)
          endif

          if (ierr == 0) then
            hist_writesnapshot = 0
          else
            hist_writesnapshot = ierr
          endif
      end function
