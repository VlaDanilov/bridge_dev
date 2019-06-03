      subroutine cross(p,i,j,r)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: i,j
      real(dp):: p(mxpart,4),r(3)
      
      r(1)=p(i,2)*p(j,3)-p(j,2)*p(i,3)
      r(2)=p(i,3)*p(j,1)-p(j,3)*p(i,1)
      r(3)=p(i,1)*p(j,2)-p(j,1)*p(i,2)
      
      return
      end
            
      function fphi(n1,n2,p)
      implicit none
      include 'types.f'
      real(dp):: fphi
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: n1,n2
      real(dp):: p(mxpart,4)
    
      fphi=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)
      fphi=fphi/sqrt(p(n1,1)**2+p(n1,2)**2)
      fphi=fphi/sqrt(p(n2,1)**2+p(n2,2)**2)
      if     (fphi > +0.9999999_dp) then
        fphi=0._dp
      elseif (fphi < -0.9999999_dp) then
        fphi=pi
      else
        fphi=acos(fphi)
      endif

      return
      end
                     
      function coslpairet(n1,n2,nm1,nm2,p)
      implicit none
      include 'types.f'
      real(dp):: coslpairet
      include 'nf.f'
      include 'mxpart.f'
      integer:: n1,n2,nm1,nm2,i
      real(dp):: p(mxpart,4),misset(4),pp(4)
      
      if (nm2 == 0) then
        do i=1,4
          misset(i)=p(nm1,i)
        enddo
      else
        do i=1,4
          misset(i)=p(nm1,i)+p(nm2,i)
        enddo
      endif
      
      do i=1,4
        pp(i)=p(n1,i)+p(n2,i)
      enddo
      
      coslpairet=pp(1)*misset(1)+pp(2)*misset(2)
      coslpairet=coslpairet/sqrt(pp(1)**2+pp(2)**2)
      coslpairet=coslpairet/sqrt(misset(1)**2+misset(2)**2)
      
      return
      end
            
      function deltar(i,j,p)
      implicit none
      include 'types.f'
      real(dp):: deltar
      include 'constants.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),phi1,phi2,etarap,dphi
      integer:: i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(p(j,1),p(j,2))
      dphi=phi1-phi2
      if (dphi > pi) dphi=twopi-dphi
      if (dphi < -pi) dphi=twopi+dphi
      deltar=(etarap(i,p)-etarap(j,p))**2+dphi**2
      deltar=sqrt(deltar)
      
      return
      end
      
      subroutine checkmaxhisto(n)
      implicit none
c--- ensure the built-in maximum number of histograms is not exceeded    
      include 'types.f'
      include 'histo.f'

      integer:: n
      
      if (n > maxhisto) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > ',maxhisto,', which is the built-in maximum.'
      write(6,*) 'To use more histograms, change the value of the'
      write(6,*) 'constant MAXHISTO in src/Inc/nplot.f then do:'
      write(6,*)
      write(6,*) ' make clean; make        to recompile from scratch.'
        write(6,*)
      stop
      endif
      
      return
      end
      
