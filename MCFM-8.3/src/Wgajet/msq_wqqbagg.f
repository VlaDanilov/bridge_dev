      subroutine wagg_a70h(j1,j2,j3,j4,j5,j6,j7,xa,xb,a70qq,a70ql)
      implicit none
      include 'types.f'

********************************************************************
* 0 -> f(p1) + f(p2) + l(p3) + vb(p4) + gam(p5) + f(p6) + f(p7)    *
* return helicity amplitudes for each channel                      *
********************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: h5,h6,h7,ic
      complex(dp):: a70qq(2,2,2,2),a70ql(2,2,2,2)
      complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)

c-----initialize
      do h7=1,2
      do h6=1,2
      do h5=1,2
      do ic=1,2
      a70qq(ic,h5,h6,h7)=czip
      a70ql(ic,h5,h6,h7)=czip
      enddo
      enddo
      enddo
      enddo

c----compute amplitude
      call xwqqagg_qq(j1,j2,j3,j4,j5,j6,j7,xa,xb,a70qq)
      call xwqqagg_ql(j1,j2,j3,j4,j5,j6,j7,xa,xb,a70ql)

      return 
      end

      subroutine wagg_m70sq(a70qq,a70ql,msq)
      implicit none
      include 'types.f'

********************************************************************
* 0 -> f(p1) + f(p2) + l(p3) + vb(p4) + gam(p5) + f(p6) + f(p7)    *
* return matrix element squared given the                          *
* helicity amplitudes from each channel                            *
********************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'new_pspace.f'
      include 'ipsgen.f'
      real(dp):: msq,m70hsq(8),CF1,CF2,fac
      real(dp):: m70hsqAA(8),m70hsqBB(8),m70hsqAB(8)
      integer:: i,j,k,ihel,ic
      integer:: h5,h6,h7
      complex(dp):: a70qq(2,2,2,2),a70ql(2,2,2,2)
      complex(dp):: m70hA(2,8),m70hB(2,8)

c-----fill m70hA and m70hB
      ihel=1
      do h7=1,2
      do h6=1,2
      do h5=1,2
      do ic=1,2
         m70hA(ic,ihel)=a70qq(ic,h5,h6,h7)
         m70hB(ic,ihel)=a70ql(ic,h5,h6,h7)
      enddo
      ihel=ihel+1
      enddo
      enddo
      enddo

c-----overcount ihel
      ihel=ihel-1

c-----square them up
c-----fac=( (sqrt(2)*e)*(sqrt(2)*g)**2*(sqrt(2)*gw/sqrt(2))**2 )**2
c-----CF1=Nc*CF**2 and CF2=-CF/2

      fac=8._dp*gwsq**2*gsq**2*esq

      CF1=16._dp/three
      CF2=-2._dp/three

      do i=1,8
         m70hsqAA(i)=fac*(CF1*(abs(m70hA(1,i))**2+abs(m70hA(2,i))**2)
     &              +CF2*two*real(m70hA(1,i)*conjg(m70hA(2,i))))
         m70hsqBB(i)=fac*(CF1*(abs(m70hB(1,i))**2+abs(m70hB(2,i))**2)
     &              +CF2*two*real(m70hB(1,i)*conjg(m70hB(2,i))))
         m70hsqAB(i)=fac*(CF1*( 2._dp*real(conjg(m70hA(1,i))*m70hB(1,i))
     &                         +2._dp*real(conjg(m70hA(2,i))*m70hB(2,i)) )
     &              +CF2*( two*real(conjg(m70hA(1,i))*m70hB(2,i))
     &                    +two*real(conjg(m70hB(1,i))*m70hA(2,i)) ))
      enddo

      do i=1,8
c         m70hsq(i)=m70hsqAA(i)+m70hsqBB(i)+m70hsqAB(i)
         if     (ipsgen == 1) then
             m70hsq(i)=m70hsqAA(i)
         elseif (ipsgen == 2) then
             m70hsq(i)=m70hsqBB(i)+m70hsqAB(i)
         else
            write(6,*) 'Parameter ipsgen should be 1 or 2'
            write(6,*) 'ipsgen = ',ipsgen
            stop
         endif
C         if (new_pspace) then
C             m70hsq(i)=m70hsqBB(i)+m70hsqAB(i)
C         else
C             m70hsq(i)=m70hsqAA(i)
C         endif
      enddo

c-----and sum them up
      msq=zip
      do i=1,8
         msq=msq+m70hsq(i)
      enddo

      return
      end

