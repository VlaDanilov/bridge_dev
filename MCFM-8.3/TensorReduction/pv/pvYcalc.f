      subroutine pvYcalc(q1,q2,q3,m1s,m2s,m3s,m4s)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'TRconstants.f'
      integer:: nu
      real(dp):: q1(4),q2(4),q3(4),q4(4),m1s,m2s,m3s,m4s,Y(4,4),pvdot

      do nu=1,4
      q4(nu)=0._dp
      enddo

      Y(1,1)=m1s
      Y(2,2)=m2s
      Y(3,3)=m3s
      Y(4,4)=m4s

      Y(1,2)=half*(m1s+m2s-pvdot(q4,q4)-pvdot(q1,q1)+two*pvdot(q1,q4))
      Y(1,3)=half*(m1s+m3s-pvdot(q4,q4)-pvdot(q2,q2)+two*pvdot(q2,q4))
      Y(1,4)=half*(m1s+m4s-pvdot(q4,q4)-pvdot(q3,q3)+two*pvdot(q3,q4))

      Y(2,3)=half*(m2s+m3s-pvdot(q1,q1)-pvdot(q2,q2)+two*pvdot(q1,q2))
      Y(2,4)=half*(m2s+m4s-pvdot(q1,q1)-pvdot(q3,q3)+two*pvdot(q1,q3))

      Y(3,4)=half*(m3s+m4s-pvdot(q2,q2)-pvdot(q3,q3)+two*pvdot(q2,q3))

c      write(6,*) 'row1',Y(1,1),Y(1,2),Y(1,3),Y(1,4)
c      write(6,*) 'row2',Y(1,2),Y(2,2),Y(2,3),Y(2,4)
c      write(6,*) 'row3',Y(1,3),Y(2,3),Y(3,3),Y(3,4)
c      write(6,*) 'row4',Y(1,4),Y(2,4),Y(3,4),Y(4,4)

      if (qlzero(Y(1,3))  
     .  .and. qlnonzero(Y(1,2)) .and. qlnonzero(Y(3,4))) then
      write(6,*) 'swapping 1 and 2'
      call pvswap(q1,q2,m2s,m3s)

      elseif (qlzero(Y(2,4))  
     . .and. qlnonzero(Y(1,4)) .and. qlnonzero(Y(2,3))) then
      call pvswap(q2,q3,m3s,m4s)
      write(6,*) 'swapping 2 and 3'

      elseif (qlzero(Y(2,4))  
     . .and. qlnonzero(Y(1,2)) .and. qlnonzero(Y(3,4))) then
      call pvswap(q2,q3,m3s,m4s)
      call pvswap(q1,q2,m2s,m3s)
      write(6,*) '231 --> 123'
      endif

c      Y(1,2)=half*(m1s+m2s-pvdot(q4,q4)-pvdot(q1,q1)+two*pvdot(q1,q4))
c      Y(1,3)=half*(m1s+m3s-pvdot(q4,q4)-pvdot(q2,q2)+two*pvdot(q2,q4))
c      Y(1,4)=half*(m1s+m4s-pvdot(q4,q4)-pvdot(q3,q3)+two*pvdot(q3,q4))

c      Y(2,3)=half*(m2s+m3s-pvdot(q1,q1)-pvdot(q2,q2)+two*pvdot(q1,q2))
c      Y(2,4)=half*(m2s+m4s-pvdot(q1,q1)-pvdot(q3,q3)+two*pvdot(q1,q3))

c      Y(3,4)=half*(m3s+m4s-pvdot(q2,q2)-pvdot(q3,q3)+two*pvdot(q2,q3))

c      write(6,*) 'new1',Y(1,1),Y(1,2),Y(1,3),Y(1,4)
c      write(6,*) 'new2',Y(1,2),Y(2,2),Y(2,3),Y(2,4)
c      write(6,*) 'new3',Y(1,3),Y(2,3),Y(3,3),Y(3,4)
c      write(6,*) 'new4',Y(1,4),Y(2,4),Y(3,4),Y(4,4)


      return
      end
