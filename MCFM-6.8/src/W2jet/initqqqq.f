      subroutine initqqqq(za,zb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'basic.f'
      integer j1,j2,j3,j4,qp(4)
      double complex bll
      data qp/1,2,5,6/
      save qp

      do j1=1,4
      do j2=1,4
      do j3=1,4
      do j4=1,4
        Lla(j1,j2,j3,j4)=0d0
        Lal(j1,j2,j3,j4)=0d0
        Rla(j1,j2,j3,j4)=0d0
        Ral(j1,j2,j3,j4)=0d0
      enddo
      enddo
      enddo
      enddo

      do j1=1,4
       do j2=1,4
        do j3=1,4
         do j4=1,4
           if ((j2.eq.j1).or.(j3.eq.j1).or.(j3.eq.j2)
     &     .or.(j4.eq.j1).or.(j4.eq.j2).or.(j4.eq.j3)) go to 1 
           Lla(j1,j2,j3,j4)=bll(qp(j1),qp(j2),qp(j3),qp(j4),3,4,za,zb)
           Lal(j3,j4,j1,j2)=-Dconjg(Lla(j1,j2,j3,j4))

           Rla(j1,j2,j3,j4)=+Dconjg(Lla(j1,j2,j3,j4))
           Ral(j3,j4,j1,j2)=-Lla(j1,j2,j3,j4)
 1         continue
         enddo
        enddo
       enddo
      enddo

      return
      end

      double complex function bll(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5,i6
c--- note that the quarks are labelled by i1 ... i4 and the
c--- W boson leptons by i5 and i6
c--- this means that i1 ... i4 take the values 1,2,5 and 6
c--- whilst i5=3, i6=4

C     cf Giele and Glover, PRD46 (92) 2008, Eqn A38.

c                    i5(lepton)                i5 
c                    |__i6(antilepton)         |__i6
c                   /                         /
c     i1___\_______/______i3          i1_____/________\_____i3
c          /    o                               o     /
c               o                +              o
c               o                               o
c    i2____\____o__________i4         i2________o_____\_____i4
c          /                                          /

      bll=
     & -four*za(i3,i5)*zb(i1,i2)
     & *(zb(i1,i6)*za(i1,i4)+zb(i2,i6)*za(i2,i4))
     & /(s(i2,i4)*s(i5,i6)*(s(i1,i2)+s(i1,i4)+s(i2,i4)))
     & +four*za(i3,i4)*zb(i1,i6)
     & *(zb(i3,i2)*za(i3,i5)+zb(i4,i2)*za(i4,i5))
     & /(s(i2,i4)*s(i5,i6)*(s(i3,i2)+s(i3,i4)+s(i2,i4)))

      return
      end
      
