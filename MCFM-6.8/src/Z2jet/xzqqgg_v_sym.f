      subroutine xzqqgg_v_sym(mqqb_ax)
      implicit none
************************************************************************
*     Author J.M.Campbell, February 2000                               *
*                                                                      *
*     Supplemental to xzqqgg_v - just calculates the axial piece       *
*     with i1 and i4 swapped wrt that routine                          *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'lc.f'
      integer i1(2),i2(2),i3(2),i4(2),i5(2),i6(2),j,lh,h2,h3,hq,h(2:3)
      double precision fac
      double complex m(2),mqqb_ax(2,2)
      double complex ml1_ax(2),ml2_ax(2)
      double complex a6treeg1,a64ax,a65ax
      character*9 st1(2,2),st3(2,2)
      data i1/1,4/
      data i2/2,3/
      data i3/3,2/
      data i4/4,1/
      data i5/6,5/
      data i6/5,6/
      data st1/'q+g-g-qb-','q+g-g+qb-','q+g+g-qb-','q+g+g+qb-'/
      data st3/'q+qb-g-g-','q+qb-g-g+','q+qb-g+g-','q+qb-g+g+'/
      
      fac=avegg*8d0*gsq**2*esq**2*cf*xn**3*ason2pi
c--- no extra factor here since colour algebra is already done in (2.12)

      do hq=1,2
      do lh=1,2
      mqqb_ax(hq,lh)=0d0
      
      if (colourchoice .le. 1) then
      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
        do j=1,2
        if (hq .eq. 1) then
        m(j)=  a6treeg1(st1(3-h(i2(j)),3-h(i3(j))),
     .     i4(1),i2(j),i3(j),i1(1),i6(lh),i5(lh),zb,za)
c--- note: this symmetry relation (including minus sign) checked numerically 
        ml1_ax(j)=-a64ax(st3(3-h(i2(j)),3-h(i3(j))),
     .     i4(1),i1(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        ml2_ax(j)=-a65ax(st3(3-h(i2(j)),3-h(i3(j))),
     .     i4(1),i1(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        else
        m(j)=  a6treeg1(st1(h(i2(j)),h(i3(j))),
     .     i4(1),i2(j),i3(j),i1(1),i5(lh),i6(lh),za,zb)
        ml1_ax(j)=a64ax(st3(h(i2(j)),h(i3(j))),
     .     i4(1),i1(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        ml2_ax(j)=a65ax(st3(h(i2(j)),h(i3(j))),
     .     i4(1),i1(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        endif
        enddo

      mqqb_ax(hq,lh)=mqqb_ax(hq,lh)+fac/xnsq*(
     .  Dconjg(m(1))*(
     .    (xn-2d0/xn)*ml1_ax(1)-2d0/xn*ml1_ax(2)+one/xn*ml2_ax(1))
     . +Dconjg(m(2))*(
     .    (xn-2d0/xn)*ml1_ax(2)-2d0/xn*ml1_ax(1)+one/xn*ml2_ax(2)))
      
      enddo
      enddo
      endif

      enddo
      enddo

      return
      end
      
