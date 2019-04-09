      subroutine TRsettensorcontrol(a)
      implicit none
      logical first
      integer a
      include 'TRtensorcontrol.f'
      data first/.true./
      save first
c--- This routine sets the value of the integrer (TRtensorcontrol)
c--- that determines how the tensor integrals are calculated

c--- a = 1 :   OV reduction, followed by PV reduction if necessary
c--- a = 2 :   OV reduction only
c--- a = 3 :   PV reduction only

      TRtensorcontrol=a      
      
      if (first) then
        if     (TRtensorcontrol .eq. 1) then
          write(6,*) '>>> TensorReduction: OV and PV routines'
        elseif (TRtensorcontrol .eq. 2) then
          write(6,*) '>>> TensorReduction: OV routines only'
        elseif (TRtensorcontrol .eq. 3) then
          write(6,*) '>>> TensorReduction: PV routines only'
        else
          write(6,*) 'Invalid call to pvtensorcontrol'
          stop 
        endif   
        first=.false.  
      endif
      
      if ((a .eq. 1) .or. (a .eq. 2)) then
        doovred=.true.
        dopvred=.false.
      elseif (a .eq. 3) then
        doovred=.false.
        dopvred=.true.
      endif
      
      return
      end
      
