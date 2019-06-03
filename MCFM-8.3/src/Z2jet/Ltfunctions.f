      function Ltm1(bx,by,x,y,xInt)
      implicit none
      include 'types.f'
      integer bx,by
      real(dp)::x,y
      complex(dp)::Ltm1,xInt(15)
      
      Ltm1=xInt(by)-xInt(bx)
      
      return
      end
      

      function Lt0(bx,by,x,y,xInt)
      implicit none
      include 'types.f'
      integer bx,by
      real(dp)::x,y
      complex(dp)::Lt0,Ltm1,xInt(15)
      
      Lt0=y/(y-x)*Ltm1(bx,by,x,y,xInt)
      
      return
      end
      

      function Lt1(bx,by,x,y,xInt)
      implicit none
      include 'types.f'
      integer bx,by
      real(dp)::x,y
      complex(dp)::Lt1,Lt0,xInt(15)
      
      Lt1=y/(y-x)*(Lt0(bx,by,x,y,xInt)+1._dp)
      
      return
      end
