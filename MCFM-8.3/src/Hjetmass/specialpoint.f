      subroutine specialpoint(boxints,triints,bubints,amp)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qqbggnames.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      real(dp)::p(mxpart,4),msq
      complex(dp)::rat(2,2),amp(2,2)
      complex(dp)::boxints(Nboxints),dcoeff(2,2,Nboxcoeffs)
      complex(dp)::triints(Ntriints),ccoeff(2,2,Ntricoeffs)
      complex(dp)::bubints(Nbubints),bcoeff(2,2,Nbubcoeffs)
      integer::h3,h4,j
      integer::i1,i2,i3,i4

      msq=mt**2

      bcoeff(2,2,b124)=
     & - 6205723._dp/115031826._dp
     & - 5725873._dp/57515913._dp*im
      bcoeff(2,2,b123)=
     & + 277804._dp/90261925._dp
     & - 7510572._dp/90261925._dp*im
      bcoeff(2,2,b12)=
     & + 21869._dp/1314950._dp
     & + 14379._dp/657475._dp*im
      bcoeff(2,2,b34)=0._dp

      bcoeff(2,2,b1234)=
     & + 1351778026._dp/39480565995._dp
     & + 6352085732._dp/39480565995._dp*im

      ccoeff(2,2,c3x4)=0._dp

      ccoeff(2,2,c3x124m2)=
     & - 217160699._dp/1290304080._dp
     & + 118792591._dp/645152040._dp*im
      ccoeff(2,2,c3x124)=
     & + 7758842323._dp/1624827360._dp
     & + 747586393._dp/812413680._dp*im
      ccoeff(2,2,c3x12)=
     & - 25548._dp/63869._dp
     & - 48336._dp/63869._dp*im
      ccoeff(2,2,c4x123m2)=
     & - 23922007._dp/185182530._dp
     & - 2274729._dp/30863755._dp*im
      ccoeff(2,2,c4x123)=
     & + 50682289._dp/10155171._dp
     & + 19168041._dp/3385057._dp*im
      ccoeff(2,2,c4x12)=
     & - 178836._dp/319345._dp
     & - 338352._dp/319345._dp*im
      ccoeff(2,2,c12x34m2)=
     & - 217259._dp/995605._dp
     & - 240138._dp/995605._dp*im
      ccoeff(2,2,c12x34)=
     & + 346813._dp/117130._dp
     & + 141483._dp/58565._dp*im
      dcoeff(2,2,d3x4x12m2)=
     & - 219097._dp/58565._dp
     & - 34854._dp/58565._dp*im
      dcoeff(2,2,d3x4x12)=
     & + 340343._dp/6890._dp
     & + 2913._dp/3445._dp*im
      dcoeff(2,2,d3x12x4m2)=
     & + 29227._dp/18785._dp
     & + 67914._dp/18785._dp*im
      dcoeff(2,2,d3x12x4)=
     & + 3647916._dp/319345._dp
     & + 2338812._dp/319345._dp*im
      dcoeff(2,2,d12x3x4m2)=
     & + 1838._dp/58565._dp
     & - 205284._dp/58565._dp*im
      dcoeff(2,2,d12x3x4)=
     & + 647._dp/689._dp
     & + 27714._dp/689._dp*im
      rat(2,2)=
     &  - 4578853._dp/17757792._dp
     & - 580543._dp/8878896._dp*im

      bcoeff(1,1,b124)=
     & + 17882701._dp/575159130._dp
     & + 6017759._dp/287579565._dp*im
      bcoeff(1,1,b123)=
     & + 1513006._dp/90261925._dp
     & + 4127508._dp/90261925._dp*im
      bcoeff(1,1,b12)=
     & - 6029._dp/1314950._dp
     & - 17811._dp/657475._dp*im
      bcoeff(1,1,b34)=0._dp
      bcoeff(1,1,b1234)=
     & - 1708291493._dp/39480565995._dp
     & - 1561995674._dp/39480565995._dp*im
      ccoeff(1,1,c3x4)=0._dp
      ccoeff(1,1,c3x124m2)=
     & + 174166633._dp/1290304080._dp
     & + 146108147._dp/645152040._dp*im
      ccoeff(1,1,c3x124)=
     & - 4539405473._dp/1624827360._dp
     & - 4989765307._dp/812413680._dp*im
      ccoeff(1,1,c3x12)=
     & + 11854._dp/63869._dp
     & + 16572._dp/63869._dp*im
      ccoeff(1,1,c4x123m2)=
     & - 40727299._dp/185182530._dp
     & + 8099253._dp/30863755._dp*im
      ccoeff(1,1,c4x123)=
     & + 214486183._dp/101551710._dp
     & - 118745451._dp/16925285._dp*im
      ccoeff(1,1,c4x12)=
     & + 82978._dp/319345._dp
     & + 116004._dp/319345._dp*im
      ccoeff(1,1,c12x34m2)=
     & + 16207._dp/199121._dp
     & + 62706._dp/199121._dp*im
      ccoeff(1,1,c12x34)=
     & - 20109._dp/23426._dp
     & - 65811._dp/11713._dp*im
      dcoeff(1,1,d3x4x12m2)=
     & + 177913._dp/58565._dp
     & + 132534._dp/58565._dp*im
      dcoeff(1,1,d3x4x12)=
     & - 251501._dp/6890._dp
     & - 158859._dp/3445._dp*im
      dcoeff(1,1,d3x12x4m2)=
     & - 17679._dp/18785._dp
     & + 2478._dp/18785._dp*im
      dcoeff(1,1,d3x12x4)=
     & - 509459._dp/319345._dp
     & - 9117462._dp/319345._dp*im
      dcoeff(1,1,d12x3x4m2)=
     & - 96878._dp/58565._dp
     & + 180996._dp/58565._dp*im
      dcoeff(1,1,d12x3x4)=
     & + 5806._dp/265._dp
     & - 13092._dp/265._dp*im
      rat(1,1)=
     &  - 157901._dp/88788960._dp
     & + 17842241._dp/44394480._dp*im
!

      bcoeff(1,2,b124)=
     & + 332840516._dp/2795193765._dp
     & - 271750712._dp/2795193765._dp*im
      bcoeff(1,2,b123)=
     & + 47294196._dp/877320925._dp
     & + 6344928._dp/877320925._dp*im
      bcoeff(1,2,b12)=
     & + 2985321231._dp/105161656600._dp
     & - 1456932201._dp/52580828300._dp*im
      bcoeff(1,2,b34)=
     & + 206192497._dp/883711400._dp
     & - 64022247._dp/441855700._dp*im
      bcoeff(1,2,b1234)=
     & - 527891785077851._dp/1214390053889100._dp
     & + 159443942952901._dp/607195026944550._dp*im
      ccoeff(1,2,c3x4)=
     & + 63807324._dp/9677005._dp
     & - 35983968._dp/9677005._dp*im
      ccoeff(1,2,c3x124m2)=
     & - 89555593._dp/502839090._dp
     & + 4699163._dp/251419545._dp*im
      ccoeff(1,2,c3x124)=
     & - 15099523141._dp/1974109020._dp
     & + 7670797931._dp/987054510._dp*im
      ccoeff(1,2,c3x12)=
     & + 1507728._dp/1935401._dp
     & + 790704._dp/1935401._dp*im
      ccoeff(1,2,c4x123m2)=
     & + 45738289._dp/192444590._dp
     & - 3798999._dp/96222295._dp*im
      ccoeff(1,2,c4x123)=
     & - 2332909523._dp/658036340._dp
     & - 592294407._dp/329018170._dp*im
      ccoeff(1,2,c4x12)=
     & + 41993112._dp/9677005._dp
     & - 35168784._dp/9677005._dp*im
      ccoeff(1,2,c12x34m2)=
     & - 12027659._dp/68286790._dp
     & + 2607789._dp/34143395._dp*im
      ccoeff(1,2,c12x34)=
     & - 3485617322623._dp/398111985700._dp
     & + 1219243038873._dp/199055992850._dp*im
      dcoeff(1,2,d3x4x12m2)=
     & - 1246477._dp/182585._dp
     & + 612414._dp/182585._dp*im
      dcoeff(1,2,d3x4x12)=
     & - 1440596756._dp/9677005._dp
     & + 1476648792._dp/9677005._dp*im
      dcoeff(1,2,d3x12x4m2)=
     & + 102241._dp/58565._dp
     & - 47262._dp/58565._dp*im
      dcoeff(1,2,d3x12x4)=
     & - 2841._dp/442._dp
     & + 951._dp/221._dp*im
      dcoeff(1,2,d12x3x4m2)=
     & + 121162._dp/182585._dp
     & - 120084._dp/182585._dp*im
      dcoeff(1,2,d12x3x4)=
     & - 450928209._dp/9677005._dp
     & - 243871362._dp/9677005._dp*im
      rat(1,2)=
     &  - 377114971._dp/6470495460._dp
     &  + 89918201._dp/3235247730._dp*im

      bcoeff(2,1,b124)=-44560972._dp/2795193765._dp
     & -259206904._dp/2795193765._dp*im
      bcoeff(2,1,b123)=114766374._dp/877320925._dp
     & -604332._dp/877320925._dp*im
      bcoeff(2,1,b12)=10081395773._dp/525808283000._dp
     & -7137214407._dp/262904141500._dp*im
      bcoeff(2,1,b34)=869926547._dp/4418557000._dp
     & -302512473._dp/2209278500._dp*im
      bcoeff(2,1,b1234)=-26121758057349853._dp/78935353502791500._dp
     & +10162835075580127._dp/39467676751395750._dp*im
      ccoeff(2,1,c3x4)=35997976._dp/9677005._dp
     & -35218968._dp/9677005._dp*im
      ccoeff(2,1,c3x124m2)=11303507._dp/38679930._dp
     & +283987._dp/19339965._dp*im
      ccoeff(2,1,c3x124)=1422852461._dp/1974109020._dp
     & +8425318801._dp/987054510._dp*im
      ccoeff(2,1,c3x12)=5925336._dp/1935401._dp
     & +1217952._dp/1935401._dp*im
      ccoeff(2,1,c4x123m2)=-6165179._dp/38488918._dp 
     & +506361._dp/19244459._dp*im
      ccoeff(2,1,c4x123)=-6762028121._dp/658036340._dp
     & -1266225561._dp/329018170._dp*im
      ccoeff(2,1,c4x12)=-2366392._dp/1935401._dp
     & -7505904._dp/1935401._dp*im
      ccoeff(2,1,c12x34m2)=-34084737._dp/341433950._dp
     & +12661683._dp/170716975._dp*im
      ccoeff(2,1,c12x34)=-10099056008573._dp/1990559928500._dp
     & +5984273261607._dp/995279964250._dp*im

      dcoeff(2,1,d3x4x12m2)=706557._dp/182585._dp
     & +709674._dp/182585._dp*im
      dcoeff(2,1,d3x4x12)=66074724._dp/1935401._dp
     & +310467048._dp/1935401._dp*im
      dcoeff(2,1,d3x12x4m2)=104399._dp/58565._dp
     & -42282._dp/58565._dp*im
      dcoeff(2,1,d3x12x4)=-19319._dp/2210._dp
     & +4221._dp/1105._dp*im
      dcoeff(2,1,d12x3x4m2)=-271802._dp/36517._dp
     & -46524._dp/36517._dp*im
      dcoeff(2,1,d12x3x4)=-1111731639._dp/9677005._dp
     & -318996198._dp/9677005._dp*im
      rat(2,1)=521247091._dp/32352477300._dp
     & +931457531._dp/16176238650._dp*im

      do h3=1,2
      do h4=1,2
    
      amp(h3,h4)=im*(
     & +(dcoeff(h3,h4,d3x4x12)+msq*dcoeff(h3,h4,d3x4x12m2))
     & *boxints(d3x4x12)
     & +(dcoeff(h3,h4,d3x12x4)+msq*dcoeff(h3,h4,d3x12x4m2))
     & *boxints(d3x12x4)
     & +(dcoeff(h3,h4,d12x3x4)+msq*dcoeff(h3,h4,d12x3x4m2))
     & *boxints(d12x3x4)

     &+(ccoeff(h3,h4,c3x124)+msq*ccoeff(h3,h4,c3x124m2))*triints(c3x124)
     &+(ccoeff(h3,h4,c4x123)+msq*ccoeff(h3,h4,c4x123m2))*triints(c4x123)
     &+(ccoeff(h3,h4,c12x34)+msq*ccoeff(h3,h4,c12x34m2))*triints(c12x34)
     &+ccoeff(h3,h4,c3x12)*triints(c3x12)
     &+ccoeff(h3,h4,c4x12)*triints(c4x12)
     &+ccoeff(h3,h4,c3x4)*triints(c3x4)

     &+bcoeff(h3,h4,b12)*bubints(b12)
     &+bcoeff(h3,h4,b34)*bubints(b34)
     &+bcoeff(h3,h4,b123)*bubints(b123)
     &+bcoeff(h3,h4,b124)*bubints(b124)
     &+bcoeff(h3,h4,b1234)*bubints(b1234)
     &+rat(h3,h4))
      enddo
      enddo
      return
      end
