      include 'parttypes.f'
      real(dp):: accumhist(maxParts,maxIps,nplot,maxnbin), finalhist(nplot,maxnbin)
      common/accumcommon/accumhist
      common/finalcommon/finalhist

