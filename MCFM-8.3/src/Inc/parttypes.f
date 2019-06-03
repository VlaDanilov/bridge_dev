      enum, bind(c)
          enumerator :: nloReal=1, nloVirt=2, nnloBelow=3,
     &          nnloVirtAbove=4, nnloRealAbove=5,
     &          snloBelow = 6, snloAbove = 7, lord=8, nloFrag=9
      endenum
      ! adjust maxParts to be the maximum number used above
      integer, parameter :: maxParts = 9
      integer, parameter :: maxIPS = 2
