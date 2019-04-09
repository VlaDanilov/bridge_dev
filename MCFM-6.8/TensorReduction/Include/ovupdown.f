      double precision up(4,4),dn(4,4)
      data up/-1d0, 0d0, 0d0, 0d0,
     &         0d0,-1d0, 0d0, 0d0,
     &         0d0, 0d0,-1d0, 0d0,
     &         0d0, 0d0, 0d0,+1d0/
      data dn/+1d0, 0d0, 0d0, 0d0,
     &         0d0,+1d0, 0d0, 0d0,
     &         0d0, 0d0,+1d0, 0d0,
     &         0d0, 0d0, 0d0,+1d0/
      save up,dn

