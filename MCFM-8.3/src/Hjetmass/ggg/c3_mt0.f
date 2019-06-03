      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = 0.1D1 / uman
      t4 = cdlogwrap(tman * t1)
      t5 = cdlogwrap(tman * t3)
      t6 = cdlogwrap(uman * t1)
      t7 = tman ** 2
      t8 = cA * t7
      t9 = LogMums + LogMumt + LogMumu
      t10 = sman * uman
      c3mt0 = t7 * ((0.2D1 / 0.3D1) * LogMuMtop + 7 * cA - (0.20D2 / 0.3
     #D1) * tr) + (0.4D1 / 0.3D1) * t7 * ((dilogc(-t1 * t2 + 1) + dilogc
     #(-t2 / tman + 1) + dilogc(-t2 * t3 + 1)) * cA - LogMuMtop * tr) + 
     #t8 * (-pi ** 2 + t4 ** 2 + t5 ** 2 + t6 ** 2) / 3 + (0.11D2 / 0.9D
     #1) * t8 * t9 - (0.20D2 / 0.9D1) * tr * (t7 * t9 + t10) - 2 * t7 * 
     #(-Kg + cF) + (0.2D1 / 0.9D1) * t10 * cA

