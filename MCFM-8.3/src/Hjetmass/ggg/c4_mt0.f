      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = 0.1D1 / uman
      t4 = cdlogwrap(tman * t1)
      t5 = cdlogwrap(tman * t3)
      t6 = cdlogwrap(uman * t1)
      t7 = LogMums + LogMumt + LogMumu
      t8 = uman ** 2
      t9 = sman * tman
      t10 = cA * t8
      c4mt0 = t8 * ((0.2D1 / 0.3D1) * LogMuMtop + 7 * cA - (0.20D2 / 0.3
     #D1) * tr) - (0.20D2 / 0.9D1) * tr * (t7 * t8 + t9) + t10 * (-pi **
     # 2 + t4 ** 2 + t5 ** 2 + t6 ** 2) / 3 + (0.4D1 / 0.3D1) * t8 * ((d
     #ilogc(-t1 * t2 + 1) + dilogc(-t2 / tman + 1) + dilogc(-t2 * t3 + 1
     #)) * cA - LogMuMtop * tr) - 2 * t8 * (-Kg + cF) + (0.11D2 / 0.9D1)
     # * t10 * t7 + (0.2D1 / 0.9D1) * t9 * cA

