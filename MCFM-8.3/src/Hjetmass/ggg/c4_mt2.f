      t1 = 0.1D1 / sman
      t2 = cdlogwrap(tman * t1)
      t3 = cdlogwrap(uman * t1)
      t4 = 0.1D1 / uman
      t5 = cdlogwrap(tman * t4)
      t6 = mH ** 2
      t7 = sman + tman + uman
      t8 = uman ** 2
      t9 = cA * t8
      t10 = t7 * LogMumu
      t11 = t7 * LogMumt
      t2 = pi ** 2 - t2 ** 2 - t3 ** 2 - t5 ** 2
      t1 = dilogc(-t1 * t6 + 1) + dilogc(-t6 / tman + 1) + dilogc(-t4 * 
     #t6 + 1)
      t3 = sman + tman
      t4 = tman * sman
      c4mt2 = (0.7D1 / 0.540D3) * cA * sman * tman * t7 + (0.7D1 / 0.60D
     #2) * Kg * t8 * t7 - (0.7D1 / 0.18D2) * tr * t8 * t7 + (0.7D1 / 0.1
     #80D3) * LogMuMtop * t8 * t7 + (0.289D3 / 0.1080D4) * t9 * t7 + (0.
     #77D2 / 0.1080D4) * t9 * (LogMums * t7 + t10 + t11) - (0.7D1 / 0.36
     #0D3) * t9 * (sman * t2 + t2 * tman + t2 * uman) + (0.7D1 / 0.90D2)
     # * t8 * (cA * (sman * t1 + tman * t1 + uman * t1) + (-sman - tman 
     #- uman) * LogMuMtop * tr) - (0.7D1 / 0.54D2) * tr * ((t3 * t8 + um
     #an * t8) * LogMums + t10 * t8 + t11 * t8 + t4 * t3) + (0.61D2 / 0.
     #135D3) * cF * t8 * t7 - (0.14D2 / 0.45D2) * t4 * tr * uman

