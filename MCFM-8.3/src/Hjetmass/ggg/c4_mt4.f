      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = dilogc(-t1 * t2 + 1)
      t4 = dilogc(-t2 / tman + 1)
      t5 = 0.1D1 / uman
      t2 = dilogc(-t2 * t5 + 1)
      t6 = cdlogwrap(uman * t1)
      t5 = cdlogwrap(tman * t5)
      t1 = cdlogwrap(tman * t1)
      t7 = t3 + t4 + t2
      t8 = tman ** 2
      t9 = uman ** 2
      t10 = uman * t9
      t11 = sman ** 2
      t12 = t11 + t8 + t9
      t13 = sman + tman
      t14 = sman * tman
      t15 = t13 * uman
      t16 = t15 * LogMuMtop
      t17 = t12 * LogMuMtop
      t18 = t15 + t14
      t19 = pi ** 2
      t1 = t1 ** 2
      t5 = t5 ** 2
      t6 = t6 ** 2
      t20 = t6 - t19 + t1 + t5
      t21 = t9 * cA
      t22 = t9 * tr
      t23 = sman + tman
      t24 = cA * t10
      t25 = t12 * LogMumu
      t26 = t12 * LogMumt
      t27 = t11 + t8
      t27 = tr * ((t27 * t9 + t9 ** 2) * LogMums + t25 * t9 + t26 * t9 +
     # t14 * t27)
      t28 = t14 * cA
      t13 = t10 * t13
      t13 = tr * (t13 * LogMums + t13 * LogMumt + t13 * LogMumu + t11 * 
     #t8)
      t29 = LogMums + LogMumt + LogMumu
      t30 = t21 * t14 * (LogMums + LogMumt)
      t31 = t14 * ((t9 * ((0.31D2 / 0.2160D4) * LogMumu - (0.13D2 / 0.25
     #20D4) * t19 + (0.28429D5 / 0.635040D6)) + t14 / 378) * cA + t9 * (
     #(0.23D2 / 0.840D3) * Kg + (0.3721D4 / 0.22680D5) * cF - (0.1411D4 
     #/ 0.9450D4) * tr) + t9 * (-(0.23D2 / 0.1260D4) * tr + (0.23D2 / 0.
     #2520D4)) * LogMuMtop)
      c4mt4 = (0.2216D4 / 0.14175D5) * cF * t10 * t23 + t9 * ((t11 * t7 
     #+ t14 * LogMuMtop + t7 * t8 + t7 * t9) * cA + t16 - t17 * tr) / 12
     #6 + t9 * (cA * (t18 * t2 + t18 * t3 + t18 * t4) - t16 * tr) / 63 +
     # t21 * (t11 * t20 + t20 * t8 + t20 * t9) / 504 + t9 * ((t1 * t18 -
     # t15 * t19 + t18 * t5 + t18 * t6) * cA + t17) / 252 + (0.1108D4 / 
     #0.14175D5) * t9 * cF * t12 + (0.167D3 / 0.7560D4) * t21 * t12 + t9
     # * Kg * t12 / 84 - (0.5D1 / 0.126D3) * t22 * t12 + (0.167D3 / 0.37
     #80D4) * t24 * t23 + Kg * t10 * t23 / 42 - (0.5D1 / 0.63D2) * tr * 
     #t10 * t23 + (0.11D2 / 0.1512D4) * t21 * (LogMums * t12 + t25 + t26
     #) - (0.5D1 / 0.378D3) * t27 + t28 * (t11 + t8) / 756 - (0.5D1 / 0.
     #189D3) * t13 + (0.43D2 / 0.15120D5) * t28 * uman * t23 - (0.1751D4
     # / 0.37800D5) * t14 * tr * uman * t23 + (0.11D2 / 0.756D3) * t24 *
     # (sman * t29 + t29 * tman) + (0.211D3 / 0.15120D5) * t30 - (0.23D2
     # / 0.756D3) * t22 * t14 * t29 + t31

