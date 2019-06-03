      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = dilogc(-t1 * t2 + 1)
      t4 = dilogc(-t2 / tman + 1)
      t5 = 0.1D1 / uman
      t2 = dilogc(-t2 * t5 + 1)
      t6 = cdlogwrap(uman * t1)
      t5 = cdlogwrap(tman * t5)
      t1 = cdlogwrap(tman * t1)
      t7 = tman ** 2
      t8 = tman * t7
      t9 = uman ** 2
      t10 = sman ** 2
      t11 = t10 + t7 + t9
      t12 = sman + uman
      t13 = cA * t8
      t14 = t3 + t4 + t2
      t15 = sman + uman
      t16 = sman * uman
      t17 = t11 * LogMuMtop
      t18 = tman * t15
      t19 = t18 * LogMuMtop
      t20 = t10 + t9
      t21 = t11 * LogMumu
      t22 = t11 * LogMumt
      t23 = t7 * cA
      t24 = (sman + tman) * uman + sman * tman
      t25 = pi ** 2
      t6 = t6 ** 2
      t1 = t1 ** 2
      t5 = t5 ** 2
      t15 = t8 * t15
      t26 = t10 * t9
      t27 = tr * t7
      t28 = t16 * cA
      t29 = -t5 + t25 - t6 - t1
      t29 = t23 * (t10 * t29 + t29 * t7 + t29 * t9)
      t2 = t7 * (cA * (t2 * t24 + t24 * t3 + t24 * t4) - t19 * tr)
      t3 = LogMums + LogMumt + LogMumu
      t4 = t23 * t16 * (LogMums + LogMumu)
      t30 = t16 * t7 * ((0.23D2 / 0.840D3) * Kg + ((0.31D2 / 0.2160D4) *
     # LogMumt - (0.13D2 / 0.2520D4) * t25 + (0.28429D5 / 0.635040D6)) *
     # cA + (0.3721D4 / 0.22680D5) * cF + (-(0.23D2 / 0.1260D4) * tr + (
     #0.23D2 / 0.2520D4)) * LogMuMtop - (0.1411D4 / 0.9450D4) * tr)
      c3mt4 = (0.2216D4 / 0.14175D5) * cF * t8 * t12 - (0.23D2 / 0.756D3
     #) * t27 * t16 * t3 + t26 * cA / 378 + (0.1108D4 / 0.14175D5) * t7 
     #* cF * t11 + (0.167D3 / 0.7560D4) * t23 * t11 + (0.211D3 / 0.15120
     #D5) * t4 + (0.43D2 / 0.15120D5) * t28 * tman * t12 + (0.167D3 / 0.
     #3780D4) * t13 * t12 + (0.11D2 / 0.1512D4) * t23 * (LogMums * t11 +
     # t21 + t22) + (0.11D2 / 0.756D3) * t13 * (sman * t3 + t3 * uman) +
     # t28 * (t10 + t9) / 756 + t7 * ((t1 * t24 - t18 * t25 + t24 * t5 +
     # t24 * t6) * cA + t17) / 252 + t7 * Kg * t11 / 84 + t7 * ((t10 * t
     #14 + t14 * t7 + t14 * t9 + t16 * LogMuMtop) * cA - t17 * tr + t19)
     # / 126 + t2 / 63 + Kg * t8 * t12 / 42 - t29 / 504 - (0.5D1 / 0.378
     #D3) * tr * ((t20 * t7 + t7 ** 2) * LogMums + t21 * t7 + t22 * t7 +
     # t16 * t20) - (0.5D1 / 0.126D3) * t27 * t11 - (0.5D1 / 0.189D3) * 
     #tr * (t15 * LogMums + t15 * LogMumt + t15 * LogMumu + t26) - (0.5D
     #1 / 0.63D2) * t8 * tr * t12 - (0.1751D4 / 0.37800D5) * t16 * tman 
     #* tr * t12 + t30

