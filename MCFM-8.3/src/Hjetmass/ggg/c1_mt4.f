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
      t8 = t7 ** 2
      t9 = tman * t7
      t10 = uman ** 2
      t11 = t10 ** 2
      t12 = uman * t10
      t13 = sman ** 2
      t14 = t13 ** 2
      t15 = sman * t13
      t16 = t14 + t8 + t11
      t17 = (sman + tman) * t10
      t18 = t13 + t7
      t19 = sman * tman
      t20 = (t9 + t15 + t17) * uman + t19 * t18
      t21 = t13 * t7
      t18 = t10 * t18 + t21
      t22 = -t7 - t13
      t23 = t10 * t22 - t21
      t7 = t13 + t7
      t24 = ((sman + tman) * t10 + t15 + t9) * uman + t19 * t7
      t7 = t10 * t7 + t21
      t5 = t5 ** 2
      t1 = t1 ** 2
      t6 = t6 ** 2
      t10 = t5 + t1 + t6
      t21 = t19 * uman
      t25 = t14 + t8 + t11
      t25 = LogMums * t25 + LogMumt * t25 + LogMumu * t25
      t26 = pi ** 2
      t27 = -t6 - t5 + t26 - t1
      t15 = (-t9 - t15 - t17) * uman + t19 * t22
      t17 = t16 * LogMuMtop
      t19 = t3 + t4 + t2
      t28 = -t3 - t4 - t2
      t9 = cA * (-sman * t12 * t1 - tman * t12 * t1 - t9 * uman * t1 - t
     #11 * t19 + t14 * t28 + t28 * t8) + cA * (t15 * t5 + t15 * t6 + t20
     # * t26 + sman * (-t13 * uman + t22 * tman) * t1) + t17 * tr
      t12 = t21 * cA
      t1 = Kg * t16 + (t1 * t18 + t18 * t5 + t18 * t6 + t23 * t26) * cA
      t5 = LogMums * t24 + LogMumt * t24 + LogMumu * t24
      t6 = LogMums * t7 + LogMumt * t7 + LogMumu * t7
      t13 = cA * (t15 * t2 + t15 * t3 + t15 * t4) + LogMuMtop * t20 * tr
      t15 = sman + tman + uman
      t22 = t12 * t15
      t28 = t21 * LogMuMtop * t15
      t29 = t21 * tr
      t30 = LogMums * t15 + LogMumt * t15 + LogMumu * t15
      t31 = (0.1D1 / 0.21D2)
      t32 = (0.1D1 / 0.42D2)
      t1 = (0.1108D4 / 0.14175D5) * t16 * cF + t31 * (Kg * t20 + (t18 * 
     #t2 + t18 * t3 + t18 * t4) * cA + t23 * LogMuMtop * tr) + LogMuMtop
     # * t24 / 63 + (0.113D3 / 0.1260D4) * cA * t24 - (0.65D2 / 0.378D3)
     # * tr * t24 + (0.73D2 / 0.540D3) * cA * t7 + Kg * t7 / 14 + t32 * 
     #(t18 * LogMuMtop + t21 * (sman * t10 + t10 * tman + uman * t10) * 
     #cA) - (0.50D2 / 0.189D3) * tr * t7 + (0.11D2 / 0.1512D4) * cA * t2
     #5 - cA * (t11 * t27 + t14 * t27 + t27 * t8) / 504 - (0.5D1 / 0.378
     #D3) * tr * t25 - t9 / 126 + (0.2D1 / 0.21D2) * t12 * (sman * t19 +
     # t19 * tman + t19 * uman) + (0.167D3 / 0.7560D4) * t16 * cA + t1 /
     # 84 + t17 / 252 - (0.5D1 / 0.126D3) * t16 * tr + (0.11D2 / 0.378D3
     #) * cA * t5 + (0.11D2 / 0.252D3) * cA * t6 - (0.10D2 / 0.189D3) * 
     #tr * t5 - (0.2D1 / 0.63D2) * t13 - (0.5D1 / 0.63D2) * tr * t6 + (0
     #.31607D5 / 0.117600D6) * t22
      c1mt4 = (0.39D2 / 0.280D3) * t21 * Kg * t15 + (0.33731D5 / 0.37800
     #D5) * t21 * cF * t15 + (0.4432D4 / 0.14175D5) * cF * t24 + (0.2216
     #D4 / 0.4725D4) * cF * t7 + (0.7D1 / 0.80D2) * t12 * t30 - (0.19D2 
     #/ 0.840D3) * t22 * t26 - t22 * LogMuMtop / 140 - (0.13D2 / 0.140D3
     #) * t28 * tr - (0.13D2 / 0.84D2) * t29 * t30 - (0.11623D5 / 0.1890
     #0D5) * t29 * t15 + (0.13D2 / 0.280D3) * t28 + t1

