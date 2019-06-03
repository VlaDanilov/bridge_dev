
      coeff2(2,2) = (    
     &   2._dp*za(p3,p4)*zb(p1,p3)*zb(p4,p6)/(s123)*(
     &      2._dp*za(p1,p2)*za(p2,p5)*zb(p1,p2)
     &    - za(p1,p5)*za(p2,p4)*zb(p1,p4)
     &    - za(p2,p4)*za(p3,p5)*zb(p3,p4) 
     &    + 2._dp*za(p2,p5)*za(p3,p4)*zb(p3,p4)
     &    )   
     & - 4._dp*zab2(p2,p5,p6,p4)*zab2(p4,p5,p6,p3)/(zab2(p4,p5,p6,p4)*s123)*(
     &     za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p4,p6)
     &    )   
     & - 4._dp*zab2(p2,p5,p6,p4)/(za(p1,p3)*zab2(p4,p5,p6,p4)*s123)*(
     &     za(p1,p2)*za(p1,p5)*za(p2,p4)*za(p3,p4)*zb(p1,p2)**2*zb(p4,p6)
     &    )
     & + zab2(p2,p5,p6,p4)/(zab2(p3,p5,p6,p4))*za(p3,p4)*zb(p4,p6)*(
     &      za(p4,p5)*zb(p1,p4)
     &    - za(p2,p5)*zb(p1,p2)
     &    )   
     & + zab2(p2,p5,p6,p4)/(zab2(p4,p5,p6,p3))*(
     &      za(p4,p5)**2*zb(p1,p3)*zb(p5,p6)
     &    )
     & + 4._dp*zab2(p3,p5,p6,p4)/(za(p1,p3)*zab2(p4,p5,p6,p4)*s123)* (
     &      za(p1,p2)*za(p2,p4)*za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p2,p3)*zb(p4,p6)
     &    )   
     & + zab2(p4,p5,p6,p1)/(zab2(p4,p5,p6,p3)) * (
     &      za(p1,p2)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)
     &    )
     & + zab2(p4,p5,p6,p4)/(zab2(p3,p5,p6,p4))*za(p3,p5)*zb(p1,p4)* (
     &      2._dp*za(p2,p4)*zb(p4,p6)
     &    + za(p2,p5)*zb(p5,p6)
     &    )   
     & - 2._dp*zab2(p4,p5,p6,p4)/(zab2(p4,p5,p6,p3)) * (
     &     za(p2,p4)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)
     &    ) 
     & - 4._dp*za(p1,p2)*za(p3,p4)*za(p5,p6)*zb(p1,p2)*zb(p4,p6)/(za(p1,p3)*zab2(p4,p5,p6,p4)*s123)*(
     &      za(p1,p2)*za(p2,p4)*zb(p1,p2)*zb(p4,p6)
     &    )   
     & - 4._dp*za(p1,p2)*za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p4,p6)*zab2(p4,p1,p3,p4)/(za(p1,p3)*zab2(p4,p5,p6,p4))
     & + 4._dp*za(p1,p2)*za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p4,p6)/(za(p1,p3))
     & - 2._dp*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)/(zab2(p3,p5,p6,p4))*s123    
     & - s123/(zab2(p3,p5,p6,p4)*zab2(p4,p5,p6,p4))*(
     &    za(p2,p4)*za(p3,p4)*za(p5,p6)*zb(p1,p4)*zb(p4,p6)**2
     &    )
     & + za(p2,p4)*za(p5,p6)*zb(p4,p6)**2._dp/(zab2(p3,p5,p6,p4))*(
     &      za(p3,p4)*zb(p1,p4)
     &    - za(p2,p3)*zb(p1,p2)
     &    )   
     & - s123/(zab2(p4,p5,p6,p3))*(
     &    za(p2,p4)*za(p4,p5)*zb(p1,p6)*zb(p3,p4)
     &    )
     & + s123/(zab2(p4,p5,p6,p3)*zab2(p4,p5,p6,p4))*(
     &      za(p2,p4)*za(p4,p5)**2*zb(p1,p4)*zb(p3,p4)*zb(p5,p6)
     &    )   
     & + za(p4,p5)*zb(p5,p6)/(zab2(p4,p5,p6,p3))*(
     &      za(p1,p2)*za(p4,p5)*zb(p1,p3)*zb(p1,p4)
     &    + 2._dp*za(p2,p4)*za(p2,p5)*zb(p1,p2)*zb(p3,p4)
     &    - za(p2,p4)*za(p3,p5)*zb(p1,p3)*zb(p3,p4)
     &    - za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p3,p4)
     &    )   
     & + 2._dp*za(p3,p4)*za(p5,p6)*zb(p1,p3)*zb(p4,p6)/(zab2(p4,p5,p6,p4)*s123)*(
     &      za(p1,p3)*za(p2,p4)*zb(p1,p3)*zb(p4,p6)
     &    - za(p1,p2)*za(p2,p4)*zb(p1,p4)*zb(p2,p6) 
     &    - za(p1,p4)*za(p2,p4)*zb(p1,p4)*zb(p4,p6)
     &    + za(p2,p3)*za(p2,p4)*zb(p2,p6)*zb(p3,p4)    
     &    - 2._dp*za(p1,p2)*za(p3,p4)*zb(p1,p4)*zb(p3,p6)
     &    - 2._dp*za(p1,p4)*za(p2,p3)*zb(p1,p3)*zb(p4,p6)
     &    + 2._dp*za(p2,p3)*za(p3,p4)*zb(p3,p4)*zb(p3,p6)
     &    + za(p2,p4)*za(p3,p4)*zb(p3,p4)*zb(p4,p6)
     &    )   
     & + 4._dp*za(p2,p4)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)/(zab2(p4,p5,p6,p4))*s123 
     & + zb(p4,p6)/(zab2(p4,p5,p6,p4))*(
     &      4._dp*za(p1,p4)*za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p1,p4)
     &    + 2._dp*za(p2,p4)**2*za(p5,p6)*zb(p1,p2)*zb(p4,p6)    
     &    - 2._dp*za(p2,p4)*za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p3,p4)
     &    + 2._dp*za(p2,p4)*za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p2,p3)
     &    )
     & + za(p2,p4)*za(p2,p5)*zb(p1,p2)*zb(p4,p6)
     &    - 2._dp*za(p2,p4)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)   
     &    - 3._dp*za(p2,p5)*za(p4,p5)*zb(p1,p4)*zb(p5,p6)
     &    ) / (4*za(p3,p4)**2*s12*s56)


      coeff2(1,1) = (    
     &    2._dp*za(p4,p5)/(s123)*(
     &      2._dp*za(p1,p2)*za(p1,p3)*zb(p1,p3)*zb(p1,p4)*zb(p1,p6) 
     &    + za(p1,p2)*za(p2,p3)*zb(p1,p3)*zb(p1,p6)*zb(p2,p4)
     &    - 2._dp*za(p1,p2)*za(p2,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,p6)    
     &    - 2._dp*za(p1,p2)*za(p3,p4)*zb(p1,p4)**2*zb(p3,p6)
     &    + 2._dp*za(p1,p3)*za(p2,p3)*zb(p1,p3)*zb(p1,p6)*zb(p3,p4)
     &    - za(p1,p4)*za(p2,p3)*zb(p1,p4)**2*zb(p3,p6)
     &    - za(p2,p3)**2*zb(p1,p2)*zb(p3,p4)*zb(p3,p6)    
     &    + 2._dp*za(p2,p3)**2*zb(p1,p3)*zb(p2,p6)*zb(p3,p4)
     &    - 2._dp*za(p2,p3)*za(p2,p4)*zb(p1,p4)*zb(p2,p3)*zb(p4,p6)
     &    - za(p2,p3)*za(p2,p4)*zb(p1,p6)*zb(p2,p4)*zb(p3,p4)
     &    + za(p2,p3)*za(p2,p5)*zb(p1,p2)*zb(p3,p4)*zb(p5,p6)    
     &    + za(p2,p3)*za(p3,p4)*zb(p1,p4)*zb(p3,p4)*zb(p3,p6)
     &    - 2._dp*za(p2,p3)*za(p3,p4)*zb(p1,p6)*zb(p3,p4)**2
     &    + za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p3,p4)*zb(p5,p6)
     &    )   
     & + zab2(p2,p5,p6,p4)/(zab2(p3,p5,p6,p4))*za(p4,p5)*(
     &      2._dp*za(p2,p3)*zb(p1,p2)*zb(p4,p6)
     &    - za(p1,p3)*zb(p1,p4)*zb(p1,p6)
     &    - za(p2,p3)*zb(p1,p4)*zb(p2,p6)
     &    - 4._dp*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     &    )   
     & - 2._dp*zab2(p4,p5,p6,p1)*za(p1,p2)*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p2,p6)*zb(p3,p4)/(zab2(p4,p5,p6,p4)*s123)
     & + zab2(p4,p5,p6,p3)/(s123)*(
     &      2._dp*za(p1,p5)*za(p2,p3)*zb(p1,p4)*zb(p1,p6)
     &    )   
     & + 2._dp*zab2(p4,p5,p6,p3)/(zab2(p4,p5,p6,p4)*s123)*zb(p1,p4)*zb(p1,p6)*(
     &      za(p1,p3)*za(p2,p3)*za(p4,p5)*zb(p3,p4)
     &    + za(p1,p4)*za(p2,p3)*za(p5,p6)*zb(p4,p6)
     &    - 2._dp*za(p1,p2)*za(p1,p3)*za(p4,p5)*zb(p1,p4)
     &    )   
     & + 2._dp*zab2(p4,p5,p6,p4)/(zab2(p4,p5,p6,p3))*za(p4,p5)*zb(p1,p3)* (
     &      2._dp*za(p2,p4)*zb(p4,p6)
     &    + za(p2,p5)*zb(p5,p6)
     &    )
     & + 4._dp*zab2(p5,p2,p3,p1)/(zab2(p4,p5,p6,p4)*s123)*(
     &      za(p2,p3)*za(p3,p4)*za(p4,p5)*zb(p3,p4)**2*zb(p5,p6)
     &    )   
     & - 8._dp*za(p1,p2)*za(p2,p4)*za(p4,p5)*zb(p1,p2)*zb(p1,p4)*zb(p1,p6)*zb(p3,p4)/(zb(p1,p3)*zab2(p4,p5,p6,p4))
     & + 1._dp/(zab2(p3,p5,p6,p4))*s123*zb(p1,p4)*zb(p4,p6)*(
     &      za(p2,p3)*za(p4,p5) 
     &    + 4._dp*za(p2,p5)*za(p3,p4) 
     &    )   
     & + za(p2,p4)*za(p3,p4)*za(p5,p6)*zb(p1,p4)*zb(p4,p6)**2._dp/(zab2(p3,p5,p6,p4)*zab2(p4,p5,p6,p4))*s123
     & - za(p2,p3)*za(p5,p6)*zb(p4,p6)/(zab2(p3,p5,p6,p4))*(
     &    za(p2,p4)*zb(p1,p2)*zb(p4,p6)
     &    + 4._dp*za(p3,p4)*zb(p1,p4)*zb(p3,p6)
     &    )   
     & + za(p2,p4)*za(p4,p5)*zb(p1,p4)*zb(p3,p6)*s123/(zab2(p4,p5,p6,p3))
     & - za(p2,p4)*za(p4,p5)**2*zb(p1,p4)*zb(p3,p4)*zb(p5,p6)*s123/(zab2(p4,p5,p6,p3)*zab2(p4,p5,p6,p4))
     & - 4._dp*za(p2,p4)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*s123/(zab2(p4,p5,p6,p4))   
     & + 2._dp*za(p4,p5)*zb(p3,p4)*zb(p5,p6)/(zab2(p4,p5,p6,p3))*(
     &      2._dp*za(p2,p3)*za(p4,p5)*zb(p1,p3)
     &    + za(p2,p4)*za(p2,p5)*zb(p1,p2) 
     &    + za(p2,p5)*za(p3,p4)*zb(p1,p3) 
     &    )   
     & + 2._dp*za(p4,p5)*zb(p1,p4)*zb(p5,p6)/(zab2(p4,p5,p6,p4)*s123)*(
     &      2._dp*za(p1,p2)*za(p3,p4)*za(p4,p5)*zb(p1,p4)*zb(p3,p4)
     &    - za(p1,p2)*za(p2,p3)*za(p4,p5)*zb(p1,p2)*zb(p3,p4)
     &    - za(p1,p4)*za(p2,p3)*za(p5,p6)*zb(p1,p3)*zb(p4,p6)
     &    )     
     & + 2._dp*za(p4,p5)/(zab2(p4,p5,p6,p4))*(
     &      2._dp*za(p2,p3)*za(p2,p4)*zb(p1,p2)*zb(p3,p4)*zb(p4,p6)
     &    - 2._dp*za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p3,p4)*zb(p5,p6)     
     &    - za(p1,p2)*za(p4,p5)*zb(p1,p4)**2*zb(p5,p6)
     &    - za(p1,p4)*za(p2,p3)*zb(p1,p3)*zb(p1,p4)*zb(p4,p6)
     &    - za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p3,p4)*zb(p5,p6)
     &    )   
     & + za(p2,p3)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)
     &    + 4._dp*za(p2,p3)*za(p4,p5)*zb(p1,p6)*zb(p3,p4)
     &    - za(p2,p4)*za(p2,p5)*zb(p1,p2)*zb(p4,p6)
     &    - 5._dp*za(p2,p4)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)     
     &    - za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p4,p6)
     &    + za(p2,p5)*za(p4,p5)*zb(p1,p4)*zb(p5,p6)
     &    ) / (4*zb(p3,p4)**2*s12*s56)


      coeff2(2,1)  = (    
     &    2._dp*za(p4,p5)/(s123)*(
     &      3._dp*za(p1,p4)*za(p2,p4)*zb(p1,p3)*zb(p1,p6)*zb(p3,p4)
     &    - za(p1,p2)*za(p2,p4)*zb(p1,p3)**2*zb(p2,p6)
     &    - 2._dp*za(p1,p4)*za(p2,p5)*zb(p1,p3)**2*zb(p5,p6)    
     &    - 4._dp*za(p2,p4)**2*zb(p1,p2)*zb(p3,p4)*zb(p3,p6)
     &    + 2._dp*za(p2,p4)**2*zb(p1,p3)*zb(p2,p6)*zb(p3,p4)
     &    - za(p2,p4)*za(p2,p5)*zb(p1,p3)*zb(p2,p3)*zb(p5,p6)
     &    + za(p2,p4)*za(p3,p4)*zb(p1,p3)*zb(p3,p4)*zb(p3,p6)
     &    )   
     & - zab2(p2,p5,p6,p4)*zab2(p4,p5,p6,p3)*za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p4,p6)/(zab2(p3,p5,p6,p4)**2)
     & + zab2(p2,p5,p6,p4)*zab2(p4,p5,p6,p3)*za(p4,p5)*zb(p1,p6)/(zab2(p3,p5,p6,p4))   
     & + zab2(p2,p5,p6,p4)/(zab2(p3,p5,p6,p4))*za(p4,p5)*zb(p3,p6)*(
     &      za(p3,p4)*zb(p1,p3)
     &    - za(p2,p4)*zb(p1,p2)
     &    )
     & - zab2(p4,p5,p6,p1)*zab2(p4,p5,p6,p3)*za(p2,p5)*zb(p4,p6)/(zab2(p3,p5,p6,p4))   
     & + 4._dp*zab2(p4,p5,p6,p1)/(za(p1,p3)*zab2(p4,p5,p6,p4)*s123)*(
     &      za(p1,p2)*za(p1,p4)*za(p2,p5)*za(p4,p5)*zb(p1,p2)*zb(p3,p4)*zb(p5,p6)
     &    )   
     & - zab2(p4,p5,p6,p1)/(zab2(p3,p5,p6,p4)**2)*s123 * (
     &    za(p2,p4)*za(p3,p5)*zb(p3,p4)*zb(p4,p6)
     &    ) 
     & + zab2(p4,p5,p6,p1)/(zab2(p3,p5,p6,p4))*za(p4,p5)*zb(p3,p4)*(
     &      3._dp*za(p2,p4)*zb(p4,p6)
     &    + za(p2,p5)*zb(p5,p6)
     &    )   
     & - 4._dp*zab2(p4,p5,p6,p3)*zab2(p4,p5,p6,p4)/(zab2(p3,p5,p6,p4)**2) * (
     &     za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     &    )
     & - zab2(p4,p5,p6,p3)/(zab2(p3,p5,p6,p4)) * (
     &    za(p2,p3)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)
     &    )   
     & + 4._dp*zab2(p4,p5,p6,p3)/(za(p1,p3)*zab2(p4,p5,p6,p4)*s123)*za(p1,p2)*za(p3,p4)*za(p4,p5)*zb(p1,p2)*zb(p3,p4)*(
     &      za(p1,p2)*zb(p1,p6)
     &    - za(p2,p3)*zb(p3,p6)
     &    )   
     & + 3._dp*zab2(p4,p5,p6,p3)/(zab2(p3,p5,p6,p4)**2)*s123 * (
     &      za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     &    )
     & - zab2(p4,p5,p6,p3)/(zab2(p3,p5,p6,p4)**2)*za(p2,p4)*za(p5,p6)*zb(p4,p6)**2*(
     &    za(p2,p3)*zb(p1,p2) 
     &    + 4._dp*za(p3,p4)*zb(p1,p4)
     &    )   
     & + 2._dp*zab2(p4,p5,p6,p3)*za(p4,p5)/(zab2(p4,p5,p6,p4)*s123)*(
     &       za(p1,p2)*za(p2,p4)*zb(p1,p3)*zb(p1,p4)*zb(p2,p6)
     &    + 2._dp*za(p1,p2)*za(p3,p4)*zb(p1,p3)*zb(p1,p6)*zb(p3,p4)
     &    + 2._dp*za(p1,p4)*za(p2,p5)*zb(p1,p3)*zb(p1,p4)*zb(p5,p6)    
     &    - za(p2,p3)*za(p2,p4)*zb(p1,p3)*zb(p2,p6)*zb(p3,p4)
     &    - 2._dp*za(p2,p3)*za(p3,p4)*zb(p1,p3)*zb(p3,p4)*zb(p3,p6)
     &    - 2._dp*za(p2,p4)*za(p2,p5)*zb(p1,p2)*zb(p3,p4)*zb(p5,p6)
     &    + za(p2,p4)*za(p2,p5)*zb(p1,p3)*zb(p2,p4)*zb(p5,p6)
     &    )   
     & + zab2(p4,p5,p6,p4)/(zab2(p3,p5,p6,p4))*za(p4,p5)*zb(p1,p6)*(
     &      za(p2,p4)*zb(p3,p4)
     &    - za(p1,p2)*zb(p1,p3)
     &    ) 
     & + 4._dp*za(p1,p2)*za(p2,p4)*za(p4,p5)*zb(p1,p2)*zb(p3,p4)/(za(p1,p3)*s123)*(
     &      2._dp*za(p1,p4)*zb(p1,p6)
     &    + za(p3,p4)*zb(p3,p6)
     &    )   
     & + 4._dp*za(p1,p2)*za(p2,p4)*za(p4,p5)*zb(p1,p2)*zb(p3,p4)*zb(p5,p6)/(za(p1,p3)*zab2(p4,p5,p6,p4)*s123)*(
     &      za(p2,p5)*za(p3,p4)*zb(p2,p3)
     &    - za(p3,p4)*za(p4,p5)*zb(p3,p4)
     &    - za(p1,p4)*za(p2,p5)*zb(p1,p2)
     &    - 2._dp*za(p1,p4)*za(p4,p5)*zb(p1,p4)
     &    )    
     & +za(p2,p4)*za(p3,p4)*za(p5,p6)*zb(p1,p3)*zb(p4,p6)**2*s123/(zab2(p3,p5,p6,p4)**2)
     & + za(p2,p4)/(zab2(p3,p5,p6,p4))*(
     &      2._dp*za(p2,p4)*za(p5,p6)*zb(p1,p2)*zb(p3,p6)*zb(p4,p6)    
     &    - za(p2,p5)*za(p4,p5)*zb(p1,p2)*zb(p3,p4)*zb(p5,p6)
     &    - za(p3,p4)*za(p5,p6)*zb(p1,p3)*zb(p3,p6)*zb(p4,p6)
     &    - za(p4,p5)**2*zb(p1,p4)*zb(p3,p4)*zb(p5,p6)
     &    )   
     & + 2._dp*za(p2,p4)*za(p4,p5)*zb(p3,p4)*zb(p5,p6)/(zab2(p4,p5,p6,p4)*s123)*(
     &      4._dp*za(p2,p4)*za(p4,p5)*zb(p1,p2)*zb(p3,p4) 
     &    - 3._dp*za(p1,p4)*za(p4,p5)*zb(p1,p3)*zb(p1,p4)
     &    - 2._dp*za(p1,p4)*za(p2,p5)*zb(p1,p2)*zb(p1,p3)     
     &    - 2._dp*za(p2,p4)*za(p2,p5)*zb(p1,p2)*zb(p2,p3) 
     &    - 2._dp*za(p2,p4)*za(p4,p5)*zb(p1,p3)*zb(p2,p4)
     &    - za(p1,p3)*za(p4,p5)*zb(p1,p3)**2
     &    - za(p3,p4)*za(p4,p5)*zb(p1,p3)*zb(p3,p4)
     &    )   
     & - 2._dp*za(p2,p4)*za(p4,p5)**2*zb(p1,p3)*zb(p3,p4)*zb(p5,p6)/(zab2(p4,p5,p6,p4))
     & + za(p2,p4)*za(p4,p5)*zb(p1,p3)*zb(p3,p6)
     &    ) / (4*za(p3,p4)*zb(p3,p4)*s12*s56)


      coeff2(1,2)  = (    
     &   2._dp*za(p2,p3)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)/(s123)*(
     &      za(p1,p5)*zb(p1,p4)
     &    - za(p3,p5)*zb(p3,p4)
     &    )
     & - 4._dp*zab2(p2,p5,p6,p4)/(zab2(p4,p5,p6,p4)*s123)*(
     &     za(p2,p3)*za(p3,p4)*za(p5,p6)*zb(p1,p2)*zb(p4,p6)**2
     &    )   
     & + 4._dp*zab2(p3,p5,p6,p4)/(zab2(p4,p5,p6,p3)**2)*za(p2,p4)*za(p4,p5)*zb(p3,p4)*zb(p5,p6)*(
     &      za(p2,p5)*zb(p1,p2)
     &    + za(p3,p5)*zb(p1,p3)
     &    )   
     & + zab2(p3,p5,p6,p4)/(zab2(p4,p5,p6,p3))*(
     &      4._dp*za(p2,p3)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)
     &    + 2._dp*za(p2,p4)*za(p2,p5)*zb(p1,p2)*zb(p4,p6)
     &    +   za(p2,p4)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)    
     &    + 2._dp*za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p4,p6)
     &    + 3._dp*za(p2,p5)*za(p4,p5)*zb(p1,p4)*zb(p5,p6)
     &    )   
     & + 4._dp*zab2(p3,p5,p6,p4)*za(p3,p4)*zb(p4,p6) /(zab2(p4,p5,p6,p4)*s123)*(
     &      za(p1,p2)*za(p2,p5)*zb(p1,p2)*zb(p1,p4)
     &    + za(p2,p3)*za(p2,p5)*zb(p1,p2)*zb(p3,p4)
     &    + za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p3,p4)
     &    )   
     & - 4._dp*zab2(p3,p5,p6,p4)/(zab2(p4,p5,p6,p4)) * (
     &     za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     &    )
     & - 8._dp*za(p1,p2)*za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p1,p4)**2*zb(p4,p6)/(zb(p1,p3)*zab2(p4,p5,p6,p4))   
     & - zab2(p4,p5,p6,p4)/(zab2(p4,p5,p6,p3))*za(p2,p5)*(
     &      2._dp*za(p2,p3)*zb(p1,p2)*zb(p4,p6)
     &    + 2._dp*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     &    + za(p3,p5)*zb(p1,p4)*zb(p5,p6)
     &    )   
     & + zb(p1,p4)/(zab2(p4,p5,p6,p3))*(
     &      za(p1,p2)*za(p2,p4)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)
     &    + za(p1,p2)*za(p3,p4)*za(p3,p5)*zb(p1,p3)*zb(p4,p6)
     &    - za(p2,p4)*za(p3,p4)*za(p5,p6)*zb(p4,p6)**2     
     &    - 2._dp*za(p1,p2)*za(p3,p5)*za(p4,p5)*zb(p1,p4)*zb(p5,p6)
     &    - 2._dp*za(p2,p3)*za(p3,p5)*za(p4,p5)*zb(p3,p4)*zb(p5,p6)
     &    + za(p2,p3)*za(p4,p5)*za(p5,p6)*zb(p4,p6)*zb(p5,p6)
     &    )   
     & + 2._dp*za(p3,p4)*zb(p1,p4)*zb(p4,p6)/(zab2(p4,p5,p6,p4)*s123)*(
     &      2._dp*za(p1,p2)*za(p2,p3)*za(p3,p5)*zb(p1,p2)*zb(p3,p4)
     &    - 2._dp*za(p1,p2)**2*za(p3,p5)*zb(p1,p2)*zb(p1,p4)    
     &    + 5._dp*za(p1,p2)*za(p2,p3)*za(p5,p6)*zb(p1,p2)*zb(p4,p6)
     &    + za(p1,p2)*za(p2,p3)*za(p5,p6)*zb(p1,p4)*zb(p2,p6)
     &    + 2._dp*za(p1,p3)*za(p2,p3)*za(p5,p6)*zb(p1,p4)*zb(p3,p6)     
     &    + za(p1,p4)*za(p2,p3)*za(p5,p6)*zb(p1,p4)*zb(p4,p6)
     &    + za(p2,p3)**2*za(p5,p6)*zb(p2,p4)*zb(p3,p6)
     &    - za(p2,p3)*za(p3,p4)*za(p5,p6)*zb(p3,p4)*zb(p4,p6)
     &    )   
     & + 2._dp*za(p3,p4)*zb(p1,p4)*zb(p4,p6)/(zab2(p4,p5,p6,p4))*(
     &      za(p2,p3)*za(p3,p5)*zb(p3,p4)
     &    - za(p1,p2)*za(p3,p5)*zb(p1,p4)
     &    - za(p1,p3)*za(p2,p5)*zb(p1,p4)
     &    - za(p2,p3)*za(p5,p6)*zb(p4,p6)
     &    )   
     & - 3._dp*za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)
     &    ) / (4*za(p3,p4)*zb(p3,p4)*s12*s56)


