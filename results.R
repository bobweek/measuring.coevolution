# bottling the significance results for variable population sizes

#     a given matrix
#
#  n1   50  500 5000  n2
# 50   [  ,    ,    ]
# 500  [  ,    ,    ]
# 5000 [  ,    ,    ]

# coev against null1
resultsc1 <- matrix(numeric(9),nc=3)
resultsc1[1,1] <- 0.000
resultsc1[1,2] <- 0.000
resultsc1[2,1] <- 0.000
resultsc1[1,3] <- 0.128
resultsc1[3,1] <- 0.000
resultsc1[2,2] <- 0.000
resultsc1[2,3] <- 0.000
resultsc1[3,2] <- 0.000
resultsc1[3,3] <- 0.000

# coev against null2
resultsc2 <- matrix(numeric(9),nc=3)
resultsc2[1,1] <- 0.968
resultsc2[1,2] <- 0.964
resultsc2[2,1] <- 0.968
resultsc2[1,3] <- 0.976
resultsc2[3,1] <- 0.916
resultsc2[2,2] <- 0.976
resultsc2[2,3] <- 0.956
resultsc2[3,2] <- 0.984
resultsc2[3,3] <- 0.960

# coev against null3
resultsc3 <- matrix(numeric(9),nc=3)
resultsc3[1,1] <- 1.000
resultsc3[1,2] <- 0.996
resultsc3[2,1] <- 1.000
resultsc3[1,3] <- 0.984
resultsc3[3,1] <- 1.000
resultsc3[2,2] <- 1.000
resultsc3[2,3] <- 1.000
resultsc3[3,2] <- 1.000
resultsc3[3,3] <- 0.996

# null1 agains null3
results13 <- matrix(numeric(9),nc=3)
results13[1,1] <- 1.000
results13[1,2] <- 0.996
results13[2,1] <- 1.000
results13[1,3] <- 0.980
results13[3,1] <- 1.000
results13[2,2] <- 1.000
results13[2,3] <- 1.000
results13[3,2] <- 1.000
results13[3,3] <- 0.996

# null2 against null3
results23 <- matrix(numeric(9),nc=3)
results23[1,1] <- 0.880
results23[1,2] <- 0.776
results23[2,1] <- 0.912
results23[1,3] <- 0.720
results23[3,1] <- 0.932
results23[2,2] <- 0.904
results23[2,3] <- 0.812
results23[3,2] <- 0.984
results23[3,3] <- 0.884


# selection strengths

# n1 x n2:
# 50x50:     B1 = 6.84e-16,  B2 = 0.000164,  A1 = 3.62e-05,  A2 = 0.000128
# 50x500:    B1 = 1.97e-24,  B2 = 4.68e-05,  A1 = 2.36e-05,  A2 = 3.58e-05
# 500x50:    B1 = 1.69e-16,  B2 = 0.000131,  A1 = 3.94e-06,  A2 = 9.97e-05
# 50x5000:   B1 = 6.91e-05,  B2 = 6.67e-07,  A1 = 1.54e-06,  A2 = 8.54e-07
# 5000x50:   B1 = 2.47e-15,  B2 = 0.000126,  A1 = 3.97e-07,  A2 = 9.63e-05
# 500x500:   B1 = 5.64e-24,  B2 = 1.64e-05,  A1 = 3.62e-06,  A2 = 1.28e-05
# 500x5000:  B1 = 9.91e-17,  B2 = 4.68e-06,  A1 = 2.36e-06,  A2 = 3.58e-06
# 5000x500:  B1 = 4.92e-35,  B2 = 1.31e-05,  A1 = 3.94e-07,  A2 = 9.97e-06
# 5000x5000: B1 = 3.49e-26,  B2 = 1.64e-06,  A1 = 3.62e-07,  A2 = 1.28e-06


#     a given matrix
#
#  v1   0.5 1.0  v2
# 0.5  [  ,    ]
# 1.0  [  ,    ]

# coev against null1
Vc1 <- matrix(numeric(9),nc=3)
Vc1[1,1] <- 0.000
Vc1[1,2] <- 0.000
Vc1[2,1] <- 0.000
Vc1[2,2] <- 0

# coev against null2
Vc2 <- matrix(numeric(9),nc=3)
Vc2[1,1] <- 0.956
Vc2[1,2] <- 0.980
Vc2[2,1] <- 0.952
Vc2[2,2] <- 0

# coev against null3
Vc3 <- matrix(numeric(9),nc=3)
Vc3[1,1] <- 1.000
Vc3[1,2] <- 1.000
Vc3[2,1] <- 1.000
Vc3[2,2] <- 1

# null1 agains null3
V13 <- matrix(numeric(9),nc=3)
V13[1,1] <- 1.000
V13[1,2] <- 1.000
V13[2,1] <- 1.000
V13[2,2] <- 1

# null2 agains null3
V23 <- matrix(numeric(9),nc=3)
V23[1,1] <- 0.924
V23[1,2] <- 0.900
V23[2,1] <- 0.832
V23[2,2] <- 0

# selection strengths

# v1 x v2:
# 0.5x0.5:  B1 = 6.32e-29,  B2 = 1.64e-05,  A1 = 3.62e-06,  A2 = 1.28e-05
# 0.5x1.0:  B1 = 2.95e-27,  B2 = 1.47e-05,  A1 = 3.8e-06,   A2 = 1.13e-05
# 1.0x0.5:  B1 = 9.72e-34,  B2 = 1.95e-05,  A1 = 3.32e-06,  A2 = 1.53e-05
# 1.0x1.0:  