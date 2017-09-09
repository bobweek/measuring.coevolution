# ACTUAL RESULTS

# UNDER PERFECT HERITABILITY FOR THE WEEVIL

# B1 = 6.4e-11,  B2 = 2.33e-09,  A1 = 1.19e-06,  A2 = 1.05e-05

# P-value against null 1 ~ 0.000
# P-value against null 2 ~ 0.012
# P-value against null 3 ~ 0.000
# P-value null 1 against null 3 ~ 0.008
# P-value null 2 against null 3 ~ 0.000

# UNDER A HERITABILITY OF 0.5 FOR THE WEEVIL

# B1 = 6.4e-11,  B2 = 2.33e-09,  A1 = 1.19e-06,  A2 = 1.05e-05

# P-value against null 1 ~ 0.000
# P-value against null 2 ~ 0.016
# P-value against null 3 ~ 0.000
# P-value null 1 against null 3 ~ 0.008
# P-value null 2 against null 3 ~ 0.000


# OLD OLD OLD \/ \/ \/ \/ \/

# bottling the significance results for variable population sizes

#     a given matrix
#
#  n1    100 1000 10000  n2
# 100   [  ,     ,     ]
# 1000  [  ,     ,     ]
# 10000 [  ,     ,     ]

# coev against null 1
resultsc1 <- matrix(numeric(9),nc=3)
resultsc1[1,1] <- 1.000
resultsc1[1,2] <- 1.000
resultsc1[2,1] <- 1.000
resultsc1[2,2] <- 1.000
resultsc1[1,3] <- 1.000
resultsc1[3,1] <- 1.000
resultsc1[2,3] <- 1.000
resultsc1[3,2] <- 1.000
resultsc1[3,3] <- 1.000

# coev against null 2
resultsc2 <- matrix(numeric(9),nc=3)
resultsc2[1,1] <- 0.972
resultsc2[1,2] <- 0.976
resultsc2[2,1] <- 0.976
resultsc2[2,2] <- 0.960
resultsc2[1,3] <- 0.976
resultsc2[3,1] <- 0.988
resultsc2[2,3] <- 0.992
resultsc2[3,2] <- 0.996
resultsc2[3,3] <- 0.996

# coev against null 3
resultsc3 <- matrix(numeric(9),nc=3)
resultsc3[1,1] <- 1.000
resultsc3[1,2] <- 1.000
resultsc3[2,1] <- 1.000
resultsc3[2,2] <- 1.000
resultsc3[1,3] <- 1.000
resultsc3[3,1] <- 1.000
resultsc3[2,3] <- 1.000
resultsc3[3,2] <- 1.000
resultsc3[3,3] <- 1.000

# null 1 against null 3
results13 <- matrix(numeric(9),nc=3)
results13[1,1] <- 0.980
results13[1,2] <- 0.972
results13[2,1] <- 0.976
results13[2,2] <- 0.984
results13[1,3] <- 0.992
results13[3,1] <- 0.976
results13[2,3] <- 0.984
results13[3,2] <- 0.984
results13[3,3] <- 0.996

# null 2 against null 3
results23 <- matrix(numeric(9),nc=3)
results23[1,1] <- 1.000
results23[1,2] <- 1.000
results23[2,1] <- 1.000
results23[2,2] <- 1.000
results23[1,3] <- 1.000
results23[3,1] <- 1.000
results23[2,3] <- 1.000
results23[3,2] <- 1.000
results23[3,3] <- 1.000


# selection strengths

# n1 x n2:
# 100x100:      B1 = 6.52e-06,  B2 = 7.37e-07,  A1 = 0.000376,  A2 = 0.000187
# 100x1000:     B1 = 7.16e-06,  B2 = 7.37e-09,  A1 = 0.000392,  A2 = 1.87e-05
# 1000x100:     B1 = 6.05e-08,  B2 = 7.61e-07,  A1 = 3.66e-05,  A2 = 0.000189
# 1000x1000:    B1 = 6.1e-08,   B2 = 7.45e-09,  A1 = 3.67e-05,  A2 = 1.88e-05
# 100x10000:    B1 = 7.27e-06,  B2 = 7.43e-11,  A1 = 0.000394,  A2 = 1.88e-06
# 10000x100:    B1 = 6.06e-10,  B2 = 7.62e-07,  A1 = 3.66e-06,  A2 = 0.000189
# 1000x10000:   B1 = 6.15e-08,  B2 = 7.44e-11,  A1 = 3.68e-05,  A2 = 1.88e-06
# 10000x1000:   B1 = 6.06e-10,  B2 = 7.47e-09,  A1 = 3.66e-06,  A2 = 1.88e-05
# 10000x10000:  B1 = 6.06e-10,  B2 = 7.45e-11,  A1 = 3.66e-06,  A2 = 1.88e-06

# REDO 100x100, 300x300 and 1000x1000!!

#     a given matrix
#     using n1xn2 = 300x300
#  v1   
# 0.5  [  ]  
# 1.0  [  ]

# coev against null1
Vc1 <- numeric(2)
Vc1[1] <- 1.000
Vc1[2] <- 0

# coev against null2
Vc2 <- numeric(2)
Vc2[1] <- 0.980
Vc2[2] <- 0

# coev against null3
Vc3 <- numeric(2)
Vc3[1] <- 1.00
Vc3[2] <- 0

# null1 agains null3
V13 <- numeric(2)
V13[1] <- 0.992
V13[2] <- 0

# null2 agains null3
V23 <- numeric(2)
V23[1] <- 1.000
V23[2] <- 0

# selection strengths

# v1:
# 0.5:  B1 = 6.08e-08,  B2 = 7.45e-09,  A1 = 3.67e-05,  A2 = 1.88e-05
# 1.0:  
