#
# Test Non-negative data
#

# Test NN-1: X
NegativeX <- X
NegativeX$X2[2, 3] <- -1
expect_error(GCTF(NegativeX, R_CP, Ranks=Ranks_CP, Beta=0))