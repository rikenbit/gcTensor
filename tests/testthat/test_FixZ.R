#
# Test Fixation of Z
#

# Test Z-1
fixAlphas <- c(TRUE, FALSE, TRUE, FALSE, FALSE)

resultGCTF <- GCTF(X=X_test, R=R_CP, Ranks=Ranks_test, Beta=1, num.iter=1)
fixedA <- resultGCTF$Z$A
fixedC <- resultGCTF$Z$C
initial_Z <- resultGCTF$Z
initial_Z$A <- initial_Z$A * 0 + 0.1
initial_Z$C <- initial_Z$C * 0 + 0.2
resultGCTF <- GCTF(X=X_test, R=R_CP, initZ=initial_Z, fix=fixAlphas, Ranks=Ranks_test, Beta=1, num.iter=10)

expect_true(all(resultGCTF$Z$A == initial_Z$A))
expect_true(all(resultGCTF$Z$B != initial_Z$B))
expect_true(all(resultGCTF$Z$C == initial_Z$C))
expect_true(all(resultGCTF$Z$D != initial_Z$D))
expect_true(all(resultGCTF$Z$E != initial_Z$E))
