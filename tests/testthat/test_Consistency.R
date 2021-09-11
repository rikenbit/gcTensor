#
# Test Consistency of Objects
#

# Test C-1: X and Z
# CP
## I
expect_all_identical(
    dim(X$X1)[1],
    dim(X$X2)[1],
    dim(out.CP_EUC$Z$A)[1],
    dim(out.CP_KL$Z$A)[1],
    dim(out.CP_IS$Z$A)[1])
## J
expect_all_identical(
    dim(X$X1)[2],
    dim(X$X3)[1],
    dim(out.CP_EUC$Z$B)[1],
    dim(out.CP_KL$Z$B)[1],
    dim(out.CP_IS$Z$B)[1])
## K
expect_all_identical(
    dim(X$X1)[3],
    dim(out.CP_EUC$Z$C)[1],
    dim(out.CP_KL$Z$C)[1],
    dim(out.CP_IS$Z$C)[1])
## M
expect_all_identical(
    dim(X$X2)[2],
    dim(out.CP_EUC$Z$D)[1],
    dim(out.CP_KL$Z$D)[1],
    dim(out.CP_IS$Z$D)[1])
## N
expect_all_identical(
    dim(X$X3)[2],
    dim(out.CP_EUC$Z$E)[1],
    dim(out.CP_KL$Z$E)[1],
    dim(out.CP_IS$Z$E)[1])

# Tucker
## I
expect_all_identical(
    dim(X$X1)[1],
    dim(X$X2)[1],
    dim(out.Tucker_EUC$Z$A)[1],
    dim(out.Tucker_KL$Z$A)[1],
    dim(out.Tucker_IS$Z$A)[1])
## J
expect_all_identical(
    dim(X$X1)[2],
    dim(X$X3)[1],
    dim(out.Tucker_EUC$Z$B)[1],
    dim(out.Tucker_KL$Z$B)[1],
    dim(out.Tucker_IS$Z$B)[1])
## K
expect_all_identical(
    dim(X$X1)[3],
    dim(out.Tucker_EUC$Z$C)[1],
    dim(out.Tucker_KL$Z$C)[1],
    dim(out.Tucker_IS$Z$C)[1])
## M
expect_all_identical(
    dim(X$X2)[2],
    dim(out.Tucker_EUC$Z$E)[1],
    dim(out.Tucker_KL$Z$E)[1],
    dim(out.Tucker_IS$Z$E)[1])
## N
expect_all_identical(
    dim(X$X3)[2],
    dim(out.Tucker_EUC$Z$F)[1],
    dim(out.Tucker_KL$Z$F)[1],
    dim(out.Tucker_IS$Z$F)[1])

# Test C-2: X and M
# Expect no error when partially masked some rows
M <- X
M$X1 <- M$X1 * 0 + 1
# Mask only rows 1:3
M$X2 <- M$X2 * 1
M$X2[1:3, 1:7] <- 0
M$X3 <- M$X3 * 0 + 1
expect_error(
    GCTF(X=X, R=R_CP, M=M, Ranks=Ranks_CP, Beta=0),
    regexp=NA)  # regexp=NA means "there should be no errors"

# Expect error when completely masked some columns
M <- X
M$X1 <- M$X1 * 0 + 1
# Mask some column
M$X2 <- M$X2 * 1
M$X2[1:4, 7] <- 0
M$X3 <- M$X3 * 0 + 1
expect_error(GCTF(X=X, R=R_CP, M=M, Ranks=Ranks_CP, Beta=0))

# Expect error when completely masked
M <- X
M$X1 <- M$X1 * 0 + 1
# Mask all rows (1:4)
M$X2 <- M$X2 * 1
M$X2[1:4, 1:7] <- 0
M$X3 <- M$X3 * 0 + 1
expect_error(GCTF(X=X, R=R_CP, M=M, Ranks=Ranks_CP, Beta=0))

# Test C-3: Rank and Z
# CP
## A
expect_all_equal(
    unlist(Ranks_CP$A),
    dim(out.CP_EUC$Z$A),
    dim(out.CP_KL$Z$A),
    dim(out.CP_IS$Z$A))
## B
expect_all_equal(
    unlist(Ranks_CP$B),
    dim(out.CP_EUC$Z$B),
    dim(out.CP_KL$Z$B),
    dim(out.CP_IS$Z$B))
## C
expect_all_equal(
    unlist(Ranks_CP$C),
    dim(out.CP_EUC$Z$C),
    dim(out.CP_KL$Z$C),
    dim(out.CP_IS$Z$C))
## D
expect_all_equal(
    unlist(Ranks_CP$D),
    dim(out.CP_EUC$Z$D),
    dim(out.CP_KL$Z$D),
    dim(out.CP_IS$Z$D))
## E
expect_all_equal(
    unlist(Ranks_CP$E),
    dim(out.CP_EUC$Z$E),
    dim(out.CP_KL$Z$E),
    dim(out.CP_IS$Z$E))

# Tucker
## A
expect_all_equal(
    unlist(Ranks_Tucker$A),
    dim(out.Tucker_EUC$Z$A),
    dim(out.Tucker_KL$Z$A),
    dim(out.Tucker_IS$Z$A))
## B
expect_all_equal(
    unlist(Ranks_Tucker$B),
    dim(out.Tucker_EUC$Z$B),
    dim(out.Tucker_KL$Z$B),
    dim(out.Tucker_IS$Z$B))
## C
expect_all_equal(
    unlist(Ranks_Tucker$C),
    dim(out.Tucker_EUC$Z$C),
    dim(out.Tucker_KL$Z$C),
    dim(out.Tucker_IS$Z$C))
## E
expect_all_equal(
    unlist(Ranks_Tucker$E),
    dim(out.Tucker_EUC$Z$E),
    dim(out.Tucker_KL$Z$E),
    dim(out.Tucker_IS$Z$E))
## F
expect_all_equal(
    unlist(Ranks_Tucker$F),
    dim(out.Tucker_EUC$Z$F),
    dim(out.Tucker_KL$Z$F),
    dim(out.Tucker_IS$Z$F))

# Test C-4: num.iter, RecError, and RelChange
# CP
expect_all_equal(
    formals(GCTF)$num.iter,
    length(out.CP_EUC$RecError)-1,
    length(out.CP_KL$RecError)-1,
    # length(out.CP_IS$RecError)-1,
    length(out.CP_EUC$RelChange)-1,
    length(out.CP_KL$RelChange)-1)
    # length(out.CP_IS$RelChange)-1)

# Tucker
expect_all_equal(
    formals(GCTF)$num.iter,
    length(out.Tucker_EUC$RecError)-1,
    length(out.Tucker_KL$RecError)-1,
    # length(out.Tucker_IS$RecError)-1,
    length(out.Tucker_EUC$RelChange)-1,
    length(out.Tucker_KL$RelChange)-1)
    # length(out.Tucker_IS$RelChange)-1)