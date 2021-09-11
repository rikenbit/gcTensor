#
# Test Output object / type
#

# Test O-1: Object
# CP
expect_identical(is.list(out.CP_EUC), TRUE)
expect_identical(is.list(out.CP_KL), TRUE)
expect_identical(is.list(out.CP_IS), TRUE)

# Tucker
expect_identical(is.list(out.Tucker_EUC), TRUE)
expect_identical(is.list(out.Tucker_KL), TRUE)
expect_identical(is.list(out.Tucker_IS), TRUE)

# Test O-2: Object Names
objNames <- c("X", "Z", "R", "Ranks", "Beta", "num.iter",
	"thr", "RecError", "RelChange")

# CP
expect_identical(names(out.CP_EUC), objNames)
expect_identical(names(out.CP_KL), objNames)
expect_identical(names(out.CP_IS), objNames)

# Tucker
expect_identical(names(out.Tucker_EUC), objNames)
expect_identical(names(out.Tucker_KL), objNames)
expect_identical(names(out.Tucker_IS), objNames)

# Test O-3: X
# CP
expect_identical(out.CP_EUC$X, X)
expect_identical(out.CP_KL$X, X)
expect_identical(out.CP_IS$X, X)

# Tucker
expect_identical(out.Tucker_EUC$X, X)
expect_identical(out.Tucker_KL$X, X)
expect_identical(out.Tucker_IS$X, X)

# Test O-4: Z
# CP
expect_identical(is.list(out.CP_EUC$Z), TRUE)
expect_identical(is.list(out.CP_KL$Z), TRUE)
expect_identical(is.list(out.CP_IS$Z), TRUE)

# Tucker
expect_identical(is.list(out.Tucker_EUC$Z), TRUE)
expect_identical(is.list(out.Tucker_KL$Z), TRUE)
expect_identical(is.list(out.Tucker_IS$Z), TRUE)

# Test O-5: R
# CP
expect_identical(out.CP_EUC$R, R_CP)
expect_identical(out.CP_KL$R, R_CP)
expect_identical(out.CP_IS$R, R_CP)

# Tucker
expect_identical(out.Tucker_EUC$R, R_Tucker)
expect_identical(out.Tucker_KL$R, R_Tucker)
expect_identical(out.Tucker_IS$R, R_Tucker)

# Test O-6: Rank
# CP
expect_identical(out.CP_EUC$Ranks, Ranks_CP)
expect_identical(out.CP_KL$Ranks, Ranks_CP)
expect_identical(out.CP_IS$Ranks, Ranks_CP)

# Tucker
expect_identical(out.Tucker_EUC$Ranks, Ranks_Tucker)
expect_identical(out.Tucker_KL$Ranks, Ranks_Tucker)
expect_identical(out.Tucker_IS$Ranks, Ranks_Tucker)

# Test O-7: Beta
# CP
expect_identical(out.CP_EUC$Beta, 0)
expect_identical(out.CP_KL$Beta, 1)
expect_identical(out.CP_IS$Beta, 2)

# Tucker
expect_identical(out.Tucker_EUC$Beta, 0)
expect_identical(out.Tucker_KL$Beta, 1)
expect_identical(out.Tucker_IS$Beta, 2)

# Test O-8: num.iter
# CP
expect_identical(out.CP_EUC$iter, formals(GCTF)$iter)
expect_identical(out.CP_KL$iter, formals(GCTF)$iter)
expect_identical(out.CP_IS$iter, formals(GCTF)$iter)

# Tucker
expect_identical(out.Tucker_EUC$iter, formals(GCTF)$iter)
expect_identical(out.Tucker_KL$iter, formals(GCTF)$iter)
expect_identical(out.Tucker_IS$iter, formals(GCTF)$iter)

# Test O-9: thr
# CP
expect_identical(out.CP_EUC$thr, formals(GCTF)$thr)
expect_identical(out.CP_KL$thr, formals(GCTF)$thr)
expect_identical(out.CP_IS$thr, formals(GCTF)$thr)

# Tucker
expect_identical(out.Tucker_EUC$thr, formals(GCTF)$thr)
expect_identical(out.Tucker_KL$thr, formals(GCTF)$thr)
expect_identical(out.Tucker_IS$thr, formals(GCTF)$thr)

# Test O-10: RecError
# CP
expect_identical(is.vector(out.CP_EUC$RecError), TRUE)
expect_identical(is.vector(out.CP_KL$RecError), TRUE)
expect_identical(is.vector(out.CP_IS$RecError), TRUE)

# Tucker
expect_identical(is.vector(out.Tucker_EUC$RecError), TRUE)
expect_identical(is.vector(out.Tucker_KL$RecError), TRUE)
expect_identical(is.vector(out.Tucker_IS$RecError), TRUE)

# Test O-11: RelChange
# CP
expect_identical(is.vector(out.CP_EUC$RelChange), TRUE)
expect_identical(is.vector(out.CP_KL$RelChange), TRUE)
expect_identical(is.vector(out.CP_IS$RelChange), TRUE)

# Tucker
expect_identical(is.vector(out.Tucker_EUC$RelChange), TRUE)
expect_identical(is.vector(out.Tucker_KL$RelChange), TRUE)
expect_identical(is.vector(out.Tucker_IS$RelChange), TRUE)

