#
# Test Input object / type
#

# Test I-1: Object Names
expect_identical(names(formals(GCTF)),
    c("X", "R", "M", "initZ", "fix", "Ranks", "Beta", "num.iter", "thr", "verbose"))

# Test I-2: X
expect_identical(as.character(formals(GCTF)$X), "")

# Test I-3: R
expect_identical(as.character(formals(GCTF)$R), "")

# Test I-4: M
expect_identical(formals(GCTF)$M, NULL)

# Test I-5: initZ
expect_identical(formals(GCTF)$initZ, NULL)

# Test I-6: fix
expect_identical(formals(GCTF)$fix, NULL)

# Test I-7: Ranks
expect_identical(as.character(formals(GCTF)$Ranks), "")

# Test I-8: Beta
expect_identical(formals(GCTF)$Beta, 1)

# Test I-9: num.iter
expect_identical(formals(GCTF)$num.iter, 30)

# Test I-10: thr
expect_identical(formals(GCTF)$thr, 1E-10)

# Test I-11: verbose
expect_identical(formals(GCTF)$verbose, FALSE)
