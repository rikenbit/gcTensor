#
# Test Monotonous Decrease of Error
#

# Test M-1: RecError
# CP
expect_identical(order(out.CP_EUC$RecError[2:31]), rev(seq(30)))
expect_identical(order(out.CP_KL$RecError[2:31]), rev(seq(30)))

# Not run: In IS case, computation may be unsteady.
# expect_identical(order(out.CP_IS$RecError[2:31]), rev(seq(30)))

# Tucker
expect_identical(order(out.Tucker_EUC$RecError[2:31]), rev(seq(30)))
expect_identical(order(out.Tucker_KL$RecError[2:31]), rev(seq(30)))

# Not run: In IS case, computation may be unsteady.
# expect_identical(order(out.Tucker_IS$RecError[2:31]), rev(seq(30)))

# Test M-2: RelChange
# first element of RelChange doesn't have to be sequentially ordered.

# CP
expect_identical(order(selectThreePoints(out.CP_EUC$RelChange[2:31])), 3:1)
expect_identical(order(selectThreePoints(out.CP_KL$RelChange[2:31])), 3:1)
# expect_identical(order(selectThreePoints(out.CP_IS$RelChange[2:31])), 3:1)

# Tucker
expect_identical(order(selectThreePoints(out.Tucker_EUC$RelChange[2:31])), 3:1)
expect_identical(order(selectThreePoints(out.Tucker_KL$RelChange[2:31])), 3:1)
# expect_identical(order(selectThreePoints(out.Tucker_IS$RelChange[2:31])), 3:1)

# # Check RecError and RelChange of Tucker EUC
# plot(out.Tucker_EUC$RecError)
# plot(out.Tucker_EUC$RelChange, log="y")