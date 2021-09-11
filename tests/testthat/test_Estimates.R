#
# Test Estimates
#
# Test E-1: CP/MF/MF
# generate initial Z
# .initZ <- .genLatentVals(X_test, R_CP, Ranks_test)
.initZ <- gcTensor:::.genLatentVals(X_test, R_CP, Ranks_test)

# normalize initial Z
# visibleIdxs <- .visibleIdxs(X_test)
# latentIdxs <- .latentIdxs(Ranks_test, .initZ, visibleIdxs)
# normalizedInitZ <- .normalizeFactors(.initZ, R_CP, latentIdxs)$normalizedFactors
visibleIdxs <- gcTensor:::.visibleIdxs(X_test)
latentIdxs <- gcTensor:::.latentIdxs(Ranks_test, .initZ, visibleIdxs)
normalizedInitZ <- gcTensor:::.normalizeFactors(.initZ, R_CP, latentIdxs)$normalizedFactors

# estimate Z
resultGCTF <- GCTF(
    X_test, R_CP, M=NULL, initZ=.initZ, fix=NULL,
    Ranks_test, Beta=1, num.iter=1000, thr=1E-10, verbose=TRUE)

# normalize estimated Z
# normalizedEstimatedZ <- .normalizeFactors(resultGCTF$Z, R_CP, latentIdxs)$normalizedFactors
normalizedEstimatedZ <- gcTensor:::.normalizeFactors(resultGCTF$Z, R_CP, latentIdxs)$normalizedFactors

for (factorName in c("A", "B", "C", "D", "E")) {
    # plot answer
    plot(1:10,
        (Z_true[[factorName]] / sqrt(colSums(Z_true[[factorName]]^2)))[1:10, 1],
        type="b", col='black', pch="+", main=factorName, xlab="", ylab="")
    # rotate initial Z
    rotation <- SolveOrthgonalProcrustes(Z_true[[factorName]],
        normalizedInitZ[[factorName]])
    rotatedEstimates <- RotateEstimates(normalizedInitZ[[factorName]], rotation)
    # plot initial Z
    lines(1:10, rotatedEstimates[1:10, 1], type="b", col='blue', pch="*")
    # rotate estimated Z
    rotation <- SolveOrthgonalProcrustes(Z_true[[factorName]], normalizedEstimatedZ[[factorName]])
    rotatedEstimates <- RotateEstimates(normalizedEstimatedZ[[factorName]], rotation)
    # plot estimated Z
    lines(1:10, rotatedEstimates[1:10, 1], type="b", col='red', pch="o")
    if (factorName == 'A') {
        legend("topleft",
               legend = c("Original", "Initial", "Final"),
               col = c("black", "blue", "red"),
               pch = c("+", "*", "o"))
    }
}

# Visualize estimated factors (not normalized)
# to check if underflow or overflow happens.
for (factorName in c("A", "B", "C", "D", "E")) {
    ts.plot(resultGCTF$Z[[factorName]], main=factorName, xlab="row", ylab="", col=1:ncol(resultGCTF$Z[[factorName]]))
    legend("topright", legend=paste("column", 1:ncol(resultGCTF$Z[[factorName]])), col=1:ncol(resultGCTF$Z[[factorName]]), lty=1, cex=1)
}