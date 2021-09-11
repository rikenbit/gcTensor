#
# For testthat/test_Consistency.R
#

# define `expect_all_*()` because `expect_*()` cannot check if multiple
# values are all the same.
expect_all_identical <- function(object, ...) {
    expects <- list(...)
    for (expect in expects) {
        expect_identical(object, expect)
    }
}
expect_all_equal <- function(object, ...) {
    expects <- list(...)
    for (expect in expects) {
        expect_equal(object, expect)
    }
}
expect_all_equivalent <- function(object, ...) {
    expects <- list(...)
    for (expect in expects) {
        expect_equivalent(object, expect)
    }
}

#
# For All
#

# Simulation Datasets
set.seed(123)

# I × J × K
X1 <- rand_tensor(modes = c(4, 5, 6))
X1 <- X1@data^2
names(dim(X1)) <- c("I", "J", "K")

# I × P
X2 <- matrix(runif(4 * 7), nrow=4, ncol=7)
names(dim(X2)) <- c("I", "M")

# J × Q
X3 <- matrix(runif(5 * 8), nrow=5, ncol=8)
names(dim(X3)) <- c("J", "N")

# Coupled Tensor/Matrix
X <- list(X1 = X1, X2 = X2, X3 = X3)

# Coupling matrix R（CP）
R_CP <- rbind(
    c(1,1,1,0,0),
    c(1,0,0,1,0),
    c(0,1,0,0,1)
)
rownames(R_CP) <- paste0("X", seq(3))
colnames(R_CP) <- LETTERS[seq(5)]

# Size of Factor matrices (CP)
Ranks_CP <- list(
    A=list(I=4, r=3),
    B=list(J=5, r=3),
    C=list(K=6, r=3),
    D=list(M=7, r=3),
    E=list(N=8, r=3))

# Coupling matrix R（Tucker）
R_Tucker <- rbind(
    c(1,1,1,1,0,0),
    c(1,0,0,0,1,0),
    c(0,1,0,0,0,1)
)
rownames(R_Tucker) <- paste0("X", seq(3))
colnames(R_Tucker) <- LETTERS[seq(6)]

# Size of Factor matrices (Tucker)
Ranks_Tucker <- list(
    A=list(I=4, p=3),
    B=list(J=5, q=4),
    C=list(K=6, r=3),
    D=list(p=3, q=4, r=3),
    E=list(M=7, p=3),
    F=list(N=8, q=4))

# CP
out.CP_EUC <- GCTF(X, R_CP, Ranks=Ranks_CP, Beta=0, verbose=FALSE)
out.CP_KL <- GCTF(X, R_CP, Ranks=Ranks_CP, Beta=1, verbose=FALSE)
out.CP_IS <- GCTF(X, R_CP, Ranks=Ranks_CP, Beta=2, verbose=FALSE)

# Tucker
out.Tucker_EUC <- GCTF(X, R_Tucker, Ranks=Ranks_Tucker, Beta=0, verbose=FALSE)
out.Tucker_KL <- GCTF(X, R_Tucker, Ranks=Ranks_Tucker, Beta=1, verbose=FALSE)
out.Tucker_IS <- GCTF(X, R_Tucker, Ranks=Ranks_Tucker, Beta=2, verbose=FALSE)

# CP
expect_identical(is.list(out.CP_EUC), TRUE)
expect_identical(is.list(out.CP_KL), TRUE)
expect_identical(is.list(out.CP_IS), TRUE)

# Tucker
expect_identical(is.list(out.Tucker_EUC), TRUE)
expect_identical(is.list(out.Tucker_KL), TRUE)
expect_identical(is.list(out.Tucker_IS), TRUE)

#
# For testthat/test_DecreaseError.R
#
selectThreePoints <- function(x){
    x1 <- x[1]
    x2 <- x[ceiling(median(seq_along(x)))]
    x3 <- rev(x)[1]
    c(x1, x2, x3)
}

#
# For testthat/test_{Estimates,FixZ}.R
#

# The test in Yilmaz and Umut (2011) is performed using the synthetic data.
# No mentions about the scale in the article.
# This implementation compares the true and estimated Zs after normalization.
set.seed(123)

I_test = 30
J_test = 30
K_test = 30
M_test = 30
N_test = 30
r_test = 5
A_test = matrix(runif(I_test*r_test), I_test, r_test); names(dim(A_test)) <- c("I", "r")
B_test = matrix(runif(J_test*r_test), J_test, r_test); names(dim(B_test)) <- c("J", "r")
C_test = matrix(runif(K_test*r_test), K_test, r_test); names(dim(C_test)) <- c("K", "r")
D_test = matrix(runif(M_test*r_test), M_test, r_test); names(dim(D_test)) <- c("M", "r")
E_test = matrix(runif(N_test*r_test), N_test, r_test); names(dim(E_test)) <- c("N", "r")
Z_true <- list(A=A_test, B=B_test, C=C_test, D=D_test, E=E_test)
X1_test <- array(rep(0, I_test*J_test*K_test), dim=c(I_test, J_test, K_test)); names(dim(X1_test)) <- c("I", "J", "K")
for (i in 1:I_test) {
    for (j in 1:J_test) {
        for (k in 1:K_test) {
            X1_test[i, j, k] <- sum(A_test[i, ] * B_test[j, ] * C_test[k, ])
        }
    }
}
X2_test <- array(rep(0, I_test*M_test), dim=c(I_test, M_test)); names(dim(X2_test)) <- c("I", "M")
for (i in 1:I_test) {
    for (m in 1:M_test) {
        X2_test[i, m] <- sum(A_test[i, ] * D_test[m, ])
    }
}
X3_test <- array(rep(0, J_test*N_test), dim=c(J_test, N_test)); names(dim(X3_test)) <- c("J", "N")
for (j in 1:J_test) {
    for (n in 1:N_test) {
        X3_test[j, n] <- sum(B_test[j, ] * E_test[n, ])
    }
}
X_test <- list(X1=X1_test, X2=X2_test, X3=X3_test)

# Visible Indices e.g. I, J, K, P, Q
# visibleIdxs_a <- .visibleIdxs(X_test)
visibleIdxs_a <- gcTensor:::.visibleIdxs(X_test)

# Latent Indices e.g. p, q, r
# latentIdxs_a <- .latentIdxs(Ranks_CP, Z_true, visibleIdxs_a)
latentIdxs_a <- gcTensor:::.latentIdxs(Ranks_CP, Z_true, visibleIdxs_a)

Ranks_test <- list(
    A=list(I=I_test, r=r_test),
    B=list(J=J_test, r=r_test),
    C=list(K=K_test, r=r_test),
    D=list(M=M_test, r=r_test),
    E=list(N=N_test, r=r_test))

#
# For testthat/test_Estimates.R
#
SolveOrthgonalProcrustes <- function(refrence, estimates) {
    # To find the correct permutation, for each of Z_alpha
    # the matching permutation between the original and estimate found
    # by solving an orthgonal Procrustes problem.
    svdResult <- svd(t(refrence) %*% estimates)
    rotation <- svdResult$u %*% t(svdResult$v)
    return(rotation)
}

RotateEstimates <- function(estimates, rotation) {
    return(estimates %*% t(rotation))
}

# These methods are used for Test E-1: CP/MF/MF
# GCTF doesn't use these methods.

.MulTensors <- function(tensorList) {
    dims <- unlist(lapply(tensorList, function(x){dim(x)}))
    dimnames <- sort(unique(names(dims)))
    dims <- dims[dimnames]
    tensorListExpand <- lapply(tensorList, function(x){
        this_dims <- dim(x)
        this_dimnames <- names(this_dims)
        missing_dimnames <- setdiff(dimnames, this_dimnames)
        missing_dims <- dims[missing_dimnames]
        new_x <- array(rep(x, prod(missing_dims)), c(this_dims, missing_dims))
        new_x <- aperm(new_x, order(names(dim(new_x))))

        return(new_x)
    })
    tensorsProd <- 1
    for (tensorExpand in tensorListExpand) {
        tensorsProd <- tensorsProd * tensorExpand
    }
    tensorsProd
}

# Normalization
# compute normalized Z and scaling tensors Lambdas
.normalizeFactors <- function(Z, R, latentIdxs) {

    # normalize all factors.
    # compute normalized factors and scaling tensors.

    # normalized factors:
    #   list. number of elements are same of Z.
    #   each elements are tensors.
    #   the shape of each elements are same of Zs.
    #
    # scaling tensors:
    #   list. number of eleme-nts are same of Z.
    #   each elements are tensors.
    #   the length of each elements are
    #   number of latent indice of each Zs elements.
    #   (for example, scaling tensor S1's shape is [p,q] if Z1 is shape[I,J,p,q])

    normalizedFactors <- list()
    scalingTensors <- list()

    # loop for each factors
    for (Alpha in seq_len(ncol(R))) {

        latentDimNamesOfZAlpha <- intersect(names(latentIdxs), names(dim(Z[[Alpha]])))
        latentIdxLocations <- which(names(dim(Z[[Alpha]])) %in% latentDimNamesOfZAlpha)
        latentDimNamesOfZAlpha <- names(dim(Z[[Alpha]]))[latentIdxLocations]

        if (length(latentDimNamesOfZAlpha) == 1) {
            scalingTensor <- as.array(apply(
                    Z[[Alpha]]^2,
                    which(names(dim(Z[[Alpha]])) %in% latentDimNamesOfZAlpha),
                sum)^(1/2))
            names(dim(scalingTensor)) <- latentDimNamesOfZAlpha
            normalizedFactor <- .MulTensors(list(Z[[Alpha]], 1 / scalingTensor))
        } else {
            scalingTensor <- sqrt(sum(Z[[Alpha]]^2))
            normalizedFactor <- Z[[Alpha]] / scalingTensor
        }

        normalizedFactors[[length(normalizedFactors) + 1]] <- normalizedFactor
        scalingTensors[[length(scalingTensors) + 1]] <- scalingTensor
    }
    # aggregate scaling tensors of factors into scaling tensors of observational tensors.

    # scaling tensors of observational tensors
    #  ... list. number of elements are same of X.
    #      each elements are tensors.
    #      the dimensions of each elements(=tensor) are
    #      all latent indeice of its related factors.

    observationalScalingTensors <- list()

    for (v in seq_len(nrow(R))) {
        scalingTensorsOfRelatedFactors <- list()
        for (Alpha in seq_len(ncol(R))) {
            if (R[v, Alpha] == 1) {
                scalingTensorsOfRelatedFactors[[length(scalingTensorsOfRelatedFactors) + 1]] <- scalingTensors[[Alpha]]
            }
        }
        observationalScalingTensors[[length(observationalScalingTensors) + 1]] <-
            .MulTensors(scalingTensorsOfRelatedFactors)
    }

    result <- list(
        normalizedFactors = normalizedFactors,
        observationalScalingTensors = observationalScalingTensors)

    names(result$normalizedFactors) <- names(Z)
    names(result$observationalScalingTensors) <- names(X)

    result
}
