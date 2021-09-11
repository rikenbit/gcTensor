GCTF <- function(X, R, M=NULL, initZ=NULL, fix=NULL, Ranks, Beta=1,
    num.iter=30, thr=1E-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkGCTF(X, R, Ranks, Beta, num.iter, thr)
    ######################################
    # Setting
    ######################################
    int <- .initGCTF(X, R, M, initZ, fix, Ranks, thr)
    M <- int$M
    Z <- int$Z
    fix <- int$fix
    visibleIdxs <- int$visibleIdxs
    latentIdxs <- int$latentIdxs
    RecError <- int$RecError
    RelChange <- int$RelChange
    X_bar <- int$X_bar
    ######################################
    # Iteration
    ######################################
    iter <- 1
    while((RelChange[iter] > thr) && (iter <= num.iter)){
      # Update each Factor Matrix
      for(Alpha in seq_len(ncol(R))){
        # Skip update if current Alpha is set fixed.
        if(fix[Alpha]) next
        # Update factor matrix
        Z[[Alpha]] <- .MU(Z, M, X, X_bar, R, Alpha, Beta, Ranks)
        # Reconstruct observational tensor
        X_bar <- .recTensor(X, Z, R, Ranks, latentIdxs)
      }
      # After Update
      iter <- iter + 1
      RecError[iter] <- .recError(X, X_bar)
      RelChange[iter] <- .relChange(iter, RecError)
      # Verbose
      if(verbose){
             cat(paste0(iter - 1, " / ", num.iter,
                " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
      }
    }
    # Output
    list(X=X, Z=Z, R=R, Ranks=Ranks, Beta=Beta,
        num.iter=num.iter, thr=thr,
        RecError=RecError, RelChange=RelChange)
}

# Argument Check Functions
.checkGCTF <- function(X, R, Ranks, Beta, num.iter, thr){
    # Sign
    .checkSign(X)
    # Common Name
    .checkCommonName(X, R, Ranks)
    # Beta
    .checkBeta(Beta)
    # Number of iteration
    .checkNumIteration(num.iter)
    # Threshold value
    .checkthr(thr)
    # Ranks
    .checkRanks(X, R, Ranks)
}

.checkSign <- function(X){
    lapply(X, function(x){
        if(!all(x >= 0)){
            stop("Specify the input X as non-negative")
        }
    })
}

.checkCommonName <- function(X, R, Ranks){
    if(!all(rownames(R) %in% names(X))){
        stop("rownames(R) must be the same as names(X)")
    }
    if(!all(colnames(R) %in% names(Ranks))){
        stop("rownames(R) must be the same as names(X)")
    }
    namesX <- unique(unlist(lapply(X, function(x){names(dim(x))})))
    namesRanks <- unique(unlist(lapply(Ranks, names)))
    if(!all(namesX %in% namesRanks)){
        stop("rownames(R) must be the same as names(X)")
    }
}

.checkBeta <- function(Beta){
    check1 <- is.numeric(Beta)
    check2 <- length(Beta) == 1
    if(!check1 || !check2){
        stop("Specify the Beta as single numerical value")
    }
}

.checkNumIteration <- function(num.iter){
    if(num.iter <= 0){
        stop("Specify num.iter as positive integar (e.g. 100)")
    }
}

.checkthr <- function(thr){
    if(thr <= 0){
        stop("Specify thr as positive real number (e.g. 1E-10)")
    }
}

.checkRanks <- function(X, R, Ranks){
    Dims <- lapply(X, function(x){dim(x)})
    lapply(seq_len(ncol(R)), function(x){
        objects <- rownames(R)[which(R[,x] == 1)]
        if(min(unlist(Dims[objects])) < min(unlist(Ranks[x]))){
            stop(paste0("Specify the Ranks[", x, "] as smaller value"))
        }
    })
}


# Initialization
.initGCTF <- function(X, R, M, initZ, fix, Ranks, thr){
    # Size of X and M
    M <- .initM(X, M)
    # Initial Latent Indices
    Z <- .initInitZ(X, R, Ranks, initZ)
    # If fix is not set, all Z may be changed
    if(is.null(fix)){
        fix <- rep(FALSE, length(Z))
    }
    # Visible Indices e.g. I, J, K, P, Q
    visibleIdxs <- .visibleIdxs(X)
    # Latent Indices e.g. p, q, r
    latentIdxs <- .latentIdxs(Ranks, Z, visibleIdxs)
    # Reconstruction Error
    RecError = c()
    # Relative Change
    RelChange = c()
    # Reconstructed Tensor
    X_bar <- .recTensor(X, Z, R, Ranks, latentIdxs)
    # First iteration
    RecError[1] <- .recError(X, X_bar)
    RelChange[1] <- thr * 10
    list(M=M, Z=Z, fix=fix, visibleIdxs=visibleIdxs, latentIdxs=latentIdxs,
        RecError=RecError, RelChange=RelChange, X_bar=X_bar)
}

# Initialization of Mask Matrices
.initM <- function(X, M){
    if(!is.null(M)){
        if(!identical(lapply(X, dim), lapply(M, dim))){
           stop("Dimension of Mask tensor M must be the same as that of X")
        }
    }else{
        M <- X
        for(i in seq_along(M)){
            M[[i]][] <- 1
        }
    }
    M
}

# Initialization of Factor Matrices
.initInitZ <- function(X, R, Ranks, initZ){
    if(is.null(initZ)){
        Z <- .genLatentVals(X, R, Ranks)
    }else{
        checkInitZ <- all(unlist(lapply(seq_along(initZ), function(x){
          identical(as.integer(dim(initZ[[x]])),
                    as.integer(unlist(Ranks[[x]])))
        })))
        if(!checkInitZ){
            stop("The size of initZ must be same as the way Ranks specified")
        }else{
          Z <- initZ
        }
    }
    Z
}

# Visible Indices
.visibleIdxs <- function(X){
    vnames <- unique(unlist(lapply(X, function(x){
       names(dim(x))
    })))
    vdims <- unlist(lapply(X, function(x){
        dim(x)
    }))
    names(vdims) <- gsub(".*\\.", "", names(vdims))
    vdims <- vdims[vnames]
    vdims
}

# Latent indices
.latentIdxs <- function(Ranks, Z, visibleIdxs){
    lnames <- unique(unlist(lapply(Ranks, function(x){
       names(x)
    })))
    lnames <- setdiff(lnames, names(visibleIdxs))
    ldims <- unlist(lapply(Z, function(x){
        dim(x)
    }))
    names(ldims) <- gsub(".*\\.", "", names(ldims))
    ldims <- ldims[lnames]
    ldims
}

# Multiplicate Update
.MU <- function(Z, M, X, X_bar, R, Alpha, Beta, Ranks){
    # Numerator
    numer <- vapply(seq_len(nrow(R)), function(v){
        R[v, Alpha] * .TensorValdFunc(M[[v]] *
            X_bar[[v]]^(-Beta) * X[[v]], Z, R, v, Alpha, Ranks)
        }, Z[[Alpha]])
    # Denominator
    denom <- vapply(seq_len(nrow(R)), function(v){
        R[v, Alpha] * .TensorValdFunc(M[[v]] *
            X_bar[[v]]^(1-Beta), Z, R, v, Alpha, Ranks)
        }, Z[[Alpha]])
    # Update rules for GCTF
    sumRange <- seq_len(length(dim(numer)) - 1)
    out <- Z[[Alpha]] * (apply(numer, sumRange, sum) /
        apply(denom, sumRange, sum))
    out
}

.checkDimNames <- function(tensorList){
    for(tensor in tensorList){
        for(dimName in names(dim(tensor))){
            # check if length of dimName is 1
            stopifnot(nchar(dimName) == 1)
        }
    }
}

.generateEquationStringLeftHandSide <- function(tensorList){
    equation <- ""
    for(iTensor in seq_along(tensorList)){
        tensor <- tensorList[[iTensor]]
        for(dimName in names(dim(tensor))){
            equation <- paste0(equation, dimName)
        }
        if(iTensor != length(tensorList)){
            equation <- paste0(equation, ",")
        }
    }
    equation
}

.generateEquationStringRightHandSide <- function(marginNames){
    equation <- ""
    for(dimName in marginNames){
        equation <- paste0(equation, dimName)
    }
    equation
}

.generateEquationString <- function(tensorList, marginNames){
    .checkDimNames(tensorList)
    lhs <- .generateEquationStringLeftHandSide(tensorList)
    rhs <- .generateEquationStringRightHandSide(marginNames)
    equation <- paste(lhs, rhs, sep="->")
    equation
}

.MulSumTensors <- function(tensorList, marginNames){
    equation <- .generateEquationString(tensorList, marginNames)
    resultTensor <- do.call("einsum", c(equation, tensorList))
    names(dim(resultTensor)) <- marginNames
    resultTensor
}

# Tensor Valued Function
# Delta(Alpha,v) function is just computing a product of tensors
# and collapses this product over indices not appearing in Z_Alpha

# e.g.
# v = 1
# Q = M[[v]] * X_bar[[v]]^(-Beta) * X[[v]]
.TensorValdFunc <- function(Q, Z, R, v, Alpha, Ranks){
    if(R[v, Alpha] == 0){
        # Z[[Alpha]] for vapply
        out <- Z[[Alpha]]
    }else{
        notAlphas <- intersect(setdiff(seq_len(ncol(R)), Alpha),
                               which(R[v, ] == 1))
        multiplied_tensor_list <-lapply(notAlphas, function(notAlpha){Z[[notAlpha]]})
        multiplied_tensor_list[[length(multiplied_tensor_list)+1]] <- Q
        out <- .MulSumTensors(multiplied_tensor_list,
            names(dim(Z[[Alpha]])))
    }
    # Reordering
    out <- .reOrder(out, Z[[Alpha]])
    names(dim(out)) <- names(dim(Z[[Alpha]]))
    out
}

# Reconstruction Error
.recError <- function(X, X_bar){
    v <- unlist(lapply(seq_len(length(X)), function(x){
       X[[x]] - X_bar[[x]]
    }))
    sqrt(sum(v * v))
}

# Reconstructed Tensor
.recTensor <- function(X, Z, R, Ranks, latentIdxs){
    # Matrix object Position
    matObject <- lapply(Z, is.matrix)
    names(matObject) <- names(Z)
    # Reconstruction
    X_bar <- .Xbar(R, Z, matObject, latentIdxs)
    # Fix the order of X_bar as same as that of X
    .reOrder(X_bar, X)
}

# Reconstructed Tensor
.Xbar <- function(R, Z, matObject, latentIdxs){
    out <- lapply(seq_len(nrow(R)), function(x){
        # Related Factor matrix/tensor
        relFactor <- names(R[x, ][which(R[x, ] == 1)])
        # Tensor object position
        l <- which(!unlist(matObject[relFactor]))
        # if Tensor is contained in Z
        if(length(l) != 0){
            # If a tensor is seed
            seed <- Z[[names(l)]]
            relFactor <- setdiff(relFactor, names(l))
        }else{
            # If a matrix is seed
            if(length(latentIdxs) == 1){
                seed <- .diagonalTensor(names(latentIdxs),
                latentIdxs, length(relFactor))
            }else{
                firstPos <- relFactor[1]
                seed <- Z[firstPos][[1]]
                relFactor <- setdiff(relFactor, firstPos)
            }
        }
        for(i in relFactor){
            seed <- .xtx(seed, Z[i][[1]])
        }
        seed
    })
    names(out) <- rownames(R)
    out
}

.diagonalTensor <- function(latentIdxs, len1, len2){
    out <- rTensor:::.superdiagonal_tensor(len2, len1)@data
    names(dim(out)) <- rep(latentIdxs, len2)
    out
}

# Re order of X_bar by the dimensional order of X
.reOrder <- function(X_bar, X){
    X_barType <- is.array(X_bar)
    XType <- is.array(X)
    if(X_barType && XType){
        X_bar <- list(X_bar)
       X <- list(X)
    }
    vapply(seq_along(X_bar), function(x){
        Xdimnames <- names(dim(X[[x]]))
        X_barorder <- unlist(lapply(Xdimnames, function(xx){
            which(xx == names(dim(X_bar[[x]])))
        }))
        X_bar[[x]] <<- aperm(X_bar[[x]], X_barorder)
        0L
    }, 0L)
    if(X_barType && XType){
        X_bar <- X_bar[[1]]
    }
    X_bar
}

# Initial Latent Indices
.genLatentVals <- function(X, R, Ranks){
    out <- lapply(names(Ranks), function(x){
       tmp <- abs(rand_tensor(modes = unlist(Ranks[[x]]))@data)
        tmp <- tmp / norm(as.matrix(tmp), "F")
        names(dim(tmp)) <- names(Ranks[[x]])
        tmp
    })
    names(out) <- colnames(R)
    out
}

# Matrix/Tensor times Matrix/Tensor
.xtx <- function(A, B){
    isTensorA <- length(dim(A)) > 2
    isTensorB <- length(dim(B)) > 2
    # A : Matrix, B : Matrix
    if(!isTensorA && !isTensorB){
       out <- .mtm(as.matrix(A), as.matrix(B))
    }
    # A : Matrix, B : Tensor
    if(!isTensorA && isTensorB){
        out <- .ttm(B, as.matrix(A))
    }
    # A : Tensor, B : Matrix
    if(isTensorA && !isTensorB){
       out <- .ttm(A, as.matrix(B))
    }
    # A : Tensor, B : Tensor
    if(isTensorA && isTensorB){
        out <- .ttt(A, B)
    }
    out
}

# Matrix times Matrix
.mtm <- function(A, B){
    # common mode
    commonMode <- intersect(names(dim(A)), names(dim(B)))
    targetA <- vapply(commonMode, function(x){
        which(names(dim(A)) == x)[1]
    }, 0L)
    targetB <- vapply(commonMode, function(x){
        which(names(dim(B)) == x)[1]
    }, 0L)
    # not common mode
    ncModeA <- setdiff(seq_len(length(dim(A))), targetA)
    ncModeB <- setdiff(seq_len(length(dim(B))), targetB)
    if((targetA == 1) && (targetB == 1)){
        out <- t(A) %*% B
    }
    if((targetA == 1) && (targetB == 2)){
        out <- t(A) %*% t(B)
    }
    if((targetA == 2) && (targetB == 1)){
        out <- A %*% B
    }
    if((targetA == 2) && (targetB == 2)){
        out <- A %*% t(B)
    }
    names(dim(out))[1] <- names(dim(A))[ncModeA]
    names(dim(out))[2] <- names(dim(B))[ncModeB]
    out
}

# Tensor times Tensor
.ttt <- function(A, B){
    # common mode
    commonMode <- intersect(names(dim(A)), names(dim(B)))
    targetA <- vapply(commonMode, function(x){
        which(names(dim(A)) == x)[1]
    }, 0L)
    targetB <- vapply(commonMode, function(x){
        which(names(dim(B)) == x)[1]
    }, 0L)
    # not common mode
    ncModeA <- setdiff(seq_len(length(dim(A))), targetA)
    ncModeB <- setdiff(seq_len(length(dim(B))), targetB)
    # Matricising
    mat_1 <- rTensor::unfold(as.tensor(A),
        row_idx = ncModeA, col_idx = targetA)@data
    mat_2 <- rTensor::unfold(as.tensor(B),
        row_idx = ncModeB, col_idx = targetB)@data
    # Matrix-times-matrix
    out <- mat_1 %*% t(mat_2)
    if(length(ncModeA) == 1){
        names(dim(out))[1] <- names(dim(A)[ncModeA])
    }
    if(length(ncModeB) == 1){
        names(dim(out))[2] <- names(dim(B)[ncModeB])
    }
    if(length(ncModeA) == 0 && length(ncModeB) > 0){
        # if ncModeA is empty
        names(ncModeB) <- names(dim(B))[ncModeB]
        out <- array(out, dim(B)[ncModeB])
    }else if(length(ncModeB) == 0 && length(ncModeA) > 0){
        # if ncModeB is empty
        names(ncModeA) <- names(dim(A))[ncModeA]
        out <- array(out, dim(A)[ncModeA])
    }else if(length(ncModeA) == 0 && length(ncModeB) == 0){
        # if both ncMode1 and ncMode2 are empty
        out <- out
    }else if(length(ncModeA) == 1 && length(ncModeB) == 1){
        # both ncMode1 and ncMode2 contains only one element
        names(dim(out)) <- c(names(dim(A))[ncModeA], names(dim(B))[ncModeB])
    }else{
        # both ncMode1 and ncMode2 contains elements
        out <- rTensor::fold(
            out, row_idx=seq_along(ncModeA),
            col_idx=seq(from=length(ncModeA)+1,
            to=length(ncModeA)+length(ncModeB)),
            modes=c(ncModeA, ncModeB))
    }
    out
}

.ttm <- function(A, B){
    # common mode
    commonMode <- intersect(names(dim(A)), names(dim(B)))
    targetA <- vapply(commonMode, function(x){
        which(names(dim(A)) == x)[1]
    }, 0L)
    targetB <- vapply(commonMode, function(x){
        which(names(dim(B)) == x)[1]
    }, 0L)
    notComModeA <- setdiff(seq_len(length(dim(A))), targetA)
    notComModeB <- setdiff(seq_len(length(dim(B))), targetB)
    A <- aperm(A, c(targetA, notComModeA))
    B <- aperm(B, c(targetB, notComModeB))
    # Tensor-times-matrix
    # need support for the case modes of A and B are same.
    out <- ttm(as.tensor(A), t(B), m=1)
    out <- out@data
    names(dim(out))[1] <- names(dim(B)[2])
    vapply(2:length(dim(out)), function(x){
        names(dim(out))[x] <<- names(dim(A)[x])
    }, "")
    out
}

# Relative Change
.relChange <- function(iter, RecError){
    abs(RecError[iter-1] - RecError[iter]) / RecError[iter]
}
