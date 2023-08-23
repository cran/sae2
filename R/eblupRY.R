eblupRY <- function(formula, D, TI, vardir, method = c("REML", "ML"),
         MAXITER = 1000, PRECISION = .1e-05, data, max.rho = .98, 
         dampening = 0.9, ...) { 
         
    if (length(method) > 1) method <- method[1]
    if (!method %in% c("REML", "ML", "MLE")) 
        stop(" method=\"", method, "\" must be \"REML\", or \"ML\"")
    if(inherits(formula, "list")) {
        NV <- length(formula)
        if(NV == 1) formula <- formula[[1]]
    } else { 
        NV <- 1
    }
    if (inherits(vardir, "list")) {
        if (missing(D)) { 
            D <- length(vardir) 
        } else {
            if (D != length(vardir)) 
                stop("the length of 'vardir' must agree with D")
        }
        lapply(vardir, FUN = function(x) { 
            if (!all(dim(x) == c(NV*TI, NV*TI))) 
                stop(paste("each element of 'vardir' must be a square matrix",
                           "with 'NV*TI' rows"))
        })
        vardir.temp <- matrix(0, nrow=D*NV*TI, ncol=D*NV*TI)
        for (d in 1:D) {
            vardir.temp[((d-1)*NV*TI+1):(d*NV*TI), ((d-1)*NV*TI+1):(d*NV*TI)] <-
                vardir[[d]]
        }
        vardir <- vardir.temp
    } else {
        if (missing(D)) stop("'D' must be specified")
        if (!is.matrix(vardir)) vardir <- as.matrix(vardir)
        if (!inherits(vardir, "matrix")) 
            stop("'vardir' must be a matrix or a list of matrices")
        if (!all(dim(vardir) == c(D*NV*TI, D*NV*TI)))
            stop("'vardir' must be a square matrix with 'D*NV*TI' rows")
    }
    if (NV == 1) {
        if (!missing(data)) {
            formuladata <- model.frame(formula, na.action = NULL, data)
        } else {
            formuladata <- model.frame(formula, na.action = NULL)
        }
        X <- model.matrix(formula, data=formuladata)
        ncolx <- dim(X)[2]                                  # ncolx added
        nformula <- nrow(X)
        if(nformula != D*TI) 
            stop("length of variables must be D * TI")
        y <- formuladata[, 1]
        result <- dynRYfit(y, X, M=D, TI=TI, NV=1, vcov_e=vardir,
                              maxiter=MAXITER, iter.tol=PRECISION, 
                              ncolx=ncolx, method=method, max.rho=max.rho, 
                              dampening=dampening, model="RY", ...)
    } else {
        depvarlist <- " "
        formuladata <- list()
        X.list <- list()
        ncolx <- 0
        X.names <- NULL
        for (nv in 1:NV) {
            formula1 <- formula[[nv]]
            depvarlist[nv] <- as.character(formula1[2])
            if (!missing(data)) {
                formuladata1 <- model.frame(formula1, na.action = NULL, data)
            } else {
                formuladata1 <- model.frame(formula1, na.action = NULL)
            }
            X1 <- model.matrix(formula1, data=formuladata1)
            formuladata[[nv]] <- formuladata1
            X.list[[nv]] <- X1
            ncolx[nv] <- dim(X1)[2]
            X.names <- append(X.names, paste0(attr(X1,"dimnames")[[2]], ".", nv))
        }
        y <- rep(0, D*TI*NV)
        X <- matrix(0, nrow=D*TI*NV, ncol=sum(ncolx))
        attr(X, "dimnames")[[2]] <- X.names
        cstart <- 1
        for (nv in 1:NV) {
            X1 <- X.list[[nv]]
            y1 <- formuladata[[nv]][, 1]
            cend <- cstart + ncol(X1) - 1
            nformula <- nrow(X1)
            if(nformula != D*TI) 
                stop("length of variables must be D * TI")
            for (m in 1:D) {
                X[((m-1)*TI*NV+(nv-1)*TI+1):((m-1)*TI*NV+nv*TI), cstart:cend] <- 
                    X1[((m-1)*TI+1):(m*TI), ]
                y[((m-1)*TI*NV+(nv-1)*TI+1):((m-1)*TI*NV+nv*TI)] <-      
                    y1[((m-1)*TI+1):(m*TI)] 
            }
            cstart <- cend + 1
        }
        result <- dynRYfit(y, X, M=D, TI=TI, NV=NV, vcov_e=vardir,
                              maxiter=MAXITER, iter.tol = PRECISION, 
                              ncolx=ncolx, method=method, max.rho=max.rho, 
                              dampening=dampening, model="RY", ...) 
        for (x in c("eblup", "eblup.mse", "eblup.g1", "eblup.g2", "eblup.g3", 
                    "est.fixed", "est.fixed.var", "eblup.wt1", "eblup.wt2")) {
            if (x %in% names(result)) {
                colnames(result[[x]]) <- depvarlist
                result[[x]] <- as.data.frame(result[[x]])
            }
        }
    }

    names(result$fit$iterations) <- NULL
    if("var.coef" %in% names(result)) {
        std.error <- sqrt(diag(result$var.coef))
        tvalue <- result$coef/std.error
        pvalue <- ifelse (tvalue < 0, 2*pnorm(tvalue), 2*(1-pnorm(tvalue)))
        result$fit$estcoef <- data.frame(beta=result$coef, std.error=std.error,
                                  tvalue=tvalue, pvalue=pvalue,
                                  row.names=names(result$coef))
        std.error <- sqrt(diag(solve(result$inf.mat)))
        result$fit$estvarcomp <- data.frame(estimate=result$delta, 
                                     std.error=std.error,
                                     row.names=names(result$delta))
    }
    if ("parm" %in% names(result)) {
        result$fit$iterations <- result$parm["num.iter"]
 	    if (method == "REML") {
            goodness <- c(result$parm["loglikelihood"], 
                          result$parm["constrained.ll"])
            names(goodness) <- c("loglike", "restrictedloglike")
        } else {
            goodness <- result$parm["loglikelihood"]
            names(goodness) <- "loglike"
        }
        result$fit$goodness <- goodness
    }
    result$model <- formula
    return(result)
}           
    
    