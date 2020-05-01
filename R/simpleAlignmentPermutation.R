rm(list = ls())

loess.as = function(x, y, degree=1, criterion=c("aicc", "gcv"), 
                    family = c("gaussian", "symmetric"), user.span=NULL, plot=FALSE, ...) {
    criterion <- match.arg(criterion)
    family <- match.arg(family)
    x <- as.matrix(x)
    
    if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
    if (!is.numeric(x)) stop("argument 'x' must be numeric!")
    if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) 
        stop("argument 'user.span' must be a numerical number!")
    if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
    if(length(y) < 3) stop("not enough observations!")
    
    data.bind <- data.frame(x=x, y=y)
    if (ncol(x) == 1) {
        names(data.bind) <- c("x", "y")
    } else { names(data.bind) <- c("x1", "x2", "y") }
    
    opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
        as.crit <- function (x) {
            span <- x$pars$span
            traceL <- x$trace.hat
            sigma2 <- sum(x$residuals^2 ) / (x$n-1)
            aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
            gcv <- x$n*sigma2 / (x$n-traceL)^2
            result <- list(span=span, aicc=aicc, gcv=gcv)
            return(result)
        }
        criterion <- match.arg(criterion)
        fn <- function(span) {
            mod <- update(model, span=span)
            as.crit(mod)[[criterion]]
        }
        result <- optimize(fn, span.range)
        return(list(span=result$minimum, criterion=result$objective))
    }
    
    if (ncol(x)==1) {
        if (is.null(user.span)) {
            fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
            span1 <- opt.span(fit0, criterion=criterion)$span
        } else {
            span1 <- user.span
        }		
        fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
    } else {
        if (is.null(user.span)) {
            fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
            span1 <- opt.span(fit0, criterion=criterion)$span
        } else {
            span1 <- user.span
        }		
        fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
    }
    if (plot){
        if (ncol(x)==1) {
            m <- 100
            x.new <- seq(min(x), max(x), length.out=m)
            fit.new <- predict(fit, data.frame(x = x.new))
            plot(x, y, col="lightgrey", xlab="x", ylab="m(x)", ...)
            lines(x.new,fit.new, lwd=1.5, ...)
        } else {
            m <- 50
            x1 <- seq(min(data.bind$x1), max(data.bind$x1), len=m) 
            x2 <- seq(min(data.bind$x2), max(data.bind$x2), len=m) 
            x.new <- expand.grid(x1=x1, x2=x2) 
            fit.new <- matrix(predict(fit, x.new), m, m) 
            persp(x1, x2, fit.new, theta=40, phi=30, ticktype="detailed", xlab="x1", ylab="x2", zlab="y", col="lightblue", expand=0.6)
        }		
    }
    return(fit)
}

findComparableRTs = function(fMeanMz, fMeanRt, fZ, fMeanIntensity, libMz, libRt, libZ) {
    # Find comparable features and library entries
    res = NULL
    for (i in 1:nFeatures) {
        if (!is.na(fZ[i])) {
            mzDiff = abs(fMeanMz[i] - libMz) / libMz * 1e6
            if (fZ[i] == 0) {
                ind = which(mzDiff < matchMzTol & !is.na(libRt))
            } else {
                ind = which(mzDiff < matchMzTol & !is.na(libRt) & fZ[i] == libZ)
            }
            if (length(ind) == 1) {
                res = rbind(res, c(fMeanRt[i], libRt[ind]))
            } else if (length(ind) > 1) {
                indMax = which.max(fMeanIntensity[ind])
                ind = ind[indMax]
                res = rbind(res, c(fMeanRt[i], libRt[ind]))
            } else {
                next
            }
        }
    }
    return (res)    # column1 = feature RT, column2 = library RT
}

rtAlignment = function(RTs, featureRT) {
    x = as.numeric(RTs[, 1])    # Feature RT
    y = as.numeric(RTs[, 1]) - as.numeric(RTs[, 2]) ## RT-shifts (featureRT - libraryRT)
    
    # Removal of outliers in y (i.e. RT-shift)
    pct = 0.2
    truncatedMean = mean(y[y >= quantile(y, pct / 2) & y <= quantile(y, (1 - pct / 2))])
    truncatedSd = sd(y[y >= quantile(y, pct / 2) & y <= quantile(y, (1 - pct / 2))])
    lL = truncatedMean - 3 * truncatedSd
    uL = truncatedMean + 3 * truncatedSd
    ind = which(y >= lL & y <= uL)
    
    mod = loess.as(x[ind], y[ind], degree = 1, criterion = "aicc",
                   control = loess.control(surface = "direct")) ## This curve represents RT-shifts as a function of comp$Rt
    newFeatureRT = as.numeric(featureRT) - predict(mod, data.frame(x = as.numeric(featureRT)))
    
    res = list(residual = abs(mod$residual), newFeatureRT = newFeatureRT)
    return (res)
}

ms2Similarity = function(fSpec, lSpec) {
    k = min(30, min(dim(fSpec)[1], dim(lSpec)[1]))
    fSpec = fSpec[order(fSpec$intensity, decreasing = T), ]
    fSpec = fSpec[1:30, ]
    lSpec = lSpec[order(lSpec$intensity, decreasing = T), ]
    lSpec = lSpec[1:30, ]
    
    mz = sort(union(fSpec$mz, lSpec$mz))
    mzArray = NULL
    val = 0
    j = 0
    for (i in 1:length(mz)) {
        if (abs(mz[i] - val) <= 0.5) {
            mzArray = rbind(mzArray, c(mz[i], j))
        } else {
            j = j + 1
            mzArray = rbind(mzArray, c(mz[i], j))
            val = mz[i]
        }
    }
    
    fIntensity = rep(0, max(mzArray[, 2]))
    lIntensity = rep(0, max(mzArray[, 2]))
    for (i in 1:dim(mzArray)[1]) {
        fInd = which(fSpec$mz == mzArray[i, 1])
        if (length(fInd) > 0) {
            j = mzArray[i, 2]
            fIntensity[j] = fIntensity[j] + sqrt(fSpec$intensity[fInd])
        }
        lInd = which(lSpec$mz == mzArray[i, 1])
        if (length(lInd) > 0) {
            j = mzArray[i, 2]
            lIntensity[j] = lIntensity[j] + sqrt(lSpec$intensity[lInd])
        }
    }
    
    den = sum(fIntensity * lIntensity)
    num1 = sum(fIntensity ^ 2)
    num2 = sum(lIntensity ^ 2)
    
    if (num1 * num2 == 0) {
        normDotProduct = 0
    } else {
        normDotProduct = den / (sqrt(num1 * num2))
    }
    return (normDotProduct)
}

findComparableSpectra = function(fMeanMz, fZ, fMS2, libMz, libZ, libMS2) {
    # Find comparable MS2 spectra of features and library entries
    featList = list()
    libList = list()
    n = 1
    for (i in 1:nFeatures) {
        if (!is.na(fZ[i]) & prod(!is.na(fMS2[[i]]))) {
            for (j in 1:length(libMz)) {
                if (!is.na(libZ[j]) & prod(!is.na(libMS2[[j]]))) {
                    if (fZ[i] > 0 & libZ[j] > 0 & fZ[i] != libZ[j]) {
                        next
                    } else {
                        mzDiff = abs(fMeanMz[i] - libMz[j]) / libMz[j] * 1e6
                        if (mzDiff < matchMzTol) {
                            featList[[n]] = fMS2[[i]]
                            libList[[n]] = libMS2[[j]]
                            n = n + 1
                        }
                    }
                }
            }
        }
    }
    res = list(feat = featList, lib = libList)
    return (res)
}

########
# Main #
########
# Parameter
mode = -1
condition = "hilic"
if (mode == -1) {
    columnName = paste0(condition, "n")
} else {
    columnName = paste0(condition, "p")
}
proton = 1.007276466812
matchMzTol = 10    # ppm
libraryFile = "/Research/Projects/7Metabolomics/Library/Metabolome_library_v3.1.1.txt"
featureFile = "/Research/Projects/7Metabolomics/Dev/Jumpm_unlabel_python/IROA_IS_NEG/.IROA_IS_NEG_fully_aligned.feature"
featureMs2Path = "/Research/Projects/7Metabolomics/Dev/Jumpm_unlabel_python/IROA_IS_NEG/MS2/"

# Load library information
libDf = read.csv(libraryFile, sep = "\t")
libRt = as.numeric(as.character(libDf[[paste0(columnName, "_rt")]]))
libRt = libRt * 60
libZ = as.numeric(as.character(libDf[[paste0(columnName, "_charge")]]))
libM = as.numeric(as.character(libDf[["monoisotopic_mass"]]))
if (mode == -1) {
    libMz = (libM - libZ * proton) / libZ
} else {
    libMz = (libM + libZ * proton) / libZ
}
libMS2 = NULL
for (i in 1:dim(libDf)[1]) {
    ms2File = as.character(libDf[[paste0(columnName, "_linkms2")]][i])
    
    
    #############################################################################################################
    # For test purpose
    ms2File = gsub("/research/rgs01/applications/hpcf/authorized_apps/proteomics_apps/jumpm/library/v02/hilicn/",
                   "/Research/Projects/7Metabolomics/Library/hilicn/", ms2File)
    #############################################################################################################    
    
    
    if (file.exists(ms2File)) {
        df = read.table(ms2File, sep = "\t")
        df = df[2:nrow(df), ]
        colnames(df) = c("mz", "intensity")
        libMS2[[i]] = df
    } else {
        libMS2[[i]] = NA
    }
}

# Load feature information
featDf = read.csv(featureFile, sep = "\t")
fMeanMz = as.numeric(as.character(featDf[["meanMz"]]))
fMeanRt = as.numeric(apply(featDf[,grep("_RT", colnames(featDf))], 1, mean))
fMeanIntensity = as.numeric(apply(featDf[,grep("_Intensity", colnames(featDf))], 1, mean))
fZ = NULL
fMS2 = NULL
for (i in 1:dim(featDf)[1]) {
    ms2File = paste0(featureMs2Path, "/f", i, ".MS2")
    if (file.exists(ms2File)) {
        df = read.table(ms2File, sep = "\t")
        fZ[i] = df[1, 2]
        df = df[2:nrow(df), ]
        colnames(df) = c("mz", "intensity")
        fMS2[[i]] = df
    } else {
        fZ[i] = NA
        fMS2[[i]] = NA
    }
}
nFeatures = dim(featDf)[1]

################
# RT-alignment #
################
# RT-alignment between observed features and library
obsRt = findComparableRTs(fMeanMz, fMeanRt, fZ, fMeanIntensity, libMz, libRt, libZ)
obsRes = rtAlignment(obsRt, fMeanRt)

# RT-alignment with permutations and the estimation of empirical null distribution
permRes = NULL
nPerms = 10
pb = txtProgressBar(1, nPerms, width = 20, style = 3)
for (i in 1:nPerms) {
    setTxtProgressBar(pb, i)
    permRt = findComparableRTs(sample(fMeanMz), fMeanRt, fZ, fMeanIntensity, libMz, libRt, libZ)
    res = rtAlignment(permRt, fMeanRt)
    permRes = c(permRes, as.numeric(res$residual))
}

# Empirical null distribution of RT-shifts
f_rt = ecdf(permRes)

########################################
# Match features and library compounds #
########################################
# res = NULL
# for (i in 1:nFeatures) {
#     if (!is.na(fZ[i]) & prod(!is.na(fMS2[[i]]))) {
#         for (j in 1:length(libMz)) {
#             if (!is.na(libZ[j]) & prod(!is.na(libMS2[[j]]))) {
#                 if (fZ[i] > 0 & libZ[j] > 0 & fZ[i] != libZ[j]) {
#                     next
#                 } else {
#                     mzDiff = abs(fMeanMz[i] - libMz[j]) / libMz[j] * 1e6
#                     if (mzDiff < matchMzTol) {
#                         # Calculation of MS2 similarity
#                         sim = ms2Similarity(fMS2[[i]], libMS2[[j]])
#                         rtDiff = abs(fMeanRt[i] - libRt[j])
#                         res = rbind(res, c(i, j, sim, rtDiff))
#                     }
#                 }
#             }
#         }
#     }
# }

# Observed MS2-based similarity
res = findComparableSpectra(fMeanMz, fZ, fMS2, libMz, libZ, libMS2)
obsMS2 = NULL
for (i in 1:length(res$feat)) {
    ms2Sim = ms2Similarity(res$feat[[i]], res$lib[[i]])
    obsMS2 = c(obsMS2, ms2Sim)
}

# Permuted MS2-based similarity (null)
pb = txtProgressBar(1, nPerms, width = 20, style = 3)
permMS2 = NULL
for (i in 1:nPerms) {
    setTxtProgressBar(pb, i)
    res = findComparableSpectra(sample(fMeanMz), fZ, fMS2, libMz, libZ, libMS2)
    for (j in 1:length(res$feat)) {
        ms2Sim = ms2Similarity(res$feat[[j]], res$lib[[j]])
        permMS2 = c(permMS2, ms2Sim)
    }
}

# Empirical null distribution of RT-shifts
f_ms2 = ecdf(permMS2)
