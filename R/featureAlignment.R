rm(list = ls())

# Input: individual feature files
# Output: fully-aligned feature file

getParams = function (paramFile) {
    lines = readLines(paramFile)
    params = list()
    for (i in 1:length(lines)) {
        lines[i] = gsub("\\s", "", lines[i])
        if (grepl("^#", lines[i]) | !grepl("=", lines[i])) {
            next
        } else {
            line = gsub("#.*", "", lines[i])
            key = unlist(strsplit(line, "="))[1]
            val = unlist(strsplit(line, "="))[2]
            params[[key]] = val
        }
    }
    return (params)
}

calibrateFeatures = function (ref, comp, params, LOG) {
    initMzTol = as.numeric(params$"tol_initial")
    sdWidth = as.numeric(params$"sd_width")
    
    ref = ref[order(ref$intensity, decreasing = TRUE), ] ## Sort by intensity
    comp = comp[order(comp$intensity, decreasing = TRUE), ] ## Sort by intensity
    
    ########################################
    ## Calibration of RT and m/z globally ##
    ########################################
    ## Calculate rtShift and mzShift in a global manner, and calibrate features in the compared sample
    cat("  Global calibration of RT and m/z is being performed\n")
    cat("  Global calibration of RT and m/z is being performed\n", file = LOG)
    res = globalCalibration(ref, comp, initMzTol)
    rtShifts = res$rtShifts ## unit of second
    mzShifts = res$mzShifts ## unit of ppm
    cat(sprintf("    Based on the matched features within +/- %i ppm\n", initMzTol))
    cat(sprintf("    The global RT-shift is %.4f second\n", median(rtShifts)))
    cat(sprintf("    The global m/z-shift is %.4f ppm\n", median(mzShifts)))
    cat("    RT and m/z of the compared sample is calibrated according to the global RT- and m/z-shift\n")
    cat(sprintf("    Based on the matched features within +/- %i ppm\n", initMzTol), file = LOG)
    cat(sprintf("    The global RT-shift is %.4f second\n", median(rtShifts)), file = LOG)
    cat(sprintf("    The global m/z-shift is %.4f ppm\n", median(mzShifts)), file = LOG)
    cat("    RT and m/z of the compared sample is calibrated according to the global RT- and m/z-shift\n", file = LOG)
    comp$RT = comp$RT - median(rtShifts)
    comp$mz = comp$mz / (1 + median(mzShifts)/ 1e6)
    
    ##########################################
    ## Calibration of RT using LOESS curve  ##
    ##########################################
    cat("  Local calibration of RT and m/z is being performed (through LOESS modeling)\n")
    cat("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows\n")
    cat("      RT- and m/z-tolerance =", sdWidth, "x dynamically estimated SD of RT- and m/z-shifts\n")
    cat("    LOESS modeling may take some time. Please be patient ...\n")
    cat("  Local calibration of RT and m/z is being performed (through LOESS modeling)\n", file = LOG)
    cat("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows\n", file = LOG)
    cat("      RT- and m/z-tolerance =", sdWidth, "x dynamically estimated SD of RT- and m/z-shifts\n", file = LOG)
    cat("    LOESS modeling may take some time. Please be patient ...\n", file = LOG)
    ## 1st round of LOESS modeling and calibration of RT
    ## When performing this LOESS modeling,
    ## the matched features are selected based on the global rt- and mz-shift
    rtSd = max(sd(rtShifts), 1e-3)
    mzSd = max(sd(mzShifts), 1e-3)
    res = localCalibration(ref, comp, rtSd, mzSd, sdWidth, "RT")
    estRtShift = predict(res$model, data.frame(x = comp$RT))
    comp$RT = comp$RT - estRtShift
    cat("    The 1st round of RT-calibration is done\n")
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n")
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n")
    cat("    The 1st round of RT-calibration is done\n", file = LOG)
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n", file = LOG)
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n", file = LOG)
    
    ## 2nd round of LOESS modeling and calibration of RT
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = localCalibration(ref, comp, rtSd, mzSd, sdWidth, "RT")
    estRtShift = predict(res$model, data.frame(x = comp$RT))
    comp$RT = comp$RT - estRtShift
    cat("    The 2nd round of RT-calibration is done\n")
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n")
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n")
    cat("    The 2nd round of RT-calibration is done\n", file = LOG)
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n", file = LOG)
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n", file = LOG)
    
    ##########################################
    ## Calibration of m/z using LOESS curve ##
    ##########################################
    ## 1st round of LOESS modeling and calibration of m/z
    ## When performing this LOESS modeling,
    ## the matched features are selected based on the global rt- and mz-shift
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = localCalibration(ref, comp, rtSd, mzSd, sdWidth, "m/z")
    estMzShift = predict(res$model, data.frame(x = comp$mz))
    comp$mz = comp$mz / (1 + estMzShift / 1e6)
    cat("    The 1st round of m/z-calibration is done\n")
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n")
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n")
    cat("    The 1st round of m/z-calibration is done\n", file = LOG)
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n", file = LOG)
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n", file = LOG)
    
    ## 2nd round of LOESS modeling and calibration of m/z
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = localCalibration(ref, comp, rtSd, mzSd, sdWidth, "m/z")
    estMzShift = predict(res$model, data.frame(x = comp$mz))
    comp$mz = comp$mz / (1 + estMzShift / 1e6)
    cat("    The 2nd round of m/z-calibration is done\n")
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n")
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n")
    cat("    The 2nd round of m/z-calibration is done\n", file = LOG)
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n", file = LOG)
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n", file = LOG)
    
    ## Output organization
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = list(calibratedFeature = comp, rtSd = rtSd, mzSd = mzSd)
    return (res)
}

globalCalibration = function(ref, comp, mzTol = 20) {
    nPeaks = round(0.05 * dim(ref)[1]) ## Number of peaks "matched" between a reference and a compared sample
    i = 1 ## row index of peaks in the reference sample. The 1st peak is the strongest peak
    j = 1 ## index of the matched peaks
    subDf = data.frame()
    while (j <= nPeaks) {
        z = ref$z[i]
        mz = ref$mz[i]
        rt = ref$RT[i]
        intensity = ref$intensity[i]
        if (z == 0) {
            ## For a reference feature with undetermined charge, consider all possible charges
            rowInd = which(comp$mz >= mz - mz * mzTol / 1e6 & comp$mz < mz + mz * mzTol / 1e6)
        } else {
            ## For a reference feature with a specific charge, only consider features with the same charge
            rowInd = which(comp$mz >= mz - mz * mzTol / 1e6 & comp$mz < mz + mz * mzTol / 1e6 & comp$z == z)
        }
        if (length(rowInd) > 0) {
            rowInd = rowInd[1]
            subDf[j, 1] = mz
            subDf[j, 2] = rt
            subDf[j, 3] = intensity
            subDf[j, 4] = comp$mz[rowInd]
            subDf[j, 5] = comp$RT[rowInd]
            subDf[j, 6] = comp$intensity[rowInd]
            comp = comp[-rowInd, ]
            j = j + 1
        }
        i = i + 1
    }
    colnames(subDf) = c("refMz", "refRT", "refIntensity", "compMz", "compRT", "compIntensity")
    
    ## Calculate the global shift and its standard deviation using
    ## 10% trimmed rt- and mz-shift values
    rtShifts = subDf$"compRT" - subDf$"refRT"
    rtInd = which(rtShifts >= quantile(rtShifts, 0.1) & rtShifts <= quantile(rtShifts, 0.9))
    rtShifts = rtShifts[rtInd]
    mzShifts = (subDf$"compMz" - subDf$"refMz") / subDf$"refMz" * 1e6
    mzInd = which(mzShifts >= quantile(mzShifts, 0.1) & mzShifts <= quantile(mzShifts, 0.9))
    mzShifts = mzShifts[mzInd]
    
    res = list(rtShifts = rtShifts, mzShifts = mzShifts)
    return (res)
}

localCalibration = function(ref, comp, rtSd, mzSd, sdWidth, cal) {
    ## For each feature the reference sample, look for the corresponding feature
    ## to be aligned in the compared samples
    subDf = data.frame()
    n = dim(ref)[1]
    if (length(rtSd) == 1) {
        rtSd = rep(rtSd, n)
    }
    if (length(mzSd) == 1) {
        mzSd = rep(mzSd, n)
    }
    rtTol = sdWidth * rtSd
    mzTol = sdWidth * mzSd
    j = 1
    for (i in 1:n) {
        z = ref$z[i]
        mz = ref$mz[i]
        rt = ref$RT[i]
        intensity = ref$intensity[i]
        rtErr = comp$RT - rt
        mzErr = (comp$mz - mz) / mz * 1e6 ## unit of ppm
        if (z == 0) {
            ## For a reference feature with undetermined charge, consider all possible charges
            rowInd = which(abs(rtErr) <= rtTol[i] & abs(mzErr) <= mzTol[i])
        } else {
            ## For a reference feature with a specific charge, only consider features with the same charge
            rowInd = which(abs(rtErr) <= rtTol[i] & abs(mzErr) <= mzTol[i] & comp$z == z)
        }
        if (length(rowInd) > 0) {
            ## When multiple peaks correspond to a peak in the reference sample
            ## choose one (with the highest intensity)
            ## Note that the data.frame(comp) is already sorted according to its intensity values
            rowInd = rowInd[1]
            subDf[j, 1] = mz
            subDf[j, 2] = rt
            subDf[j, 3] = intensity
            subDf[j, 4] = comp$mz[rowInd]
            subDf[j, 5] = comp$RT[rowInd]
            subDf[j, 6] = comp$intensity[rowInd]
            j = j + 1
        }
    }
    colnames(subDf) = c("refMz", "refRT", "refIntensity", "compMz", "compRT", "compIntensity")
    
    if (cal == "RT") {
        ## Build a LOESS model for the RT-shifts between the reference and compared samples
        compRt = subDf$"compRT"
        refRt = subDf$"refRT"
        rtShifts = compRt - refRt
        if (sum(rtShifts == 0) == length(rtShifts)) {
            rtShifts = 1e-6 * rnorm(length(rtShifts))
        }
        mod = loess.as(compRt, rtShifts, degree = 1, criterion = "aicc",
                       control = loess.control(surface = "direct"))
        
        ## Calculate a new rt-tolerance using trimming
        compRt = compRt - mod$fitted ## Calibration according to the LOESS model
        rtShifts = (compRt - refRt)
        rtInd = which(rtShifts >= quantile(rtShifts, 0.1) & rtShifts <= quantile(rtShifts, 0.9))
        modRtSd = loess.as(compRt[rtInd], rtShifts[rtInd] ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynRtSd = sqrt(pmax(0, predict(modRtSd, data.frame(x = ref$RT))))
        statRtSd = sd(rtShifts[rtInd])
        
        ## Calculate a new mz-tolerance using trimming
        ## Sometimes, the variation of m/z-shifts cannot be captured when trimming is applied
        ## So, the trimming is not used for m/z-shifts
        # mzInd = which(mzShifts >= quantile(mzShifts, 0.1) & mzShifts <= quantile(mzShifts, 0.9))
        mzShifts = (subDf$"compMz" - subDf$"refMz") / subDf$"compMz" * 1e6
        if (sum(mzShifts == 0) == length(mzShifts)) {
            mzShifts = 1e-6 * rnorm(length(mzShifts))
        }
        modMzSd = loess.as(subDf$"compMz", mzShifts ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynMzSd = sqrt(pmax(0, predict(modMzSd, data.frame(x = ref$mz))))
        statMzSd = sd(mzShifts)
    } else if (cal == "m/z") {
        ## Build a LOESS model for the mz-shifts between the reference and compared samples
        compMz = subDf$"compMz"
        refMz = subDf$"refMz"
        mzShifts = (compMz - refMz) / compMz * 1e6
        if (sum(mzShifts == 0) == length(mzShifts)) {
            mzShifts = 1e-6 * rnorm(length(mzShifts))
        }
        mod = loess.as(compMz, mzShifts, degree = 1, criterion = "aicc",
                       control = loess.control(surface = "direct"))
        
        ## Calculate a new mz-tolerance using trimming
        ## Sometimes, the variation of m/z-shifts cannot be captured when trimming is applied
        ## So, the trimming is not used for m/z-shifts
        # mzInd = which(mzShifts >= quantile(mzShifts, 0.1) & mzShifts <= quantile(mzShifts, 0.9))
        compMz = compMz * (1 + mod$fitted / 1e6) ## Calibration according to the LOESS model
        mzShifts = (compMz - refMz) / compMz * 1e6
        if (sum(mzShifts == 0) == length(mzShifts)) {
            mzShifts = 1e-6 * rnorm(length(mzShifts))
        }
        modMzSd = loess.as(subDf$"compMz", mzShifts ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynMzSd = sqrt(pmax(0, predict(modMzSd, data.frame(x = ref$mz))))
        statMzSd = sd(mzShifts)
        
        ## Calculate a new rt-tolerance using trimming
        rtShifts = (subDf$"compRT" - subDf$"refRT")
        if (sum(rtShifts == 0) == length(rtShifts)) {
            rtShifts = 1e-6 * rnorm(length(rtShifts))
        }
        rtInd = which(rtShifts >= quantile(rtShifts, 0.1) & rtShifts <= quantile(rtShifts, 0.9))
        modRtSd = loess.as(subDf$"compRT"[rtInd], rtShifts[rtInd] ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynRtSd = sqrt(pmax(0, predict(modRtSd, data.frame(x = ref$RT))))
        statRtSd = sd(rtShifts[rtInd])
    }
    
    res = list(model = mod, dynRtSd = dynRtSd, statRtSd = statRtSd,
               dynMzSd = dynMzSd, statMzSd = statMzSd)
    return (res)
}

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

matchFeatures = function (ref, comp, refName, compName, rtSd, mzSd, params, LOG) {
    sdWidth = as.numeric(params$"sd_width")
    rescue = as.numeric(params$"rescue")

    ## Sort features by their intensities (to consider high intensity features first)
    ref = ref[order(ref$intensity, decreasing = T), ]
    comp = comp[order(comp$intensity, decreasing = T), ]
    
    n = dim(ref)[1]
    nc = dim(comp)[1]
    cat(paste0("  ", refName, ": ", n, " features (reference run)\n"))
    cat(paste0("  ", compName, ": ", nc, " features (compared run)\n"))
    cat(paste0("  ", refName, ": ", n, " features (reference run)\n"), file = LOG)
    cat(paste0("  ", compName, ": ", nc, " features (compared run)\n"), file = LOG)
    
    dfAligned = data.frame()
    dfUnaligned = data.frame()
    if (length(rtSd) == 1) {
        rtSd = rep(rtSd, n)
    }
    if (length(mzSd) == 1) {
        mzSd = rep(mzSd, n)
    }
    rtSd[rtSd == 0] = min(rtSd[rtSd > 0]) ## Prevent zero-tolerance of RT
    mzSd[mzSd == 0] = min(mzSd[mzSd > 0]) ## Prevent zero-tolerance of m/z
    rtTol = sdWidth * rtSd
    mzTol = sdWidth * mzSd
    
    j = 1
    mzShifts = NULL
    rtShifts = NULL
    for (i in 1:n) {
        mz = ref$mz[i]
        rt = ref$RT[i]
        intensity = ref$intensity[i]
        charge = ref$z[i]
        rtErr = comp$RT - rt
        mzErr = (comp$mz - mz) / mz * 1e6
        rowInd = which(abs(rtErr) <= rtTol[i] & abs(mzErr) <= mzTol[i])
        
        ## When there is/are matched feature(s)
        if (length(rowInd) > 0) {
            ## Check charge states of the matched features
            ## 0 charge state can be matched to any charge state
            if (charge == 0) {
                rowInd = rowInd[1]
                dfAligned = rbind(dfAligned, cbind(ref[i, ], comp[rowInd, ]))
                comp = comp[-rowInd, ]
            } else {
                if (comp[rowInd[1], ]$z == 0) { ## When the compared features has charge 0, it is okay
                    rowInd = rowInd[1]
                    dfAligned = rbind(dfAligned, cbind(ref[i, ], comp[rowInd, ]))
                    comp = comp[-rowInd, ]
                } else {
                    chargeInd = which(comp$z[rowInd] == charge)
                    if (length(chargeInd) > 0) {
                        rowInd = rowInd[chargeInd[1]]
                        dfAligned = rbind(dfAligned, cbind(ref[i, ], comp[rowInd, ]))
                        comp = comp[-rowInd, ]
                    }
                }
            }
        }
    }
    colnames(dfAligned) = c(paste0("ref_", names(ref)), paste0("comp_", names(comp)))
    nAligned = dim(dfAligned)[1]
    cat(paste0("    ", nAligned, " features are aligned between runs\n\n"))
    cat(paste0("    ", nAligned, " features are aligned between runs\n\n"), file = LOG)
    
    
    
    ## For public release, "rescue" function is silenced (too complidated)
    
    
        
    ######################################
    ## Rescue some "unaligned" features ##
    ######################################
    # if (rescue == 1) {
    #     rtTolUnit = as.numeric(unlist(strsplit(params$"rt_tolerance_unit", ",")))
    #     mzTolUnit = as.numeric(unlist(strsplit(params$"mz_tolerance_unit", ",")))
    #     rtTol = as.numeric(unlist(strsplit(params$"rt_tolerance_value", ",")))
    #     mzTol = as.numeric(unlist(strsplit(params$"mz_tolerance_value", ",")))
    #     for (i in 1:length(rtTolUnit)) {
    #         res = rescueFeatures(ref, comp, dfAligned, rtSd, mzSd, rtTolUnit[i], rtTol[i], mzTolUnit[i], mzTol[i])
    #         ref = res$ref
    #         comp = res$comp
    #         dfAligned = res$dfAligned
    #     }
    # }
    return (dfAligned)
}

rescueFeatures = function (ref, comp, dfAligned, rtSd, mzSd, rtTolScale, rtTol, mzTolScale, mzTol) {
    nUnaligned = dim(comp)[1]
    cat(paste0("    There are ", nUnaligned, " unaligned features\n"))
    cat("    The features satisfying the following conditions are going to be reviewed for further alignment\n")
    cat(paste0("    There are ", nUnaligned, " unaligned features\n"), file = LOG)
    cat("    The features satisfying the following conditions are going to be reviewed for further alignment\n", file = LOG)
    
    ## 1. Reduce the "reference" run by selecting "unaligned" features
    ##    rtTol and mzTol should be also reduced
    nRescue = 0
    indUnaligned = setdiff(ref$index, dfAligned$"ref_index")
    ind = which(ref$index %in% indUnaligned)
    ref = ref[ind, ]
    rtSd = rtSd[ind]
    mzSd = mzSd[ind]
    
    ## 2. Apply the selection criteria
    ##    - Absolute intensity level: hard-coded (grand median intensity of aligned features)
    ##    - Intensity-ratio between the ref- and comp-runs: hard-coded, within 95% of the ratios of aligned features)
    ##    - RT- and m/z-shifts should be within specified tolerances (e.g. 10SD or 10ppm)
    resIntLevel = median(as.matrix(dfAligned[, grep("intensity", colnames(dfAligned))]))
    resRatioPct = 95 ## Hard-coded
    intRatioAligned = log(dfAligned$"comp_intensity", 2) - log(dfAligned$"ref_intensity", 2)
    lRatio = quantile(intRatioAligned, (1 - resRatioPct / 100) / 2)
    uRatio = quantile(intRatioAligned, 1 - (1 - resRatioPct / 100) / 2)
    cat(paste0("    - Intensity higher than ", resIntLevel, " (median intensity of aligned features)\n"))
    cat(paste0("    - Ratio of intensities within ", resRatioPct, "% of ratios from aligned features\n"))
    cat(paste0("    - Intensity higher than ", resIntLevel, " (median intensity of aligned features)\n"), file = LOG)
    cat(paste0("    - Ratio of intensities within ", resRatioPct, "% of ratios from aligned features\n"), file = LOG)
    if (rtTolScale == 1) { # times of SD unit
        cat(paste0("    - RT-shifts within ", rtTol, " x SD of estimated RT-shifts from aligned features\n"))
        cat(paste0("    - RT-shifts within ", rtTol, " x SD of estimated RT-shifts from aligned features\n"), file = LOG)
        rtTol = rtTol * rtSd
    } else if (rtTolScale == 2) { # second unit
        cat(paste0("    - RT-shifts less than ", rtTol, " seconds\n"))
        cat(paste0("    - RT-shifts less than ", rtTol, " seconds\n"), file = LOG)
        rtTol = rep(rtTol, length(ind))
    } else {
        cat("  - WARNING: check your parameter for RT-tolerance unit. It should be either 1 or 2\n")
        cat("  - Due to the incorrect RT-tolerance unit parameter, the rescue step is skipped\n")
        cat("  - WARNING: check your parameter for RT-tolerance unit. It should be either 1 or 2\n", file = LOG)
        cat("  - Due to the incorrect RT-tolerance unit parameter, the rescue step is skipped\n", file = LOG)
        return (dfAligned)
    }
    if (mzTolScale == 1) { # times of SD unit
        cat(paste0("    - m/z-shifts within ", mzTol, " x SD of estimated m/z-shifts from aligned features\n"))
        cat(paste0("    - m/z-shifts within ", mzTol, " x SD of estimated m/z-shifts from aligned features\n"), file = LOG)
        mzTol = mzTol * mzSd
    } else if (mzTolScale == 2) { # ppm unit
        cat(paste0("    - m/z-shifts less than ", mzTol, " ppm\n"))
        cat(paste0("    - m/z-shifts less than ", mzTol, " ppm\n"), file = LOG)
        mzTol = rep(mzTol, length(ind))
    } else {
        cat("  - WARNING: check your parameter for m/z-tolerance unit. It should be either 1 or 2\n")
        cat("  - Due to the incorrect m/z-tolerance unit parameter, the rescue step is skipped\n")
        cat("  - WARNING: check your parameter for m/z-tolerance unit. It should be either 1 or 2\n", file = LOG)
        cat("  - Due to the incorrect m/z-tolerance unit parameter, the rescue step is skipped\n", file = LOG)
        return (dfAligned)
    }
    
    for (i in 1:dim(ref)[1]) {
        intRatio = log(comp$intensity, 2) - log(ref$intensity[i], 2)
        mzShift = (comp$mz - ref$mz[i]) / comp$mz * 1e6
        rtShift = comp$RT - ref$RT[i]
        selInd = which(comp$intensity > resIntLevel &
                           intRatio >= lRatio & intRatio <= uRatio &
                           abs(mzShift) <= mzTol[i] &
                           abs(rtShift) <= rtTol[i])
        if (length(selInd) > 0) {
            ## Check charge states of the matched features
            charge = ref$z[i]
            if (charge == 0) {
                selInd = selInd[1]
                rescued = cbind(ref[i, ], comp[selInd, ])
                names(rescued) = names(dfAligned)
                dfAligned = rbind(dfAligned, rescued)
                comp = comp[-selInd, ]
                nRescue = nRescue + 1
            } else {
                if (comp$z[selInd[1]] == 0) {
                    selInd = selInd[1]
                    rescued = cbind(ref[i, ], comp[selInd, ])
                    names(rescued) = names(dfAligned)
                    dfAligned = rbind(dfAligned, rescued)
                    comp = comp[-selInd, ]
                    nRescue = nRescue + 1
                } else {
                    chargeInd = which(comp$z[selInd] == charge)
                    if (length(chargeInd) > 0) {
                        selInd = selInd[chargeInd[1]]
                        rescued = cbind(ref[i, ], comp[selInd, ])
                        names(rescued) = names(dfAligned)
                        dfAligned = rbind(dfAligned, rescued)
                        comp = comp[-selInd, ]
                        nRescue = nRescue + 1
                    }
                }
            }
        } ## end if (length(selInd) > 0) 
    } ## end for
    cat(paste0("    Through the rescue procedure, ", nRescue, " features are additionally aligned\n\n"))
    cat(paste0("    Through the rescue procedure, ", nRescue, " features are additionally aligned\n\n"), file = LOG)
    res = list(ref = ref, comp = comp, dfAligned = dfAligned)
    return (res)
}

findMatchedFeatures = function(matched, calibrated, files, samples) {
    ## Caution
    ## Input arguments files and samples should have the same order as defined in the parameter file
    nFiles = length(matched)
    files = files[1:nFiles]
    
    refRunNo = NULL
    for (i in 1:length(matched)) {
        if (is.null(dim(matched[[i]]))) {
            refRunNo = i
        }
    }
    
    ## Extract the reference run's indexes of fully- and partially matched features
    fullInd = NULL
    partialInd = NULL
    for (i in 1:length(matched)) {
        if (i == refRunNo) {
            next
        } else {
            if (is.null(fullInd)) {
                fullInd = matched[[i]][, 1]
                partialInd = matched[[i]][, 1]
            } else {
                fullInd = intersect(fullInd, matched[[i]][, 1])
                partialInd = union(partialInd, matched[[i]][, 1])
            }
        }
    }
    fullInd = sort(unique(fullInd))
    partialInd = setdiff(partialInd, fullInd)
    partialInd = sort(unique(partialInd))
    
    ###############################
    ## 1. Fully-matched features ##
    ###############################
    features = list()
    runNo = setdiff(seq(1, length(matched)), refRunNo)[1]
    nCols = ncol(matched[[runNo]]) / 2
    for (i in 1:length(matched)) {
        if (i == refRunNo) {
            tmp = matched[[runNo]][order(matched[[runNo]][, 1]), ]
            tmp = tmp[tmp[, 1] %in% fullInd, ]
            features = c(features, list(tmp[, 1:nCols]))
        } else {
            tmp = matched[[i]][order(matched[[i]][, 1]), ]
            tmp = tmp[tmp[, 1] %in% fullInd, ]
            features = c(features, list(tmp[, (nCols + 1):(2 * nCols)]))
        }
    }
    features = as.data.frame(features, check.names = F)
    
    ## Manipulation of column names
    header = NULL
    for (i in 1:nFiles) {
        filename = basename(files[i])
        filename = gsub(".\\d+.feature", "", filename)
        header = c(header, rep(filename, nCols))
    }
    colnames(features) = gsub("ref_", "", colnames(features))
    colnames(features) = gsub("comp_", "", colnames(features))
    colnames(features) = paste0(header, "_", colnames(features))
    fullFeatures = features
    
    ####################################
    ## 2. Partially-matched features  ##
    ####################################
    features = data.frame(matrix(NA, nrow = length(partialInd), ncol = nCols * nFiles))
    for (i in 1:length(partialInd)) {
        for (j in 1:length(matched)) {
            if (j == refRunNo) {
                next
            } else {
                ind = which(matched[[j]][, 1] == partialInd[i])
                if (length(ind) == 1) {
                    features[i, ((j - 1) * nCols + 1) : (j * nCols)] = matched[[j]][ind, (nCols + 1) : (2 * nCols)]
                    features[i, ((refRunNo - 1) * nCols + 1) : (refRunNo * nCols)] = matched[[j]][ind, 1 : nCols]
                } else if (length(ind) > 1) {
                    stop ('Error in partially-matched feature identification')
                }
            }
        }
    }
    
    ## Manipulation of column names
    header = NULL
    for (i in 1:nFiles) {
        filename = basename(files[i])
        filename = gsub(".\\d+.feature", "", filename)
        header = c(header, rep(filename, nCols))
    }
    colnames(features) = rep(colnames(calibrated[[1]]), nFiles)
    colnames(features) = paste0(header, "_", colnames(features))
    partialFeatures = features
    
    ############################
    ## 3. Un-matched features ##
    ############################
    unmatchedFeatures = list()
    for (i in 1:nFiles) {
        filename = basename(files[i])
        filename = gsub(".\\d+.feature", "", filename)
        if (i == refRunNo) {
            unInd = unique(setdiff(calibrated[[i]]$index, union(fullInd, partialInd)))
            # unmatchedFeatures[[fileName]] = calibrated[[i]][calibrated[[i]]$index %in% unInd, ]
        } else {
            unInd = setdiff(calibrated[[i]]$index, matched[[i]][, (nCols + 1)])
            # unmatchedFeatures[[fileName]] = calibrated[[i]][calibrated[[i]]$index %in% unInd, ]
        }
        unmatchedFeatures[[filename]] = calibrated[[i]][calibrated[[i]]$index %in% unInd, ]
    }
    
    ##############################
    ## Alignment/match summary  ##
    ##############################
    cat("  Alignment/matching summary\n")
    cat("    After alignment/feature matching, fully-, partially- and un-aligned features are as follows\n")
    cat("    Filename\t\tFully-aligned\tPartially-aligned\tun-aligned\n")
    cat("  Alignment/matching summary\n", file = LOG)
    cat("    After alignment/feature matching, fully-, partially- and un-aligned features are as follows\n", file = LOG)
    cat("    Filename\t\tFully-aligned\tPartially-aligned\tun-aligned\n", file = LOG)
    for (i in 1:nFiles) {
        filename = basename(files[i])
        filename = gsub(".\\d+.feature", "", filename)
        nFull = length(fullInd)
        nPartial = sum(!is.na(partialFeatures[[paste0(filename, "_index")]]))
        nUnmatched = dim(unmatchedFeatures[[filename]])[1]
        cat(paste0("    ", filename, "\t\t", nFull, "\t\t", nPartial, "\t\t", nUnmatched, "\n"))
        cat(paste0("    ", filename, "\t\t", nFull, "\t\t", nPartial, "\t\t", nUnmatched, "\n"), file = LOG)
    }
    cat("\n")
    cat(paste0("    There are ", nrow(partialFeatures), " partially-aligned features\n"))
    cat("\n", file = LOG)
    cat(paste0("    There are ", nrow(partialFeatures), " partially-aligned features\n"), file = LOG)
    colInd = grep("index", colnames(partialFeatures))
    nRuns = apply(!is.na(partialFeatures[, colInd]), 1, sum)
    for (i in (length(files) - 1):2) {
        cat(paste0("    Partially-aligned features over ", i, " runs = ", sum(nRuns == i), "\n"))
        cat(paste0("    Partially-aligned features over ", i, " runs = ", sum(nRuns == i), "\n"), file = LOG)
    }
    
    ## Depending on the parameter, some "partially-aligned" features are merged to "fully-aligned" ones
    pctFullyAlignment = params$"pct_full_alignment"
    if (pctFullyAlignment < 100) {
        colInd = grep("index", colnames(partialFeatures))
        nRuns = apply(!is.na(partialFeatures[, colInd]), 1, sum)
        rowInd = which(nRuns >= ceiling(pctFullyAlignment / 100 * length(files)))	  
        if (length(rowInd) > 0) {
            fullFeatures = rbind(fullFeatures, partialFeatures[rowInd, ])
            fullInd = c(fullInd, partialInd[rowInd])
            partialFeatures = partialFeatures[-rowInd, ]
            partialInd = partialInd[-rowInd]
            cat("\n")
            cat("    According to the parameter setting,", length(rowInd), "partially-aligned features are regarded as fully-aligned\n")
            cat("\n", file = LOG)
            cat("    According to the parameter setting,", length(rowInd), "partially-aligned features are regarded as fully-aligned\n", file = LOG)
        } else {
            cat("\n")
            cat("    According to the parameter setting, no feature is added to the set of fully-aligned ones\n")
            cat("\n", file = LOG)
            cat("    According to the parameter setting, no feature is added to the set of fully-aligned ones\n", file = LOG)
        }
    }
    
    cat("\n")
    cat("\n", file = LOG)

    ############
    ## Output ##
    ############
    res = list(fullInd = fullInd, fullFeatures = fullFeatures, 
               partialInd = partialInd, partialFeatures = partialFeatures,
               unmatchedFeatures = unmatchedFeatures)
    return (res)
}

quantifyFeatures = function(fullFeatures, params, LOG) {
    expr = fullFeatures[, grep("intensity", colnames(fullFeatures))]
    
    ## Calculation and printout of loading-biases (filter out the most variable features +/- 10%)
    cat("  Loading-bias summary\n")
    cat("  Loading-bias summary\n", file = LOG)
    colInd = grep("intensity", colnames(expr))
    sampleName = colnames(expr)[colInd]
    sampleName = gsub("_intensity", "", sampleName)
    lexpr = expr[, colInd] ## Feature intensity table for loading-bias calculation
    nFeatures = dim(lexpr)[1]
    nSamples = dim(lexpr)[2]
    sampleMeans = rowMeans(lexpr, na.rm = T)
    lexpr = log2(lexpr / sampleMeans)
    rowInd = seq(1, nFeatures, by = 1)
    for (i in 1:nSamples) {
        rowInd = intersect(rowInd, which(lexpr[, i] < quantile(lexpr[, i], 0.9, na.rm = T) & lexpr[, i] > quantile(lexpr[, i], 0.1, na.rm = T)))
    }
    meanIntensity = round(as.numeric(2 ^ colMeans(lexpr[rowInd, ]) * 100), 2)
    sdVal = as.numeric(apply(lexpr[rowInd, ], 2, sd))
    sdIntensity = round(as.numeric(((2 ^ sdVal - 1) + (1 - 2 ^ (-sdVal))) / 2 * 100), 2)
    semIntensity = round(as.numeric(sdIntensity / sqrt(length(rowInd))), 2)
    cat("    Samplename\tMean[%]\tSD[%]\tSEM[%]\t#features\n")
    cat("    Samplename\tMean[%]\tSD[%]\tSEM[%]\t#features\n", file = LOG)
    for (i in 1:nSamples) {
        cat("   ", sampleName[i], "\t", meanIntensity[i], "\t", sdIntensity[i], "\t", semIntensity[i], "\t", length(rowInd), "\n")
        cat("   ", sampleName[i], "\t", meanIntensity[i], "\t", sdIntensity[i], "\t", semIntensity[i], "\t", length(rowInd), "\n", file = LOG)
    }
    
    ## Normalization (i.e. loading-bias correction) based on trimmed-mean intensities
    expr = log2(expr)
    rowSel = NULL
    if (params$"skip_loading_bias_correction" == "0") {
        ## Parameters for normalization
        intensityThreshold_normalization = quantile(as.matrix(expr), 0.2, na.rm = T) ## 10% percentile of overall intensities
        trimPct_normalization = 0.1 ## 10% lowest and 10% highest intensities are to be trimmed
        rowSel = expr > intensityThreshold_normalization ## Pre-filtering based on the intensity level
        for (i in 1:ncol(expr)) {
            rowSel[, i] = (expr[, i] > quantile(expr[, i], trimPct_normalization, na.rm = T) &
                               expr[, i] < quantile(expr[, i], 1 - trimPct_normalization, na.rm = T))
        }
        rowInd = which(apply(rowSel, 1, sum) == ncol(rowSel))
        meanIntensity = colMeans(expr[rowInd, ])
        normFactor = meanIntensity - mean(meanIntensity)
        expr = expr - rep(normFactor, each = nrow(expr))
    }
    
    ## Replace NAs in expr with the half of grand minimum intensity
    rowInd = which(apply(is.na(expr), 1, sum) > 0)
    if (length(rowInd) > 0) {
        grandMin = min(expr, na.rm = T)
        for (i in 1:length(rowInd)) {
            colInd = which(is.na(expr[rowInd[i], ]))
            rMean = mean(as.numeric(expr[rowInd[i], ]), na.rm = T)
            expr[rowInd[i], colInd] = (grandMin - 1) + runif(length(colInd)) * 1e-3 ## For numerical stability, small random number is added
        }    
    }
    
    colInd = grep("intensity", colnames(fullFeatures))
    for (i in 1:length(colInd)) {
        ## Replace intensity values in res$fullFeatures with the normalized intensity values
        fullFeatures[, colInd[i]] = 2 ^ expr[, i]
    }
    
    return (fullFeatures)
}


#################
##### Main ######
#################
args = commandArgs(trailingOnly = TRUE)

# For test,
args[1] = "U:/Research/Projects/7Metabolomics/JUMPm/IROAsamples/IROA_IS_NEG_1.1.dev.feature"
args[2] = "U:/Research/Projects/7Metabolomics/JUMPm/IROAsamples/IROA_IS_NEG_2.1.dev.feature"
args[3] = "U:/Research/Projects/7Metabolomics/JUMPm/IROAsamples/IROA_IS_NEG_3.1.dev.feature"
args[4] = "../jumpm_negative_desktop.params"
args[5] = nullfile()

##################################
## Load features and parameters ##
##################################
nArgs = length(args)


## To-do: argument check function


nFeatureFiles = nArgs - 2 ## The last argument is a parameter file
features = list()
for (i in 1:nFeatureFiles) {
    df = read.table(args[i], header = T, sep = "\t", row.names = NULL,
                    stringsAsFactors = F, comment.char = "", check.names = F)
    colnames(df) = gsub(" ", "", colnames(df))
    colnames(df) = gsub("/", "", colnames(df))
    features[[i]] = df

}
paramFile = args[nArgs - 1]
params = getParams(paramFile)
LOG = args[nArgs]

# Assume that "features[[i]]" is a data.frame with the following column names
# 'index' = nominal feature index
# 'mz' = m/z value of a feature
# 'z' = charge of the feature (0 for undetermined)
# 'MS1ScanNumber' = representative MS1 scan number of the feature
# 'minMS1ScanNumber' = minimum MS1 scan number of the feature (spanned)
# 'maxMS1ScanNumber' = maximum MS1 scan number of the feature (spanned)
# 'RT' = representative RT (should correspond to MS1ScanNumber)
# 'minRT' = minimum RT
# 'maxRT' = maximum RT
# 'intensity' = intensity of the feature
# 'SN' = signal-to-noise ratio
# 'PercentageTF' = percentage of true feature (= feature width over RT/scan)

############################################
## Feature calibration (global and local) ##
############################################
calibratedFeatures = list()
if (length(features) > 1) {
    cat("\n\n  Feature calibration\n")
    cat("  ===================\n")
    cat("\n\n  Feature calibration\n", file = LOG)
    cat("  ===================\n", file = LOG)

    # Select a reference run
    if (params$"reference_feature" == "0") {
        ## A run whose median intensity of top 100 features is set to a reference
        refNo = 1
        refIntensity = 0
        for (i in 1:length(features)) {
            tmpIntensity = median(sort(features[[i]]$intensity, decreasing = T)[1:100])
            if (tmpIntensity >= refIntensity) {
                refNo = i
                refIntensity = tmpIntensity
            }
        }
    } else {
        refNo = which(args == params$"referene_feature")
    }

    ## Calibration of features
    rtSd = list()
    mzSd = list()
    for (i in 1:length(features)) {
        if (i == refNo) {
            calibratedFeatures[[i]] = features[[i]]
            rtSd[[i]] = NA
            mzSd[[i]] = NA
        } else {

            cat(paste0("  ", basename(args[i]), " is being aligned against the reference run (it may take a while)\n"))
            cat(paste0("  ", basename(args[i]), " is being aligned against the reference run (it may take a while)\n"), file = LOG)
            res = calibrateFeatures(features[[refNo]], features[[i]], params, LOG)
            calibratedFeatures[[i]] = res$calibratedFeature
            rtSd[[i]] = res$rtSd
            mzSd[[i]] = res$mzSd
        }
    }
    cat("\n")
    cat("\n", file = LOG)
} else {
    cat("Since a single file is used, the feature alignment is skipped\n")
}

################################
## Summary of the calibration ##
################################
cat("  Calibration summary\n")
cat("    After calibration, RT- and m/z-shifts of each run (against a reference run) are centered to zero\n")
cat("    Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows\n")
cat("    Filename\t\t# features\tSD RT-shift[second]\tSD m/z-shift[ppm]\n")
cat("  Calibration summary\n", file = LOG)
cat("    After calibration, RT- and m/z-shifts of each run (against a reference run) are centered to zero\n", file = LOG)
cat("    Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows\n", file = LOG)
cat("    Filename\t\t# features\tSD RT-shift[second]\tSD m/z-shift[ppm]\n", file = LOG)
for (i in 1:nFeatureFiles) {
    filename = basename(args[i])
    filename = gsub(".\\d+.feature", "", filename)
    nFeatures = nrow(features[[i]])
    if (i == refNo) {
        meanRtSd = "NA (reference)"
        meanMzSd = "NA (reference)"
    } else {
        meanRtSd = sprintf("%.6f", mean(rtSd[[i]]))
        meanMzSd = sprintf("%.6f", mean(mzSd[[i]]))
    }
    cat(paste0("    ", filename, "\t\t", nFeatures, "\t\t", meanRtSd, "\t", meanMzSd, "\n"))
    cat(paste0("    ", filename, "\t\t", nFeatures, "\t\t", meanRtSd, "\t", meanMzSd, "\n"), file = LOG)
}
cat("\n")
cat("\n", file = LOG)

###################################################################
## Identification of fully-aligned features for further analysis ##
###################################################################
cat("  Feature alignment\n")
cat("  =================\n")
cat("  Feature alignment\n", file = LOG)
cat("  =================\n", file = LOG)
matchedFeatures = list()
for (i in 1:nFeatureFiles) {
    if (i == refNo) {
        matchedFeatures[[i]] = NA
    } else {
        matchedFeatures[[i]] = matchFeatures(calibratedFeatures[[refNo]], calibratedFeatures[[i]],
                                             basename(args[refNo]), basename(args[i]),
                                             rtSd[[i]], mzSd[[i]], params, LOG)
    }
}

## Fully-, partially- and un-matched/aligned features
res = findMatchedFeatures(matchedFeatures, calibratedFeatures, args)


##############################
## Processing quantity data ##
##############################
cat("  Quantity information\n")
cat("  ====================\n")
cat("  Quantity information\n", file = LOG)
cat("  ====================\n", file = LOG)
fullFeatures = res$fullFeatures
fullFeatures = quantifyFeatures(fullFeatures, params, LOG)



