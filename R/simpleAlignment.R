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

args = commandArgs(trailingOnly = TRUE)
# args[1] = "../../Dev/JUMPm_unlabel_python/refRt.txt"
# args[2] = "../../Dev/JUMPm_unlabel_python/compRt.txt"
# args[3] = "../../Dev/JUMPm_unlabel_python/compRt_new.txt"
rts = read.table(args[1], sep = "\t", row.names = NULL, stringsAsFactors = F, comment.char = "", check.names = F)
rts = as.matrix(rts)
rowInd = which(!is.na(rts[,2]))
refRt = rts[rowInd, 2]
compRt = rts[rowInd, 1]
compRt_new = rts[, 1]
x = as.numeric(compRt)
y = as.numeric(compRt) - as.numeric(refRt) ## RT-shifts

## Removal of outliers in y (i.e. RT-shift)
pct = 0.2
truncatedMean = mean(y[y >= quantile(y, pct / 2) & y <= quantile(y, (1 - pct / 2))])
truncatedSd = sd(y[y >= quantile(y, pct / 2) & y <= quantile(y, (1 - pct / 2))])
lL = truncatedMean - 3 * truncatedSd
uL = truncatedMean + 3 * truncatedSd
ind = which(y >= lL & y <= uL)

mod = loess.as(x[ind], y[ind], degree = 1, criterion = "aicc",
               control = loess.control(surface = "direct")) ## This curve represents RT-shifts as a function of comp$Rt
compRt_new = as.numeric(compRt_new) - predict(mod, data.frame(x = as.numeric(compRt_new)))
write.table(file = "alignment_result.txt", compRt_new, sep = "\t", row.names = F, col.names = F, quote = F)
write.table(file = "alignment_residual.txt", abs(mod$residuals), sep = "\t", row.names = F, col.names = F, quote = F)
