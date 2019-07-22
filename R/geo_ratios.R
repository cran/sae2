geo_ratios <-
function(data, geocode, numerators, denominators, geonames, new.names,
         designvars) 
{
    if (missing(designvars)) {
        full.list <- unique(c(geocode, numerators, 
                       denominators, "Year"))
    } else {
        full.list <- unique(c(geocode, numerators, 
                       denominators, designvars, "Year"))
    }   
    if( any(!full.list %in% names(data)) ) {
        not.found <- full.list[!full.list %in% names(data)]
        stop(paste0("  Variables not in data frame: \n", 
               paste0(not.found, collapse=", ")))
    }
    if(length(unique(numerators)) != length(numerators))
            stop("Numerators must not contain duplicates")
    if(!missing(new.names)) {
        if(length(new.names) != length(numerators)) 
            stop(paste("Length of new.names =", length(new.names),
                 "must agree with length of numerators =", 
                 length(numerators)))
        if(length(unique(new.names)) != length(new.names))
            stop("new.names must not contain duplicates") 
    }
    if(length(numerators) > 1 & length(denominators) == 1) 
        denominators <- rep(denominators, length(numerators))
    data <- data[, full.list]
    var.list <- unique(c(numerators, denominators))
    if (length(var.list) > 1) {
        formula <- paste0("cbind(", paste0(var.list, collapse=", "), ")")
    } else {
        formula <- var.list
    }
    formula <- formula(paste(formula, "~ Year +", geocode))
    xx <- aggregate(formula, data=data, FUN=sum)
#
# different syntax for 1 vs 2+ numerators
#
    if(length(numerators) == 1) {
        result <- c(xx[, numerators])/xx[, denominators]
        result <- data.frame(result)
        if(!missing(new.names)) {
            names(result) <- new.names
        } else {
            names(result) <- numerators
            new.names <- numerators
        }
    } else {
        result <- as.matrix(xx[, numerators]) 
        for (i in 1:length(numerators))
            result[, i] <- result[, i]/xx[, denominators[i]]
        if(!missing(new.names)) {
            dimnames(result)[[2]] <- new.names
        } else {
            dimnames(result)[[2]] <- numerators
            new.names <- numerators
        }
        result <- data.frame(result)
    }
    result[, geocode] <- xx[, geocode]
    result$Year <- xx$Year
    if(!missing(geonames)) 
        result <- merge(geonames, result, by=geocode, all.x=TRUE)
    result <- result[order(result[, geocode], result$Year), ]
    if(!missing(designvars)) {
        result.2 <- result
        added.name <- paste0(full.list[which.max(nchar(full.list))], ".temp")
        data[, added.name] <- paste0(data$Year, data[, geocode])
        result.2[, added.name] <- paste0(result.2$Year, result.2[, geocode])
# 
# result.2 is based on result but includes added.name for matching
#
        denom.list <- unique(denominators)
        if (length(denom.list) > 1) {
            formula <- paste0("cbind(", paste0(denom.list, collapse=", "), ")")
        } else {
            formula <- denom.list
        }
        formula <- formula(paste(formula, "~", added.name))
        result.3 <- aggregate(formula, data=data, FUN=sum)
#
# result.3 has denominator totals by added.name
#
        designvars <- designvars[!designvars %in% c("Year", geocode)]
        if (length(designvars) == 0)
           stop(paste(" designvars must contain design information beyond",
                      "year and geocode"))
        new.list <- unique(c(numerators, denominators))
        if (length(new.list) > 1) {
            formula <- paste0("cbind(", paste0(new.list, collapse=", "), ")")
        } else {
            formula <- new.list
        }
        formula <- formula(paste(formula, "~ ", 
                      paste0(c(added.name, "Year", geocode, designvars), 
                               collapse=" + ")))
        xx <- aggregate(formula, data=data, FUN=sum)
#
# xx has numerator and denominator totals by added name, year, geocode, and
# design variables
#
        z <- matrix(nrow=dim(xx)[1], ncol=length(new.names))
        for (i in 1:length(new.names)) {           
            z[, i] <- c(( xx[, numerators[i]] - 
                result.2[match(xx[, added.name], 
                               result.2[, added.name]), new.names[i]] *
                    xx[, denominators[i]]) /
                result.3[match(xx[, added.name], 
                               result.3[, added.name]), denominators[i]])
        }   
        z <- data.frame(z)
        names(z) <- new.names
        z <- cbind(z, xx[, c("Year", geocode, designvars)])
        return(list(estimates=result, linear.subs=z))
    } else {
        return(list(estimates=result))
    }
}

