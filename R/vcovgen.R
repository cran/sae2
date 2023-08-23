vcovgen <-
function(linear.subs, year.list, geocode, designvars)
{
   var.list <- names(linear.subs)
   var.list <- var.list[!var.list %in% c("Year", geocode, designvars)]
   geo.list <- sort(unique(linear.subs[, geocode]))
   var.year <- NULL
   for (i in 1:length(var.list)) {
       var <- var.list[i]
       var.year <- append(var.year, paste0(var, ".", year.list))
   }
   vcovlist <- list()
   indx <- 0
   for (geo in geo.list) {
      subs <- linear.subs[linear.subs$Year %in% year.list & 
                          linear.subs[, geocode] == geo, ]
      X <- matrix(0, nrow=nrow(subs), 
                     ncol=length(year.list) * length(var.list))
      for (i in 1:length(var.list)) {
         var <- var.list[i]
         for (j in 1:length(year.list)) {
            k <- (i-1) * length(year.list) + j
            year <- year.list[j]
            X[subs$Year == year, k] <- subs[subs$Year == year, var]
         }
      }
      dimnames(X)[[2]] <- as.list(var.year)
      X <- data.frame(X)
      X[, designvars] <- subs[, designvars]
      ids.f <- formula(paste("~", designvars[2]))
      str.f <- formula(paste("~", designvars[1]))
      dsgn <- svydesign(ids=ids.f, 
                     strata=str.f,
                     data=X,
                     weights=~1, nest=TRUE)
      formula1 <- make.formula(var.year)
      indx <- indx + 1
      vcovlist[[indx]] <- vcov(svytotal(formula1, design=dsgn))
      if(!inherits(vcovlist[[indx]], "matrix")) {
         vcovlist[[indx]] <- matrix(as.numeric(vcovlist[[indx]]), nrow=1, ncol=1)
      } 
   }
   names(vcovlist) <- geo.list
   return(vcovlist)      
}
