glme = function(fixed,
                data = sys.frame(sys.parent()),
                random = pdSymm(eval(as.call(fixed[-2]))),
                correlation = NULL,
                weights = NULL,
                subset,
                method = "GM",
                na.action = na.fail,
                control = list(),
                contrasts = NULL,
                keep.data = TRUE)
{
  data_name <- deparse(substitute(data))
  Call <- X2 <- X1 <- value <- NULL
  # libraries used
  # library(nlme)
  # library(stringr)
  # library(reshape)
  # library(Matrix)
  # library(tidyverse)
  # library(dplyr)
  REML <- method == "REML"
  data$ones <- 1
  # ML and REML
  if(method == "REML"| method == "ML"){
    # Code below is for GM only; otherwise just use lme
    out <- nlme::lme(fixed,
                     data,
                     random,
                     correlation,
                     weights,
                     subset,
                     method,
                     na.action,
                     control,
                     contrasts,
                     keep.data)
    return(out)
  }
  # Proceed to GM if not ML/REML;
  # Estimate Fized Effects by GLSE and Random Effects one at a time
  GM <- method == "GM"
  data <- data.frame(data)
  N <- nrow(data)
  data$Intercept <- rep(1, N)
  reSt <- reStruct(random, FALSE, data = NULL)
  lmeSt <- lmeStruct(reStruct = reSt, corStruct = correlation,
                     varStruct = varFunc(weights))
  groups <- getGroupsFormula(reSt)
  mfArgs <- list(formula = asOneFormula(formula(lmeSt), fixed, groups),
                 data = data, na.action = na.action)
  if (!missing(subset)) mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2L]]
  mfArgs$drop.unused.levels <- TRUE
  dataMix <- do.call(model.frame, mfArgs)
  origOrder <- row.names(dataMix)
  # Work with Fixed Effects and Random Effects as detected by lme

  X <- model.frame(fixed, dataMix)
  X <- unlist(model.matrix(fixed, data = X)) # Fixed part of mixed model design matrix
  y <- eval(fixed[[2L]], dataMix)

  temp <- as.character(random)
  ztemp <- strsplit(temp, ",")
  ztmp <- ztemp
  if (ztemp[[1]] == "~") nrands <- 1 else nrands <- length(ztemp)
  Gvars <- rep("", nrands)
  Zvars <- rep("", nrands)

  for (r in 1:nrands){
    tem <- strsplit(temp,"~")[[r]]
    if (nrands > 1)
      tem <- tem[tem != ""]
    else tem <- ztemp[[2]]
    tem2 <- strsplit(tem, "[|]")

    gtmp <- tem2[[1]][2]
    ztmp <- tem2[[1]][1]
    Gvars[r] <- trimws(gtmp) # Drop spaces
    ztmp <- trimws(ztmp)
    if (ztmp == "1") ztmp <- "Intercept"
    Zvars[r] <- ztmp
  }

  # Handle with and with Intercept cases
  Zonly <- Zvars
  Intercepts <- rep("Intercept", nrands)

  for (r in 1:nrands){
    tem <- strsplit(Zvars[r], "[+]")[[1]]
    if (length(tem) > 1) {
      Zonly[r] <- trimws(tem[2])
      if(trimws(tem[1])	== "-1") Intercepts[r] <- ""
    }
  }
  zonly <- Zonly


  # Create Dummy matrix to track active variales and create Zonly
  IntMat <- matrix(0, nrands, 2)
  # Matrix to track random effects with or without intercepts
  colnames(IntMat) <- c("Intercept", "Z")
  rownames(IntMat) <- Gvars
  IntMat[Intercepts == "Intercept", 1] <- 1
  IntMat[Intercepts == "Intercept", 1] <- 1
  IntMat[Zonly == "Intercept", 1] <- 1
  IntMat[Zonly == "Intercept", 2] <- 0
  IntMat[Zonly != "Intercept", 2] <- 1
  intmat <- IntMat

  # Also get y and X variable name
  temp <- as.character(fixed)
  tmp <- strsplit(temp, "~")
  yname <- tmp[[2]][1]
  tmpp <- strsplit(tmp[[3]][1], "[+]")
  xtmp <- strsplit(tmpp[[1]], " /") # drop spaces
  Xnames <- NULL
  for (i in 1:length(xtmp)){
    tmp = unlist(xtmp[[i]])
    tmp = unlist(strsplit(tmp," "))
    if (length(tmp) > 1){
      if(tmp[1] == "") tmp <- tmp[2]
      else tmp <- tmp[1]
    }
    Xnames = c(Xnames, tmp)
  }

  # Run LSE first
  zPart <- NULL
  for (r in 1:nrands){
    if (IntMat[r, 1] == 1){
      itmp <- paste(Intercepts[r], ":", Gvars[r], sep = "")
      zPart <- paste(zPart, "+", itmp, sep = "") # factor part of formula
    }
    if (IntMat[r, 2] == 1){
      ztmp <- paste(Zonly[r], ":", Gvars[r], sep = "")
      zPart <- paste(zPart, "+", ztmp, sep = "") # factor part of formula
    }
  }
  ftmp <- as.character(fixed)
  lmformula <- paste(yname, "~", ftmp[3], zPart)
  rlmXZ <- lm(lmformula, data = data, singular.ok = T, x = T, y = T)

  srlmXZ <- summary(rlmXZ)

  Estimates <- rlmXZ$coef
  Estimates <- Estimates[!is.na(Estimates)]

  estimates <- Estimates


  # Do without intercept  to get more non-NA Intercepts and to compare fo now
  nZs <- length(Zvars)
  zlen <- length(Gvars)

  Allxxinvu <- srlmXZ$cov
  goodcols <- colnames(Allxxinvu)
  if(nZs - length(goodcols) > 5){
    print("Extreme Collinearity. Drop some Random Effects or their intercepts")
    return(0)
  }
  ZEstimates <- Estimates[goodcols]
  e <- rlmXZ$df.resid # U degrees of freedom
  allXZ <- rlmXZ$x
  Znames <- names(ZEstimates)
  Zdat <- allXZ[,Znames]
  Gvars <- unique(Gvars)
  zlen <- length(Gvars)
  ztemp <- list(NULL, TRUE, 1:zlen)
  Intertemp <- list(NULL, 1:zlen)

  for (r in 1:zlen){
    Intertmp <-  Znames[grep("Intercept", Znames)]

    if (IntMat[r, 1] == 1){
      Inter <- Intertmp[grep(Gvars[r], Intertmp)]
      Intertemp[[r]] <- Inter
    }
    # Check adequacy of data
    if (length(Intertemp[[r]]) > 1){
      if (length(Intertemp[[r]]) < 2){
        print("Data is not appropriate even for regular Regression. Make sure grouping variable is a factor. Run RANOVA instead.")
        return(0)
      }
    }
    if (IntMat[r, 2] == 1 ){
      ztmpp <- unique(Zonly)
      zl <- length(ztmpp)
      if (zl == 1)
        ztmp <- Znames[grep(ztmpp, Znames)]
      else ztmp <- Znames[grep(ztmpp[r], Znames)]

      ztemp[[r]] <- ztmp[-1] # zString[-1]
    }
    # Handle case of only intercept
    if (IntMat[r, 1] == 1 & IntMat[r, 2] == 0 ){
      ztemp[[1]] <- Intertemp[[1]]
    }
  } # end of r randpmozed variable
  if (1 == 0) { # For testing ony
    allxxinvu <- Allxxinvu
    zestimates <- ZEstimates
    znames <- Znames
    zdat <- Zdat
    ztemp <- ztemp
    inttmp <- Intertemp
  }

  # Now separte out estimates for each random effect
  uval1 <- list(NULL, TRUE, 1:zlen) # Intercepts, formerly IntrandEsts
  uval2 <- list(NULL, TRUE, 1:zlen) # random effects, randEsts
  Zfixed = matrix(0, zlen, 2)
  colnames(Zfixed) <- c("Intercept", "Z")
  for (r in 1:zlen) {
    # Handle Intercepts
    if (IntMat[r, 1] == 1){
      tmp1 <- ZEstimates[Intertemp[[r]]]
      Zfixed[r, 1] <- mean(tmp1, na.rm = T)
      uval1[[r]] <- tmp1 - Zfixed[r, 1]
    }
    else tmp1 <- 0
    # Handle Random Effects
    if (IntMat[r, 2] == 1){
      tmp2 <- ZEstimates[ztemp[[r]]]
      Zfixed[r, 2] <- mean(tmp2, na.rm = T)
      uval2[[r]] <- tmp2 - Zfixed[r, 2]
    }
    else tmp2 <- 0
  }

  rownames(Zfixed) <- Gvars #zVariables

  # Get Allxxinvu
  # Now find utilde for each random effect

  u.tilde1 <- list(NULL, TRUE, 1:zlen) # Intercepts, formerly IntrandEsts
  u.tilde2 <- list(NULL, TRUE, 1:zlen) # for random effects
  lambda1  <- list(NULL, TRUE, 1:zlen)
  lambda2  <- list(NULL, TRUE, 1:zlen)
  xxinvu1  <- list(NULL, TRUE, 1:zlen)
  xxinvu2  <- list(NULL, TRUE, 1:zlen)
  for (r in 1:zlen) {
    # Handle intercepts of random effects
    if(IntMat[r, 1] == 1){
      xxinvu1[[r]] <- Allxxinvu[Intertemp[[r]], Intertemp[[r]]]
      eout1 <- eigen(xxinvu1[[r]])
      lambda1[[r]] <- eout1$val
      P1 <- eout1$vec
      u.tilde1[[r]] <- t(P1) %*% uval1[[r]] # for ssrho or ssalp
    }
    # Handle Random Effects
    if(IntMat[r, 2] == 1){
      xxinvu2[[r]] <- Allxxinvu[ztemp[[r]], ztemp[[r]]]
      eout2 <-eigen(xxinvu2[[r]])
      lambda2[[r]] <- eout2$val
      P2 <- eout2$vec
      u.tilde2[[r]] <- t(P2) %*% uval2[[r]] # us for ssrho or ssalp
    }
  } # end of r in 1:zlen

  SSE <- sum(srlmXZ$resid ^ 2)
  MSE <- SSE / rlmXZ$df

  # Estimate Variance Components by Solving the equation E(V * Se/Usa(siga^2 *U/se >U| Condition)=1
  nz <- 2
  e <- rlmXZ$df.resid # U degrees of freedom
  allx <- rlmXZ$x
  rhos <- matrix(0, zlen, nz)
  alps <- matrix(0, zlen, nz)

  M <- 10000
  zlen1 <- zlen + 1
  rout <- matrix(0, zlen1, 2)
  colnames(rout) <-c("SD of Intercept","SD of Random Effect")
  rownames(rout) <- c(Gvars, "Residual")

  # Create Zu matrix for all random effects;
  # We cannot use kronecker method since frequencies may be diffeent for some grpups
  zeros <- rep(0, N)
  ones <- rep(1, N)
  grpLengths <- rep(0, zlen)


  # Collect Zu matrics, ranks and ng by group;
  # keep Zu by each random effect and all random effects
  ZU1 <- list(NULL, TRUE, 1:zlen) # For intercep (z if no intercep)
  ZU2 <- list(NULL, TRUE, 1:zlen)
  GRPNames <- list(NULL, TRUE, 1:zlen)


  amat <- matrix(0, zlen, 2)
  rgmat <- matrix(0, zlen, 2)
  ZUall <- NULL
  ZonlyDat <- data[,Zonly]

  # Handle case with and without intercept

  getDummy = function(grpVec){
    grpNames <- names(table(grpVec))
    gN <- length(grpNames)
    n <- length(grpVec)
    ones <- rep(1, n)
    grpMat <- matrix(0, n, gN)
    for (g in 1: gN) {
      grpMat[grpVec == grpNames[g], g] <- 1
    }
    colnames(grpMat) <- grpNames
    return(grpMat)
  }

  for (r in 1:zlen){
    Zu <- NULL # for matix of all Xs in intercept and X order; restart counting


    # (i) there is no Z intercept (1),
    # (ii) only 1 is randomize; i.e only grpup means with out a covariate

    gVec <- data[,Gvars[r]]
    gMat <- getDummy(gVec)
    for (i in 1:2) { # i= 1 or Intecept if any;
      if (IntMat[r, i] == 1){
        if (i == 1) {
          ZUall <- cbind(ZUall, gMat)
          amat[r, 1] <- qr(gMat)$rank		# qr(Zu1)$rank
          rgmat[r, 1]<- ncol(gMat)			# ncol(Zu1)
        }
        if (i == 2){
          if (is.vector(ZonlyDat)) Zvals <- ZonlyDat
          else if (ncol(ZonlyDat) > 1) 	Zvals <- ZonlyDat[,r]
          else Zvals <- ZonlyDat
          Zu <- Zvals * gMat
          amat[r, 2] <- qr(Zu)$rank
          rgmat[r, 2] <- ncol(Zu)

          cnam <- colnames(gMat)
          colnames(Zu) <- paste(cnam, zonly[r], sep="-")
          ZUall <- cbind(ZUall, Zu)
        }
      }
    } # end of i
  }	# end of r


  zzlen <- ncol(ZUall)

  # compute rhos and valps
  rhos <- matrix(0, zlen, 2)
  valps <- matrix(0, zlen, 2)

  Zk = 0
  for (r in 1:zlen) Zk <- Zk + amat[r, 1] + amat[r, 2]
  SIGrho <- matrix(0, Zk, Zk)

  for (r in 1:zlen){
    set.seed(r)
    for (i in 1:2){

      if (i == 1 & IntMat[r, i] == 1){
        xxinvu <- xxinvu1[[r]]
        utilde <- u.tilde1[[r]]
        lambda <- lambda1[[r]]
      }
      if (i == 2 & IntMat[r, i] == 1){
        xxinvu <- xxinvu2[[r]]
        utilde <- u.tilde2[[r]]
        lambda <- lambda2[[r]]
      }

      if (IntMat[r, i] == 1){
        a <- amat[r, i] # qr(ztmp)$rank # This aj in paper
        U <- rchisq(M, e)
        U <- e * U / mean(U) # mean correct
        V <- rchisq(M, a)
        V <- a * V / mean(V)


        if (i==1){
          utilde <- u.tilde1[[r]]
          SSA <- c(t(utilde) %*% solve(xxinvu) %*% utilde) # SSA needed for each Varinace componeneny when rho=0
        }
        if (i==2){
          utilde <- u.tilde2[[r]]
          SSA <- c(t(utilde) %*% solve(xxinvu) %*% utilde) # SSA needed for each Varinace componeneny when rho=0
        }

        utilde <- utilde
        xxinvu <- xxinvu

        UVpart <- SSA / V - SSE / U
        UVpart <- UVpart[UVpart >= 0]
        U <- U[UVpart >= 0]
        V <- V[UVpart >= 0]
        U <- U[UVpart <= SSA/a]
        V <- V[UVpart <= SSA/a]

        ssvalp=function(valp, U, V, u.tilde, se, lambda){
          rhoU <- valp * U / se
          L <- length(lambda)
          srhoU <- 0
          for (l in 1:L)
            srhoU <- srhoU + u.tilde[l] ^ 2 / (rhoU + lambda[l])
          del <- V * se / srhoU - U
          return(mean(del))
        }

        mindel = function(valp, U, V, u.tilde , se, lambda){
          # function to calculate ssrho-cd
          del = (ssvalp(valp, U, V, u.tilde, se, lambda))^2
          return(del)
        }

        getalp = function(U, V, u.tilde , se, lambda){
          # get the variance component
          galp = optimize(mindel, lower = 0, upper = 1000, U = U, V = V, u.tilde = u.tilde, se = se, lambda = lambda, tol = 0.00001)
          return(galp$min)
        }

        valp <- getalp(U, V, utilde , se = SSE, lambda )
        rout[r, i] <- sqrt(valp)
        valps[r, i] <- valp

        # Estimate variance ratio
        Wae <- rf(M, a, e) # generate f random numbers
        # Wae <- (e/(e-2)) * Wae / mean(Wae)
        good <- Wae >= ((SSE / e) / (SSA / a))
        Wae <- Wae[good]
        cd <- a * (SSE / e) * mean(Wae)

        # ssrho for regression for estimating var ratio
        ssrho = function(rho, u.tilde, lambda){
          return(sum(u.tilde ^ 2/(rho + lambda)))
        }

        # Estimate variance comp or rationb by minimizing the expected difference
        mindiff1 = function(rho, cd, u.tilde, lambda){
          # function to calculate ssrho-cd
          diff = (ssrho(rho, u.tilde, lambda) - cd) ^ 2
          return(diff)
        }

        getrho = function(cd, u.tilde, lambda){
          # Get Variance Ratio
          grho = optimize(mindiff1, lower = 0, upper = .999, cd = cd, u.tilde = u.tilde, lambda = lambda ,tol = 0.00001)
          return (grho$min)
        }

        # Call grho function and estimate rho
        rhohat<-getrho(cd, utilde, lambda) # This is variance ratio
        rhos[r, i] <- rhohat
      }
    } # End of i loop
  } # End of r loop

  rhos<-rhos
  # Need to provide BLUPS for all factor levels and so do not use active cols from here on
  for (r in 1:zlen){
    for (i in 1:2){
      rhohat <- rhos[r, i]
      if (r==1 & i==1) rstart <- 1
      rgi <- rgmat[r, i]
      if (IntMat[r, i]==1){
        rend <- rstart + rgi -1
        rgi <- rend - rstart + 1
        SIGrho[rstart:rend, rstart:rend] <- diag(rhohat, rgi)
        rstart <- rend + 1
      } # End of if (IntMat[r,i]==1)
    } # End of i loop
  } # End of r loop
  rout[zlen1, 2] <- sqrt(MSE)

  if (1==0) { # For testing
    lambda <- lambda
    rhos <-rhos
    alps <- valps
    sigrho <- SIGrho
  }
  Zu <- ZUall # Gneralized method can make inferences on all random effects
  SIG <- MSE * diag(1, N)
  G <- diag(1, N)+ Zu %*% SIGrho %*% t(Zu)
  g <- G
  GI = solve(G)
  cG <- chol(G)
  cGI <- solve(t(cG))
  # Save orginal y and X since we need it to compute del belwo
  ySave <- y
  XSave <- X
  gy <- cGI %*% y
  gX <- cGI %*% X
  dat <- cbind(gy, gX)

  colnames(dat)[1:2] = c("y", "(Intercept)")
  if (1==1){ # Sanity Check: Direct Computation gives same estimate as glse below
    XtX <- t(X) %*% GI %*% X
    BetaHat <- solve((t(X) %*% GI %*% X)) %*% (t(X) %*% GI %*% y)
    # This is correct and exactly te same as lme reuslts
  }

  # Transformed X and Y for glse
  y <- dat[, 1]
  X <- dat[,-1]
  glse <- lm(y ~ - 1 + X)
  names(glse$coefficients) <- gsub('^.|$.', '', names(glse$coefficients))
  glse <- glse

  # Now Cpmpute BLUPs of random effects
  BetaHat <- glse$coef
  del <- ySave - XSave %*% BetaHat

  zzSigI <- t(Zu) %*% Zu + solve(SIGrho)
  uHat <- solve(zzSigI) %*% t(Zu) %*% del
  if (1 == 1){ #For testing only
    betahat <- BetaHat
    uhat <- uHat
  }
  otList <- list("StdDev of Random Effects, Predictors and Estimates of Fixed Effects")
  if(nz > 1){
    BLUPs <- list( "Predictors for "= paste(as.character(Znames[1]), "and", as.character(Znames[2])))
  } else{
    BLUPs <- list( "Predictors for " = as.character(Znames))
  }

  RAND1 <- list(0, TRUE, 1:zlen)
  RAND2 <- list(0, TRUE, 1:zlen)
  RAND <- list(NULL, TRUE, 1:2)

  rstart <- 1
  for (r in 1:zlen){
    cname <- as.character(Zonly[r])
    rg <- rgmat[r, 1]
    if (rg == 0) rg <- rgmat[r, 2]
    Rmat <- (matrix(0, rg, 2))
    colnames(Rmat) <- c("(Intercept)", cname)
    datvec <- data[,Gvars[r]]
    grpNames <- names(table(datvec))
    rownames(Rmat) <-	grpNames # t(GRPNames[[r]])

    for (i in 1:2){
      rgi <- rgmat[r, i]
      if (IntMat[r, i] == 1){
        rend <- rstart + rgi -1
        rgi <- rend - rstart + 1
        uhat = uHat[rstart:rend]
        Rmat[,i] <- uhat
        rstart <- rend + 1
      } # End of if (IntMat[r,i]==1)
    } #end of i
    RAND[[r]] <- Rmat
  } # end of r

  otList$Fixed <- list(glse$coef, summary(glse))
  otList$Pred <- list("Predictors of Random Effects", RAND)

  result <- list()
  coefficients <- list()

  result$fixed$coef <- glse$coef
  result$fixed$summary <- summary(glse)

  result$sd <- list("StdDev of Random Effects", rout)

  result$coefficients$fixed <- otList$Fixed[[1]]
  result$coefficients$random <- otList$Pred[[2]][[1]]

  if(method == "GM"){
    cat("Generalized linear mixed-effects model fit \n")
    cat("Data:", data_name, "\n")
    cat("Fixed:", deparse(fixed), "\n")
    print(result$fixed[[1]])
    cat("\nRandom effects: \n")
    cat("Formula:", deparse(random), "\n")
    stddev <- select(subset.data.frame(melt(result$sd[[2]]), X2 == "SD of Random Effect"), X1, value)
    randdev <- t(stddev$value)
    colnames(randdev) <- stddev$X1
    rownames(randdev) <- "StdDev"
    print(randdev)
    cat("\n Number of Observations: ", dim(data)[1], "\n")
  }
}
