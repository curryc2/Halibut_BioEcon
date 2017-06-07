# R-scripts for the FFP package.
# LIST OF FUNCTIONS.
# 	getAgeSchedules
# 	explogis
#	getSelectivities
# 	eqModel
#   yieldEquivalence

#NOTE: Funtions originally written by Steve Martell



#' Calculate the equivalent yield by sequentially removing each fishery.
#' @param ffp List of model parameters and age-schedule information.
#' @export
yieldEquivalence <- function(ffp)
{

		G 	<- length(ffp$MP$pYPR)
		x   <- rep(1,length=G)
		D 	<- rbind(rep(1,length=G),1-diag(1,G))

		fn 	<- function(x){
			ffp$MP$pscLimit <- NA
			if(ffp$HP$type=="YPR"){
				ak  	<- x * ffp$MP$pYPR
				ffp$MP$pYPR <- ak / sum(ak)
			}
			if(ffp$HP$type=="MPR"){
				ak  	<- x * ffp$MP$pMPR
				ffp$MP$pMPR <- ak / sum(ak)
			}

			ffp 	<- getFstar(ffp)
			cat("\n Fstar = ",ffp$HP$fstar,"\n")
			EM 		<- eqModel(ffp)
			EM$x    <- x
			return(EM)
		}

		XX 	<- apply(D,1,fn)
		df  <- ldply(XX,data.frame)
		Y   <- matrix(df$ye,ncol=G,byrow=TRUE)
		y  	<- Y[1,]
		M   <- Y[-1,]
		E  	<- t(t(M)-y)/y
		E[E==-1] <- NA
		out <- list(E=E,Y=Y)
		return(out)
}



#' Calculate the corresponding Fspr given SPR target
#' @param ffp List of model parameters and age-schedule information.
#' @export
getFstar <- function(ffp)
{
	if( any(!is.na(ffp$MP$pscLimit)) ){
		cat("\n PSC LIMIT\n")
		ffp <- getFsprPSC(ffp)
	}
	else{
		cat("\n NO PSC LIMIT\n")
		ffp <- getFspr(ffp)	
	}
	return(ffp)
}


#' Find the F* that is consistent with SPR target and PSC limits.
#' @param ffp List of model parameters and age-schedule information.
#' @details  This algorithm finds the value of F* that is consistent with
#' the sprTarget and one or more fixed PSC limits.  Must also determine
#' what type of allocation is being used (YPR or MPR) for setting allocations.
#' @export
getFsprPSC <- function(ffp)
{
	with(ffp,{
		ak    <- switch(HP$type,YPR=MP$pYPR,MPR=MP$pMPR,FPR=MP$pFPR)
		bGear <- !is.na(MP$pscLimit)
		iGear <- which(!is.na(MP$pscLimit))
		pk    <- ak[!bGear]/sum(ak[!bGear])

		getFootPrint <- function(phi)
		{
			tmp        <- ak
			tmp[bGear] <- phi[-1]
			tmp[!bGear]<- (1-sum(tmp[bGear]))*pk
			return(tmp)
		}

		# Objective function
		fn <- function(phi)
		{
			ffp$HP$fstar <- exp(phi[1])
			# ffp$MP$pYPR  <- getFootPrint(phi)
			tmp  <- getFootPrint(phi)
			if(HP$type=="YPR") ffp$MP$pYPR = tmp
			if(HP$type=="MPR") ffp$MP$pMPR = tmp
			if(HP$type=="FPR") ffp$MP$pFPR = tmp

			EM       <- run.ffp(ffp)
			spr  	 <- EM$spr
			psc      <- EM$ye
			
			ofn   	 <- (spr-ffp$HP$sprTarget)^2 + 
						10*sum(na.omit((psc-ffp$MP$pscLimit)^2))
			return(ofn)
		}
		params <- c(log(HP$fstar),ak[iGear])
		fit    <- nlminb(params,fn)
		print(fn(fit$par))

		ffp$HP$fstar <- exp(fit$par[1])
		if(ffp$HP$type=="YPR") ffp$MP$pYPR <- getFootPrint(fit$par)
		if(ffp$HP$type=="MPR") ffp$MP$pMPR <- getFootPrint(fit$par)	
		return(ffp)
	})
}

#' Find the F* that is consistent with SPR target given allocations.
#' @param ffp List of model parameters and age-schedule information.
#' @export
getFspr <- function(ffp)
{
	fn <- function(log.fspr)
	{
		ffp$HP$fstar <- exp(log.fspr)
		spr          <- run.ffp(ffp)$spr
		ofn          <- (spr-ffp$HP$sprTarget)^2
		return(ofn)
	}
	fit <- nlminb(log(ffp$HP$fstar),fn)
	ffp$HP$fstar <- exp(fit$par)
	return(ffp)

}

#' Run the Equilibrium Model
#' @param MP A list object containing the list of model parameters and selectivity parameters.
#' @export
run.ffp <- function(ffp)
{
	# Compute age-schedule information.
	if(with(ffp,!exists("ageSc"))) ffp <- getAgeSchedules(ffp)

	# Get length-based selectivity coefficients.
	if(with(ffp,!exists("selex"))) ffp <- getSelectivities(ffp)
	
	# Run the equilibrium model.
	EM    <- eqModel(ffp)
	return(EM)
}

#' Run profile over fstar for given allocations
#' @param sprTarget Target SPR reference point for FSPR calculations
#' @param fmax Maximum value for Fstar to profile over.
#' @export
run.prf <- function(ffp,sprTarget=0.4,fmax=0.45,model=NULL)
{
	ftry <- seq(0,fmax,length=100)
	ffp$HP$sprTarget = sprTarget
	fn <- function(ftry)
	{
		ffp$HP$fstar <- ftry
		em <- run.ffp(ffp)
		return(em)
	}
	
	runs <- lapply(ftry,fn)
	df   <- plyr::ldply(runs,data.frame)
	if(!is.null(model)) df$model <- model
	class(df) <- c("prf",class(df))
	return(df)
}

#' Plot yield profile
#' @export
#' @method plot prf
plot.prf <- function(df)
{
	p <- ggplot(df,aes(fstar,ye))
	p <- p + geom_area(aes(fill=gear),alpha=0.50)
	p <- p + labs(x="F*",y="Equilibrium yield")
	print(p)
}

# 
# EQUILIBRIUM MODEL
# eqModel <- function(theta,selex,type="YPR")
# 
#' Equilibrium age- sex-structured model for multiple fleets.
#' @param ffp List of model parameters and age-schedule information.
#' @param ...
#' @export
eqModel <- function(ffp)
{
	if( with(ffp,!exists("ageSc")) )
	{
		ffp <- getAgeSchedules(ffp)
		# with(ffp,print(theta$H))
	}

	if( with(ffp,!exists("selex")) )
	{
		ffp <- getSelectivities(ffp)
		# with(ffp,print(selex))
	}

	# with(c(ffp$theta,ffp$ageSc,ffp$selex,ffp$HP,ffp$MP),{
	with(ffp, {
		# ngear <- ffp$selex$ngear
		nage  <- theta$A
		ngear <- dim(MP)[1]
		nsex  <- theta$H
		va    <- as.array(selex)
		wa    <- ageSc$wa
		fa    <- ageSc$fa
		lz  <- matrix(1/nsex,nrow=nsex,ncol=nage)
		za  <- matrix(0,nrow=nsex,ncol=nage)
		qa  <- array(0,dim=c(nsex,nage,ngear))
		pa  <- array(0,dim=c(nsex,nage,ngear))
		ra  <- array(0,dim=c(nsex,nage,ngear))
		dlz <- array(0,dim=c(nsex,nage,ngear))

		# Survivorship under fished conditions at fstar
		fbar <- HP$fstar
		lambda <- rep(1.0,length=ngear)
		for(iter in 1:(ngear+1))
		{
			# total mortality and survival rates
			fe <- fbar * lambda
			# browser()
			for(h in 1:nsex)
			{
				# print(fe)
				if(dim(va)[3] > 1){
					fage   <- rowSums(fe*va[h,,])
				}
				else if(dim(va)[3] == 1){
					fage   <- fe * va[h,,]
				}
				za[h,] <- ageSc$mx[h,] + fage
			}
			sa <- exp(-za)
			oa <- 1.0 - sa

			# per recruit yield & numbers for gear k
			for(k in 1:ngear)
			{
				pa[,,k] <- va[,,k] * oa / za
				qa[,,k] <- va[,,k] * wa * oa / za
				ra[,,k] <- va[,,k] * fa * oa / za
			}

			#  survivorship
			for(j in 2:nage)
			{
				lz[,j] <- lz[,j-1] * sa[,j-1]
				if(j == nage)
				{
					lz[,j] <- lz[,j] / oa[,j]
				}

				# partial derivatives
				# for(k in 1:ngear)
				# {
				# 	dlz[,j,k] <- sa[,j-1]*(dlz[,j-1,k]-lz[,j-1]*va[,j-1,k])
				# 	if(j == A)
				# 	{
				# 		dlz[,j,k] <- dlz[,j,k]/oa[,j] - lz[,j-1]*sa[,j-1]*va[,j,k]*sa[,j]/oa[,j]^2
				# 	}
				# }
			}


			# Fmultipliers for fstar based on allocations
			qp    <- switch(HP$type,YPR=qa,MPR=pa,FPR=ra)
			ak    <- switch(HP$type,YPR=MP$pYPR,MPR=MP$pMPR,FPR=MP$pFPR)
			phi.t <- 0
			for(h in 1:nsex)
			{
				phi.t <- phi.t + as.vector(lz[h,] %*% qp[h,,])
				# phi.t <- phi.t + as.vector((lz[h,]*fa[h,]) %*% qp[h,,])
			}
			
			
			lam.t  <- ak / (phi.t/sum(phi.t))
			lambda <- lam.t / sum(lam.t)
			# cat(iter," lambda = ",lambda,"\n")
		}

		# incidence functions
		phi.e  <- sum(lz*fa)
		phi.q  <- phi.m <- phi.r <- dphi.e <- dre <- 0
		# dphi.q <- matrix(0,ngear,ngear) 
		for(h in 1:nsex)
		{
			# dphi.e <- dphi.e + as.vector(fa[h,] %*% dlz[h,,])
			phi.q  <- phi.q + as.vector(lz[h,] %*% qa[h,,])
			phi.m  <- phi.m + as.vector(lz[h,] %*% pa[h,,])
			phi.r  <- phi.r + as.vector(lz[h,] %*% ra[h,,])
			# derivatives for yield equivalence
			# for(k in 1:ngear)
			# {
			# 	for(kk in 1:ngear)
			# 	{
			# 		va2 <- va[h,,k] * va[h,,kk]
			# 		dqa <- va2*wa[h,]*sa[h,]/za[h,] - va2*oa[h,]*wa[h,]/za[h,]^2
			# 		dpq  <- as.double(qa[h,,k] %*% dlz[h,,k] + lz[h,] %*% dqa)
			# 		dphi.q[k,kk] <- dphi.q[k,kk] + dpq
			# 	}
			# }
		}
		spr   <- (phi.e/ageSc$phi.E)
		ispr  <- 1.0 / spr


		# equilibrium recruitment & sensitivity
		re    <- max(0,ageSc$ro*(theta$kappa-ispr)/(theta$kappa-1))
		be    <- re * phi.e
		# dre   <- ro * phi.E * dphi.e / (phi.e^2 * (kappa-1))
		
		# yield per recuit and yield
		ypr   <- fe * phi.q
		ye    <- re * ypr

		# mortality per recruit and mortality
		mpr   <- fe * phi.m
		me    <- re * mpr

		# spawning capital per recruit
		bpr   <- fe * phi.r
		
		# halibut per ton (metric for catch composition)
		hpt   <- floor(2204.62/(ypr/mpr))

		# fishery footprint
		footprint <- mpr/sum(mpr)*(1-HP$sprTarget)
		

		# Jacobian for yield
		# dye <- matrix(0,nrow=ngear,ncol=ngear)
		# for(k in 1:ngear)
		# {
		# 	for(kk in 1:ngear)
		# 	{
		# 		dye[k,kk] = fe[k]*re*dphi.q[k,kk] + fe[k]*phi.q[k]*dre[kk]
		# 		if(k == kk)
		# 		{
		# 			dye[k,kk] = dye[k,kk] + re*phi.q[k]
		# 		}
		# 	}
		# }
		# dye   <- re * as.vector(phi.q * diag(1,ngear)) + fe * phi.q * dre + fe * re * dphi.q
		# dye   <- re * (phi.q) + fe * phi.q * dre + fe * re * dphi.q
		
		# Yield equivalence
		# v  <- sqrt(diag(dye))
		# M  <- dye / (v %o% v)
		# print(v %o% v)
		# print(t(matrix(v,4,4)))

		# print(M)
		# cat("ye\n",ye,"\n")


		# Equilibrium Model output
		out <- list(
		            "gear" = MP$sector,
		            "fstar" = fbar,
		            "fe"  = fe,
		            "ye"  = ye,
		            "me"  = me,
		            "be"  = be,
		            "re"  = re,
		            "spr" = spr,
		            "ypr" = ypr,
		            "mpr" = mpr,
		            "bpr" = bpr,
		            "hpt"  = hpt,
		            "footprint" = footprint
		            )

		return(out)
	})
}


# 
# EXPONENTIAL LOGISTIC FUNCTION FOR SELECTIVITY
# 
#' Expontnential logistic function
#' @details From Thompson, 1994 paper for a simple 3 parameter dome-shape selectivity function.
#' @param   x vector of length-intervals
#' @param   mu length-at-50\%  selectivity
#' @param   sd standard deviation
#' @param   g shape parameter (0-1)
#' @export
explogis <- function(x,mu,sd,g)
{
	s <- (1/(1-g))*((1-g)/g)^g
	p <- plogis(x,mu,sd)*exp(g*(mu-x)/sd)
	return(s*p)
}


# 
# GET SELECTIVITIES 
# 
#' Calculate age-specific selectivity coefficients for fishery given a length-at-age and size-limits
#' @param ffp a list object representing selectivity parameters, size-limits, and discard mortality rtes.
#' @details  Based on the joint probability of capture and retention.
#' @export
getSelectivities <- function(ffp)
{
	if( with(ffp,!exists("ageSc")) )
	{
		ffp <- getAgeSchedules(ffp)
	}

	with (c(ffp),{
		nage  <- theta$A
		nsex  <- theta$H
		ngear <- dim(MP)[1]
		la    <- ageSc$la
		va <- array(0,dim=c(nsex,nage,ngear))
		all.ra <- array(0,dim=c(nsex,nage,ngear))
		all.sc <- array(0,dim=c(nsex,nage,ngear))
		for(k in 1:ngear)
		for(h in 1:nsex)
		{

			sc <- explogis(la[h,],MP$slx1[k],MP$slx2[k],MP$slx3[k])
			
			# sc <- explogis(la[h,],ffp$MP$slx1[k],ffp$MP$slx2[k],ffp$MP$slx3[k])
			
			sc[MP$slx4[k]:nage] <- sc[MP$slx4[k]-1]
			
			# sc[ffp$MP$slx4[k]:nage] <- sc[ffp$MP$slx4[k]-1]
      
			#CURRY:
			all.sc[h,,k] <- sc
			
			ra <- plogis(la[h,],MP$slim[k],0.1*la[h,])- plogis(la[h,],MP$ulim[k],0.1*la[h,])
			#CURRY:
			all.ra[h,,k] <- ra
			
			da <- (1-ra)*MP$dmr[k]
			va[h,,k] <- sc*(ra+da)
		}
    
# 		print('sc')
# 		print(sc)
#     print('ra')
#     print(ra)
#     print('da')
#     print(da)
		
		ffp$selex <- va
		ffp$probCap <- all.sc
		ffp$probRetain <- all.ra
		return(ffp)
	})
}


# 
# AGE SCHEDULE INFORMATION
# 
#' Calculate age-schedule information given life-history parameters.
#' @param ffp A list object containing model dimensions (age-sex), 
#'        stock parameters (Natural mortlaity, unfished biomass, steepness)
#'        and life-history parameters for growth \& maturity.
#' @export
getAgeSchedules <- function(ffp)
{
	with(ffp$theta,{
		# Length-at-age, weight-at-age, mature weight-at-age
		# vonb <- function(age,linf,k,to) {return( linf*(1-exp(-k*(age+to))) )}
		vonb <- function(age,linf,k,to) {return( linf*(1-exp(-k*(age-to))) )}
		la <- sapply(age,vonb,linf=linf,k=vbk,to=to)
		wa <- a*la^b #Error here
		ma <- sapply(age, plogis, location=ahat, scale=ghat);  ma[,1:7] <- 0
		fa <- ma * wa

		# Age-dependent natural mortality (Lorenzen)
		getM <- function(age,vbk,m,cm){
			t1 <- exp(vbk*(age+2))-1
			t2 <- exp(vbk*(age+1))-1
			sa <- (t1/t2)^(-m/vbk)
			mx=m*((log(t1)-log(t2))/vbk)^cm
			return(mx)
		}
		mx <- sapply(age,getM,m=m,vbk=vbk,cm=cm)

		# Survivorship at unfished conditions
		lx <- matrix(0,H,A)
		for(i in age)
		{
			if(i == min(age))
			{
				lx[,i] <- 1.0/H
			}
			else
			{
				lx[,i] <- lx[,i-1] * exp(-mx[,i-1])
			}
			if(i == A)
			{
				lx[,i] <- lx[,i] / (1-exp(-mx[,i]))
			}
		}
		phi.E <- sum(lx*fa)
		ro    <- bo/phi.E

		# List objects to return
		ffp$ageSc <- list(ro=ro,la=la,wa=wa,ma=ma,fa=fa,lx=lx,mx=mx,phi.E=phi.E)
		# theta <- c(theta,ageSc)
		return(ffp)
	})
}



