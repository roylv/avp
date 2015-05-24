





define.region <- function( db, chrstr, interval, dist.bp ) {

	      	 region <- list()

		 chr.bp.vec <- db$additive$genome$bp[db$additive$genome$chromosome == chrstr]
              	 peak.bp <- chr.bp.vec[interval]

              	 ind.left <- (interval - 1)
              	 diff.left <- (peak.bp - chr.bp.vec[ind.left])
              
		
	
		while ( diff.left < dist.bp ) {
              	       ind.left = ( ind.left - 1 )
              	       diff.left = (peak.bp - chr.bp.vec[ind.left])
              	 }
              	 

              	 ind.right <- (interval + 1)
              	 diff.right <- (chr.bp.vec[ind.right] - peak.bp)
              	 while ( diff.right < dist.bp ) {
              	        ind.right = ( ind.right + 1 )
              		diff.right = (chr.bp.vec[ind.right] - peak.bp)
                 }
              	 

              	 Nloc <- ( ind.right - ind.left + 1 )
              	 

	         region$ind.left <- ind.left
		 region$ind.right <- ind.right
		 region$Nloc <- Nloc
		 region$chr.bp <- chr.bp.vec

		 return(region)
}



extract.coeffs <- function( d, phen.data ) {

	         peak.coeff <- list()

		 #genetic.model = as.formula( "Surv(phen.data$x,phen.data$death) ~ d")
		 genetic.model = as.formula( "phen.data$x ~ d")
		 fit.genetic = lm( as.formula(genetic.model), data=phen.data)
#genetic.model = paste( as.character(null.model), " + d")
print(summary(fit.genetic))
		 peak.coeff$beta <- fit.genetic$coeff[2:9]
		 peak.coeff$const <- fit.genetic$coeff[1]
		 #peak.coeff$resid <- ( fit.genetic$y[,1] - fit.genetic$linear.predictors )
		 peak.coeff$resid = resid(fit.genetic)

#	myLP <- ( d %*% fit.genetic$coeff[2:9] ) + fit.genetic$coeff[1]

	     	 #print(summary(phen.data$death))
#		 print("ex.coef")
		 return( peak.coeff )	 
}




sim.locus <- function( n, chr.list, peak.coeff, d.locus, db.use, region.stats, sim.int ) {

	  pb=peak.coeff$beta
	  pb[8]=0	
	  locus.lp <- ( d.locus %*% pb ) + peak.coeff$const
	  #n=1000	
	  max.logP <- numeric(n)
	  max.interval <- numeric(n)
	  max.bp <- numeric(n)

	  for ( i in 1:n ) {

	      cat("sim ", i, "\n")


	      resid.perm <- sample( peak.coeff$resid, replace=FALSE )

	      logT <- ( locus.lp + resid.perm )
	      T <- exp(logT)
	      #censored <- (T < 28)
	      #T[ T > 28 ] = 28

	      logP.vec <- numeric(region.stats$Nloc)

	      for ( ind in region.stats$ind.left:region.stats$ind.right ) {
		  
	      	  d <- chr.list[ind][[1]]$mat
		  rownames(d) <- db$subjects
		  d = d[db.use,]
		  #fit.genetic = survreg( Surv(T, censored) ~ d )
		 
		  fit.genetic = lm( T~1 + d )
				
		  #fit.null = survreg( Surv(T, censored) ~ 1 )
		  fit.null = lm(T ~ 1 )
		  a = anova(fit.null, fit.genetic)
		  logP.vec[ (ind - region.stats$ind.left + 1) ] <- -log10(a[2,6])
	       }

	       max.logP[i] <- max(logP.vec)
	       max.int <- which( logP.vec == max.logP[i] )
	       max.interval[i] <- ( max.int + region.stats$ind.left - 1 )
	       max.bp[i] <- region.stats$chr.bp[ max.interval[i] ]
	    }

	    sim.int.bp <- region.stats$chr.bp[ sim.int ]
	    diff.interval <- ( max.interval - sim.int )
	    diff.bp <- ( max.bp - sim.int.bp )

	    locus.df <- data.frame(max.logP=max.logP, max.interval=max.interval, max.bp=max.bp, diff.interval=diff.interval, diff.bp=diff.bp)

	    return(locus.df)
}



sim.qtl <- function(n, db, chr.list, interval, chrstr="chr8", dist.bp=5000000, data.file, phen, sim.int) {

	#library(survival)
        library(MASS)

        #data.file <- "scan_1192.txt"
        #phen = "BVTV"
        #surv.censored = 28

        phen.data <- read.delim(file=data.file, header=TRUE)
        
	err=0
        for ( p in c( "CC.Line", phen )) {
	    if ( is.null( phen.data[[p]] ) ) {
               warning( "file ", data.file, " has no column ", p, "\n")
               err = err+1
            }
        }
        if ( err > 0)  return(NULL)

	agg.phen.data = aggregate( phen.data[[phen]], by=list(CC.Line=phen.data$CC.Line), FUN="median" )
        #agg.phen.data$death = agg.phen.data$x < surv.censored

        int = intersect( db$subjects, unique(agg.phen.data$CC.Line))
        phen.data = agg.phen.data[agg.phen.data$CC.Line %in% int,]

        db.use = match( phen.data$CC.Line, db$subjects)


#	extract information from actual QTL

	d.peak <- chr.list[interval][[1]]$mat
        rownames(d.peak) <- db$subjects
	d.peak = d.peak[db.use,]
        
	
	peak.coeff <- extract.coeffs( d.peak, phen.data )

	print(summary(peak.coeff))



#	select interval to simulate peak

	#sim.int <- 548



#	for selected interval, run simulation

	region.stats <- define.region( db, chrstr, sim.int, 2*dist.bp )

	d.int <- chr.list[sim.int][[1]]$mat
	rownames(d.int) <- db$subjects
	d.int = d.int[db.use,]
 
	locus.df <- sim.locus(n=n, chr.list, peak.coeff, d.int, db.use, region.stats, sim.int )

	print(summary(locus.df))

	#par( mfrow=c(1,2) )
 	#pdf("cisim4.pdf")
	#plot(locus.df$diff.bp/1.0e6, locus.df$max.logP, cex=0.75, ylim=c(0, (max(locus.df$max.logP) + 0.5) ) )
	#truehist(locus.df$diff.bp/1.0e6, h=1)
	#dev.off()
	return(locus.df)
}
