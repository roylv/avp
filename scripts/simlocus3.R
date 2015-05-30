
select.locus <- function( db, chrstr, interval, dist.bp) {

	     chr.bp.vec <- db$additive$genome$bp[db$additive$genome$chromosome == chrstr]
             peak.bp <- chr.bp.vec[interval]

             ind.left <- (interval - 1)
             diff.left <- (peak.bp - chr.bp.vec[ind.left])
             while ( diff.left < dist.bp ) {
             	   ind.left = ( ind.left - 1 )
                   diff.left = (peak.bp - chr.bp.vec[ind.left])
             }
#             print(ind.left)

             ind.right <- (interval + 1)
             diff.right <- (chr.bp.vec[ind.right] - peak.bp)
             while ( diff.right < dist.bp ) {
             	   ind.right = ( ind.right + 1 )
              	   diff.right = (chr.bp.vec[ind.right] - peak.bp)
             }
#             print(ind.right)

	     interval.vector <- ind.left:ind.right

	     return(interval.vector)
}


define.region <- function( db, chrstr, interval, dist.bp ) {

	      region <- list()

	      chr.bp.vec <- db$additive$genome$bp[db$additive$genome$chromosome == chrstr]
              peak.bp <- chr.bp.vec[interval]

              ind.left <- (interval - 1)
              diff.left <- (peak.bp - chr.bp.vec[ind.left])
              while ( (diff.left < dist.bp) && (ind.left > 1) ) {
              	        ind.left = ( ind.left - 1 )
              		    diff.left = (peak.bp - chr.bp.vec[ind.left])
              }
              print(ind.left)

              ind.right <- (interval + 1)
              diff.right <- (chr.bp.vec[ind.right] - peak.bp)
              while ( diff.right < dist.bp ) {
              	        ind.right = ( ind.right + 1 )
              		    diff.right = (chr.bp.vec[ind.right] - peak.bp)
              }
              print(ind.right)

              Nloc <- ( ind.right - ind.left + 1 )
              print(Nloc)

	      region$ind.left <- ind.left
	      region$ind.right <- ind.right
	      region$Nloc <- Nloc
	      region$chr.bp <- chr.bp.vec

	      return(region)
}



extract.coeffs <- function( d, phen.data ) {

	       peak.coeff <- list()

	       genetic.model = as.formula( "phen.data$x ~ d")
	       fit.genetic = lm( as.formula(genetic.model), data=phen.data)
	       fit.genetic$coeff[9]=0
	       peak.coeff$beta <- fit.genetic$coeff[2:9]
               peak.coeff$beta[8] = 0

	       peak.coeff$const <- fit.genetic$coeff[1]
	       peak.coeff$log.hap <- ( d %*% fit.genetic$coeff[2:9] )
	       peak.coeff$pheno <- phen.data$x
	       #peak.coeff$log.pheno <- fit.genetic$y[,1]
	       peak.coeff$log.pheno <- fit.genetic$effects

	       return( peak.coeff )	 
}




sim.locus <- function( chr.list, peak.coeff, d.locus, db.use, region.stats, sim.int, loc ) {

	  locus.acc.factors <- exp( d.locus %*% peak.coeff$beta )
	  resid.vector <- exp( peak.coeff$log.pheno - peak.coeff$log.hap )

	  max.logP <- numeric(100)
	  max.interval <- numeric(100)
	  max.bp <- numeric(100)
	  sim.interval <- numeric(100)
	  sim.bp <- numeric(100)

	  sim.int.bp <- region.stats$chr.bp[ sim.int ]

	  for ( i in 1:1000 ) {

	      cat("locus ", loc, "sim.int = ", sim.int, " sim ", i, "\n")

	      resid.perm <- sample( resid.vector, replace=FALSE )

	      T <- ( locus.acc.factors * resid.perm )
	      #T <- ( T + 4.999 )
	      #censored <- (T < 28)
	      #T[ T > 28 ] = 28

	      logP.vec <- numeric(region.stats$Nloc)

	      for ( ind in region.stats$ind.left:region.stats$ind.right ) {

	      	  d <- chr.list[ind][[1]]$mat
		  rownames(d) <- db$subjects
		  d = d[db.use,]

		  fit.genetic = lm( T ~ d )
		  fit.null = lm( T ~ 1 )
		  a = anova(fit.null, fit.genetic)

		  logP.vec[ (ind - region.stats$ind.left + 1) ] <- -log10(a[2,6])#7])
	       }

	       max.logP[i] <- max(logP.vec)
	       max.int <- which( logP.vec == max.logP[i] )
	       max.interval[i] <- ( max.int + region.stats$ind.left - 1 )
	       max.bp[i] <- region.stats$chr.bp[ max.interval[i] ]
	       sim.interval[i] <- sim.int
	       sim.bp[i] <- sim.int.bp

#	       cat("max logP = ", max.logP[i], "max interval = ", max.interval[i], "\n")
	   }

	   diff.interval <- ( max.interval - sim.int )
	   diff.bp <- ( max.bp - sim.int.bp )

	    locus.df <- data.frame(max.logP=max.logP, max.interval=max.interval, max.bp=max.bp, diff.interval=diff.interval, diff.bp=diff.bp, sim.interval=sim.interval, sim.bp=sim.bp)

	    return(locus.df)
}



sim.qtl <- function( db, chr.list, interval, chrstr="chr8", dist.bp=5000000, data.file, pheno ) {

	library(survival)
        library(MASS)

#        data.file <- "asper.phenotypes.new.txt"
#        surv.time = "Surv.day"
#        surv.censored = 28

#        phen.data <- read.table(file=data.file, header=TRUE)
	phen.data <- read.delim(file=data.file, header=TRUE)
        err=0
        for ( p in c( "CC.Line", pheno )) {
	    if ( is.null( phen.data[[p]] ) ) {
               warning( "file ", data.file, " has no column ", p, "\n")
               err = err+1
            }
        }
        if ( err > 0)  return(NULL)

	agg.phen.data = aggregate( phen.data[[pheno]], by=list(CC.Line=phen.data$CC.Line), FUN="median" )
        #agg.phen.data$death = agg.phen.data$x < surv.censored

        int = intersect( db$subjects, unique(agg.phen.data$CC.Line))
        phen.data = agg.phen.data[agg.phen.data$CC.Line %in% int,]

        db.use = match( phen.data$CC.Line, db$subjects)


#	extract information from actual QTL

	d.peak <- chr.list[interval][[1]]$mat
        rownames(d.peak) <- db$subjects
	d.peak = d.peak[db.use,]

	peak.coeff <- extract.coeffs( d.peak, phen.data )


#	run simulation for actual peak

	sim.int = interval

	region.stats <- define.region( db, chrstr, sim.int, 2*dist.bp )

        d.int <- chr.list[sim.int][[1]]$mat
        rownames(d.int) <- db$subjects
        d.int = d.int[db.use,]

        all.df <- sim.locus( chr.list, peak.coeff, d.int, db.use, region.stats, sim.int, 1 )



	interval.vector <- select.locus( db, chrstr, interval, dist.bp )

	for ( loc in 2:100 ) {

#	    select interval to simulate peak

	    interval.vector.perm <- sample( interval.vector, replace=FALSE )
	    sim.int <- interval.vector.perm[1]


#	    for selected interval, run simulation

	    region.stats <- define.region( db, chrstr, sim.int, 2*dist.bp )

	    d.int <- chr.list[sim.int][[1]]$mat
	    rownames(d.int) <- db$subjects
	    d.int = d.int[db.use,]

	    locus.df <- sim.locus( chr.list, peak.coeff, d.int, db.use, region.stats, sim.int, loc )

	   all.df <- rbind(all.df, locus.df) 
	 }


#	    par( mfrow=c(1,2) )
#	    plot(locus.df$diff.bp/1.0e6, locus.df$max.logP, cex=0.75, ylim=c(0, (max(locus.df$max.logP) + 0.5) ), xlab="Mb" )
#	    truehist(locus.df$diff.bp/1.0e6, h=0.5, xlab="Mb")


	    write.table( all.df, file=paste("simCI", chrstr, interval,"BVTV21", "txt", sep="."), row.names=FALSE, quote=FALSE )

#	    return(all.df)
}
