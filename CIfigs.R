
calc.CIs <- function( data.df, low=0.025, high=0.975, obsN=100000 ) {

	 low.index <- obsN*low
	 high.index <- obsN*high
	 cat("low index = ", low.index, "high index = ", high.index, "\n")

	 diff.bp.srt <- sort(data.df$diff.bp)

	 cat("CI lower = ", diff.bp.srt[low.index]/1.0e6, "CI higher = ", diff.bp.srt[high.index]/1.0e6, "\n")

	 thres <- 3.827

	 


}



calc.int.CIs <- function( data.df, low=0.025, high=0.975, obsN=100000 ) {

         low.index <- obsN*low
         high.index <- obsN*high
         cat("low index = ", low.index, "high index = ", high.index, "\n")

         diff.int.srt <- sort(data.df$diff.interval)

         cat("CI lower = ", diff.int.srt[low.index], "CI higher = ", diff.int.srt[high.index], "\n")
}



all.plot <- function( data.df, pdffile=NULL ) {

	 library(MASS)

	 calc.pts <- seq(0, 13, 0.1)
	 obsN <- length(calc.pts)
	 cdf.diff.Mb <- numeric(obsN)

	 for ( i in 1:obsN ) {
	     cdf.diff.Mb[i] <- length( which( abs(data.df$diff.bp)/1.0e6 <= calc.pts[i] ) )/100000
	 }

	 Nint <- length(unique(data.df$sim.interval))

	 if ( ! is.null(pdffile) ) pdf(pdffile, width=10, height=15 )

	 par( mfrow=c(3,2) )

         truehist(data.df$diff.bp/1.0e6, h=0.5, xlab="Distance from true location (Mb)", col="blue" )

	 boxplot( data.df$max.logP ~ data.df$diff.interval, cex=0.5, ylim=c(0,(max(data.df$max.logP)+0.5)), xlab="Distance from true location (Number of intervals)", ylab="Peak logP")

	 plot( calc.pts, cdf.diff.Mb, type="l", ylim=c(0,1), xlab="Distance from the true location (Mb)", ylab="CDF" )

	 hist( data.df$sim.interval, xlab="Distribution of simulated locations", main="", col="blue", breaks=Nint )

	 boxplot( data.df$max.interval ~ data.df$sim.interval, cex=0.5, xlab="Simulated location", ylab="Location of peak" )
#         abline(h=true.qtl)

	 boxplot(data.df$diff.bp/1.0e6 ~ data.df$sim.interval, cex=0.5, xlab="Simulated location", ylab="Distance of peak from simulated location")

       	 if ( ! is.null(pdffile) ) dev.off()
}


err.dist.hist.plot <- function( data.df, pdffile=NULL ) {

		   if ( ! is.null(pdffile) ) pdf(pdffile, width=7.5, height=5)

		   truehist( data.df$diff.bp/1.0e6, h=0.5, xlab="Mb", ylab="Frequency", col="grey" )

		   if ( ! is.null(pdffile) ) dev.off()
}


err.distVSmax.logP.plot <- function( data.df, pdffile=NULL ) {

                   if ( ! is.null(pdffile) ) pdf(pdffile, width=7.5, height=7.5)

		   plot(abs(data.df$diff.bp/1.0e6), data.df$max.logP, cex=0.5, ylim=c(0,(max(data.df$max.logP)+0.5)), xlab="Distance from true location (Mb)", ylab="Peak logP")

                   if ( ! is.null(pdffile) ) dev.off()
}



simVSpeak.interval.box.plot <- function( data.df, true.qtl=155, pdffile=NULL ) {

                   if ( ! is.null(pdffile) ) pdf(pdffile, width=15, height=5)

		   boxplot( data.df$max.interval ~ data.df$sim.interval, cex=0.5, xlab="Simulated location", ylab="Location of peak" )
		   abline(h=true.qtl)

                   if ( ! is.null(pdffile) ) dev.off()
}


err.distVSsim.interval.box.plot <- function( data.df, pdffile=NULL ) {

                   if ( ! is.null(pdffile) ) pdf(pdffile, width=15, height=5)

		   boxplot(data.df$diff.bp/1.0e6 ~ data.df$sim.interval, cex=0.5, xlab="Simulated location", ylab="Distance of peak from simulated location")

                   if ( ! is.null(pdffile) ) dev.off()
}