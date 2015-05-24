#bv=read.delim("scan_1192.txt")
#con=read.delim("scan_ConnD14_408.txt")
#tnh=read.delim("scan_TbTh601.txt")

#wb=1452
#wc=8602
#wn=5123
#wh=4535
#wbl=548
#wnl=154
#whl=185
#wcl=286


#phen.data=bv
#phen="BVTV"
regHer=function(db, phen, phen.data, chr, w){
	agg.phen.data = aggregate( phen.data[[phen]],
                           by=list(CC.Line=phen.data$CC.Line),
                           FUN="median" )

	int = intersect( db$subjects, unique(agg.phen.data$CC.Line))
	phen.data = agg.phen.data[agg.phen.data$CC.Line %in% int,]
	db.use = match( phen.data$CC.Line, db$subjects)

	fit.null=lm(x~1,data=phen.data, y=TRUE)
	e=db$additive$chr[[chr]][[w]]$mat
	d=e[db.use,]	
	fit.genetic=lm(x~d, data=phen.data)
	a=anova(fit.null, fit.genetic)
	cat("Hr2=",(a[1,2]-a[2,2])/a[1,2], "Pv=",a[2,6], "\n", sep="\t")
	return(a)
}






