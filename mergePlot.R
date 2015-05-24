pdf("merge4Loci3(2).pdf")
library(plotrix, lib="~/R")
load("phenAll.RData")
#tnmrg=read.delim("Tb.NPl.merge/Tb.NPl.chr8:4500000-44400000.merge.txt", sep=" ")
#thmrg=read.delim("Tb.Th3D.merge/Tb.Th3D.chr7:62200000-75300000.merge.txt", sep=" ")
#bmrg=read.delim("BVTV.merge/BVTV.chr2:1.22e+08-130700000.merge.txt", sep=" ") 
#cmrg=read.delim("Conn.D.merge/Conn.D.chr14:69300000-91500000.merge.txt", sep=" ")
#spn=strsplit(as.character(tnmrg$sdp), ".", fixed=T)
#spn=sapply(spn, function(x) as.numeric(x));spnidx=which(sapply(1:ncol(spn), function(x) sum(spn[,x] %in% 3:10)>0));tnmrg$multiallelic=0; tnmrg$multiallelic[spnidx]=1
#sph=strsplit(as.character(thmrg$sdp), ".", fixed=T)
#sph=sapply(sph, function(x) as.numeric(x));sphidx=which(sapply(1:ncol(sph), function(x) sum(sph[,x] %in% 3:10)>0)); thmrg$multiallelic=0; thmrg$multiallelic[sphidx]=1
#spc=strsplit(as.character(cmrg$sdp), ".", fixed=T)
#spc=sapply(spc, function(x) as.numeric(x));spcidx=which(sapply(1:ncol(spc), function(x) sum(spc[,x] %in% 3:10)>0)); cmrg$multiallelic=0; cmrg$multiallelic[spcidx]=1
#spb=strsplit(as.character(bmrg$sdp), ".", fixed=T)
#spb=sapply(spb, function(x) as.numeric(x));spbidx=which(sapply(1:ncol(spb), function(x) sum(spb[,x] %in% 3:10)>0)); bmrg$multiallelic=0; bmrg$multiallelic[spbidx]=1



wn=5123
w=1452
wc=8602
wh=4535
thresh=4.5
lwd=0.3
pch=19
intb=which(bmrg$bp>bv$from.bp[w-5] & bmrg$bp<bv$to.bp[w+1])
indb=which(bmrg$logP.merge[intb]>thresh)
intc=which(cmrg$bp>con$from.bp[wc-9] & cmrg$bp<con$to.bp[wc+15])
indc=which(cmrg$logP.merge[intc]>thresh)
inth=which(thmrg$bp>th$from.bp[wh-3] & thmrg$bp<th$to.bp[wh+6])
indh=which(thmrg$logP.merge[inth]>thresh)
intn=which(tnmrg$bp>tn$from.bp[wn-10] & tnmrg$bp<tn$to.bp[wn+6])
indn=which(tnmrg$logP.merge[intn]>thresh)
intb2=which(bmrg$bp>bv$from.bp[w-44-5] & bmrg$bp<bv$to.bp[w-44+5])
indb2=which(bmrg$logP.merge[intb2]>thresh)

 col=1#"#660033"
 colMrg="#CCCCCC"
 colMul="#CC0000"
 par(mfrow=c(3,2))#)
 plot(main="Conn.D", cex.main=0.93,cmrg$bp[intc][indc]/1e6, cmrg$logP.merge[intc][indc], ylim=c(0,8.5),lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(72.42,max(cmrg$bp[intc][indc]/1e6)))
 ind=which(cmrg$logP.merge[spcidx]>thresh & cmrg$bp[spcidx]>73000000)
# ind2=which(cmrg$bp>72420000)
# ind3=which(ind %in% ind2)
 points(cmrg$bp[spcidx][ind]/1e6, cmrg$logP.merge[spcidx][ind], pch=18, col=colMul)
 lines(cmrg$bp[intc][indc]/1e6, cmrg$logP.interval[intc][indc], ylim=c(0,8.5), col=col, lwd=2)
	rect(75237306/1e6,0.5,75288508/1e6,1.1, col=3)
	rect(73194749/1e6,0.5,73330793/1e6,1.1, col=4)
	rect(74754871/1e6,0.5,74952722/1e6,1.1, col=5)
	rect(74727348/1e6,1.1,74754650/1e6,1.7, col=6)


	legend("topleft", legend=c("MA", "BA","Cpb2", "Rb1", "Lrch1", "Esd"), col=c(colMul, colMrg, 3,4,5,6), pch=c(18,1,15,15,15,15), cex=0.71)
	title(main="A", adj=0)
 plot(main="BV/TV", cex.main=0.93, bmrg$bp[intb][indb]/1e6, bmrg$logP.merge[intb][indb], ylim=c(0,8.5), lwd=lwd, col=colMrg, xlim=c(129.63,130.8), pch=pch, xlab="Mb", ylab="logP")

 ind=which(bmrg$logP.merge[spbidx]>thresh & bmrg$bp[spbidx]>129800000)
 points(bmrg$bp[spbidx][ind]/1e6, bmrg$logP.merge[spbidx][ind], pch=18, col=colMul)

 lines(bmrg$bp[intb][indb]/1e6, bmrg$logP.interval[intb][indb], ylim=c(0,8.5), col=col, lwd=2)
        rect(130.576459,0.5,130.587523,1.1, col=3)
        rect(130.571252,1.1,130.582053,1.7, col=1)
	rect(130.27,0.5,130.28,1.1,col=2)
	legend("topleft", legend=c("MA", "BA","Oxt","Avp","Nop56"), col=c(colMul,colMrg,3,1,2), pch=c(18,1,15,15,15),cex=0.71)
	title(main="B", adj=0)
 plot(main="Tb.N", cex.main=0.93,tnmrg$bp[intn][indn]/1e6, tnmrg$logP.merge[intn][indn], ylim=c(0,8.5), lwd=lwd, col=colMrg, xlim=c(36.1,max(tnmrg$bp[intn][indn]/1e6)), pch=pch, xlab="Mb",ylab="logP")
 ind=which(tnmrg$logP.merge[spnidx]>thresh & tnmrg$bp[spnidx]>36.8*1e6)
 points(tnmrg$bp[spnidx][ind]/1e6, tnmrg$logP.merge[spnidx][ind], pch=18, col=colMul)
 lines(tnmrg$bp[intn][indn]/1e6, tnmrg$logP.interval[intn][indn], ylim=c(0,8.5), col=col, lwd=2)
 ablineclip(h=4.9, x1=min(tnmrg$bp[intn][indn]/1e6),x2=max(tnmrg$bp[intn][indn]/1e6), lty=2, col="gray", lwd=1.7)
	rect(37.517382,0.5,38.666439,1.1,col="orange")
	rect(40.274167,0.5,40.313316,1.1,col="coral")
	rect(40.487547,0.5,40.520789,1.1,col="8")
	legend("topleft", legend=c("MA", "BA", "Sgcz","Fgf20","Cnot7"), col=c(colMul, colMrg, "orange","coral","8"), pch=c(18,1,15,15,15), cex=0.71)
	title(main="C", adj=0)
 plot(main="Tb.Th", cex.main=0.93, thmrg$bp[inth][indh]/1e6, thmrg$logP.merge[inth][indh], ylim=c(0,9.5), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(63.37, max(thmrg$bp[inth][indh]/1e6)))
 ind=which(thmrg$logP.merge[sphidx]>thresh & thmrg$bp[sphidx]>64.1*1e6)
 points(thmrg$bp[sphidx][ind]/1e6, thmrg$logP.merge[sphidx][ind], pch=18, col=colMul)

 lines(thmrg$bp[inth][indh]/1e6, thmrg$logP.interval[inth][indh], ylim=c(0,9.5), col=col, lwd=2)
	 rect(65856837/1e6,0.5,66055280/1e6,1.1,col="blue")
	 rect(64501706/1e6,0.5,64753870/1e6,1.1,col="brown")#Apba2
	 rect(64346758/1e6,1.1,64374095/1e6,1.7,col=101)#fan1
#	 rect(64392645/1e6,1.1,64412121/1e6,1.7,col=102)#mcee
	 rect(64376576/1e6,0.5,64392268/1e6,1.1,col=106) #mphosph10

	 rect(64867052/1e6,0.5,64872997/1e6,1.1,col=107)#Ndnl2

	# rect(65296165/1e6,0.5,65371239/1e6,1.1,col=570)#Tjp1


	
	 legend("topleft", legend=c("MA", "BA", "Pcsk6", "Apba2", "Fan1","Mphosph10","Ndnl2"), col=c(colMul, colMrg, "blue", "brown",101,106,107), pch=c(18,1,15,15,15,15,15,15,15), cex=0.71)
	 title(main="D", adj=0)
#	 abline(v=c(64.5,64.75,64.86,64.87, 64.2, 64.3, 64.34, 64.37))
#	 abline(v=c(64.39, 64.41), col="red")
 plot(main="BV/TV", cex.main=0.93,bmrg$bp[intb2][indb2]/1e6, bmrg$logP.merge[intb2][indb2], ylim=c(0,8.5), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(121.9,max(bmrg$bp[intb2][indb2]/1e6)))
 ind=which(bmrg$logP.merge[spbidx]>thresh & bmrg$bp[spbidx]>121.9*1e6)
 points(bmrg$bp[spbidx][ind]/1e6, bmrg$logP.merge[spbidx][ind], pch=18, col=colMul)

 lines(bmrg$bp[intb2][indb2]/1e6, bmrg$logP.interval[intb2][indb2], ylim=c(0,8.5), col=col, lwd=2)
	 rect(122142790/1e6,0.5,122155639/1e6,1.1,col="purple")
	legend("topleft", legend=c("MA", "BA", "B2m"), col=c(colMul,colMrg,"purple"), pch=c(18,1,15), cex=0.71)
	title(main="E", adj=0)

 dev.off()
