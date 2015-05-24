source("scripts/simlocus.R")
source("/Net/mus/data/www/preCC/BONE/LEVY/mf/7jan15/trab/march15/happy.preCC.R")
if(is.null(condensed.db)) condensed.db = load.condensed.database("/Net/mus/data/www/preCC/MEGA_MUGA/Mar2015.MEGA+MDA+MUGA/CONDENSED/")
db=condensed.db
#loc=matrix(c("chr7", 185,"scan_TbTh601.txt","chr2", 548, "scan_1192.txt","chr8",154, "scan_TbTh601.txt" ,"chr14",286,"scan_ConnD14_408.txt"),nrow=4,ncol=3, byrow=T)
runsim=function(chs){
  if(chs=="th"){
	 
		 w=185
		 n=100
		 phen="Tb.Th3D"
		 data.file="scan_TbTh601.txt"
		 chr.list=db$additive$chr$chr7
		 dfth=replicate(1000,sim.qtl(sim.int=w, n=n, data.file=data.file, phen=phen, db=db, chr.list=chr.list, interval=w, chrstr=7, dist.bp=5000000))
		 save(dfth, file="cisimTh.RData")
	}
  else if(chs=="bv"){
	         w=548
                 n=100
                 phen="BVTV"
                 data.file="scan_1192.txt"
                 chr.list=db$additive$chr$chr2
                 dfbv=replicate(1000,sim.qtl(sim.int=w, n=n, data.file=data.file, phen=phen, db=db, chr.list=chr.list, interval=w, chrstr=2, dist.bp=5000000))
                 save(dfbv, file="cisimbBv.RData")
	}
  else if(chs=="tbn"){
                 w=154
                 n=100
                 phen="Tb.NPl"
                 data.file="scan_TbTh601.txt"
                 chr.list=db$additive$chr$chr8
                 dftbn=replicate(1000,sim.qtl(sim.int=w, n=n, data.file=data.file, phen=phen, db=db, chr.list=chr.list, interval=w, chrstr=8, dist.bp=5000000))
                 save(dftbn, file="cisimbTbn.RData")
        }

  else if(chs=="con"){
                 w=286
                 n=100
                 phen="Conn.D"
                 data.file="scan_ConnD14_408.txt"
                 chr.list=db$additive$chr$chr14
                 dfcon=replicate(1000,sim.qtl(sim.int=w, n=n, data.file=data.file, phen=phen, db=db, chr.list=chr.list, interval=w, chrstr=14, dist.bp=5000000))
                 save(dfcon, file="cisimbCon.RData")
        }
}







