options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
        stop("Wrong parameters. Run as: Rscript generateDB.R <directory with WGKS> <output file>", call.=FALSE)
}

dirname = args[1]
outfile = args[2]

message("Reading WGKS list...")
filelist = dir(dirname)

message("Creating database...")
db = data.frame()
for (f in filelist) {
	ff = paste(dirname,"/",f,sep="")
	ffdata <- as.matrix(read.table(ff,header=F,skip=3,row.names=1,sep="\t",check.names=FALSE))[,3]
	print(paste(ff,"read in..."))
	rown = rownames(ffdata)
	ffdata2 = sprintf("%.3f", ffdata)
	rownames(ffdata2) = rown
	db = merge(db,ffdata,by=0,all=T)
	row.names(db)=db$Row.names
	db$Row.names=NULL
	colnames(db)[dim(db)[2]]=gsub(".fa.*","",f)
}

message("Writing database...")
write.table(db, file=outfile, col.names=T, quote=F, sep="\t")
