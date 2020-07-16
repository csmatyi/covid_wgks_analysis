options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
        stop("Wrong parameters. Run as: Rscript generateDB.R <WGKS database> <input file>", call.=FALSE)
}

wgksdb = args[1]
infile = args[2]

message("Reading WGKS DB...")
db = as.matrix(read.table(wgksdb,header=T,row.names=1,sep="\t"))

message("Adding query to database...")
data = as.matrix(read.table(infile,header=F,skip=3,row.names=1,sep="\t",check.names=FALSE))[,3]
db_data = merge(db,data,by=0,all=T)
row.names(db_data)=db_data$Row.names
db_data$Row.names=NULL
colnames(db_data)[dim(db_data)[2]]=gsub(".fa.*","",infile)

message("Calculating correlation matrix...")
zcor = cor(db_data,use="complete.obs")
zcor2 = zcor[1:dim(zcor)[2]-1,1:dim(zcor)[2]-1]
zcor2v = zcor2[upper.tri(zcor2)]
xv = zcor[1:dim(zcor)[2]-1,dim(zcor)[2]]
p = t.test(zcor2v,xv)$p.value

print(paste("P-value is",p))
print(paste("Mean DB PCC is",mean(zcor2v)))
print(paste("Mean Query PCC is",mean(xv)))
