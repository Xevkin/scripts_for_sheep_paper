#to be run on the ancError output of the angsd -doAncError 1
#produces a tmp.txt file with all possible transitions removed (i.e. those involving any combination of the perfect genome, test genome, and outgroup)
#the output was then used in the input of estError.R script provided with angsd

args = commandArgs(trailingOnly=TRUE)

a<-read.table(args[1])

a[,c(3,8,11,12,13,14,15,18,23,29,34,39,41,42,43,44,45,49,51,52,53,54,55,56,61,66,71,77,81,82,83,84,85,87,92,97)]<-0

write.table(file="tmp.txt",a,row.names=F,col.names=F)
