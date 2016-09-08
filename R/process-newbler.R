#--functions for analysing the Newbler 454ContigGraph.txt file--
#--by RW on 060916 using modified version of previous functions I have written

extract.contig.data<-function()
{
 #--assumes there is a file in ./ called "454ContigGraph.txt"
 fileCheck<-list.files()
 if(!("454ContigGraph.txt"%in%fileCheck))
 {
  stop("No 454ContigGraph.txt here!")
 }

 ind<-system("grep -n ^C 454ContigGraph.txt",intern=T)[1]
 res<-strsplit(ind,":",fixed=T)[[1]][1]
 res<-as.numeric(res)-1

 system(paste("head -n",as.character(res),"454ContigGraph.txt > tmp.txt"))
 contigSummData<-utils::read.table(file="tmp.txt",header=F,sep="\t",stringsAsFactors=F)
 system("rm tmp.txt")
 contigSummDataCut<-contigSummData[,3:4]
 colnames(contigSummDataCut)<-c("len","cov")
 rownames(contigSummDataCut)<-contigSummData[,2]
 return(contigSummDataCut)
}

add.gc<-function(contigSummData,fastaFile="454LargeContigs.fna")
{
 #--import contig sequences and calculate GC-content
 contigs<-Biostrings::readDNAStringSet(filepath=fastaFile)
 contigs.1merMat<-Biostrings::oligonucleotideFrequency(contigs,1)
 contigs.GC<-rowSums(contigs.1merMat[,c("G","C")])/rowSums(contigs.1merMat)
 newContigNames<-sapply(strsplit(names(contigs)," ",fixed=T),FUN=function(x){x[1]})
 names(contigs.GC)<-newContigNames
 #--find common set...
 common<-intersect(rownames(contigSummData),names(contigs.GC))
 #--make a new version of incorporates GC-data
 summData<-data.frame(contigSummData[common,],gc=contigs.GC[common])
 return(summData)
}


