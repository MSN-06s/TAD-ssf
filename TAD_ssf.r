options(warn=-1)
library(dplyr)
library(tidyverse)

#boundary.1<-args[1]
#boundary.2<-args[2]
#ref<-args[3]

print.usage <- function() {
  cat ("Perform TAD boundaries differential analysis.")
  cat("\n")
  cat("Usage:")
  cat("\n")
  cat("$ ./TAD_ssf.r boundary_file_1 boundary_file_2 reference_file gap_threshold number_of_chromosomes")
  cat("\n")
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  print.usage()
  stop()
}
boundary.1 <- args[1]
boundary.2 <- args[2]
#boundaries produced by Hicexploerer etc.
ref <- args[3]
#Matrix bins distribution on Chromosomes produced by hicpro
gap.threshold <- args[4]
#Threshold of num. of bins TAD shift to be identified as a shift
num.chromosomes <- args[5]
#Number of chromosomes except X/Y, for example, when analysis human data,num.chromosomes=22

s1<-read.table(sprintf("%s",boundary.1))
s2<-read.table(sprintf("%s",boundary.2))
re<-read.table(sprintf("%s",ref))
                                             
num.boundary.1<-nrow(s1)
num.boundary.2<-nrow(s2)
resolution<-re[1,3]-re[1,2]+0
                                             
cat (sprintf("The resolution of analysis is %s",resolution))
cat("\n")
cat (sprintf("The number of boundaries of sample1 is %s",num.boundary.1))
cat("\n")
cat (sprintf("The number of boundaries of sample2 is %s",num.boundary.2))
cat("\n")

s1$chr<-c('chr')
s2$chr<-c('chr')

s1$V1<-str_c(s1$chr,s1$V1)
s1<-s1[,-7]
s2$V1<-str_c(s2$chr,s2$V1)
s2<-s2[,-7]
re<-re[,-4]
s1<-s1[,-c(4,5,6)]
s2<-s2[,-c(4,5,6)]
colnames(s1)<-c('chr','start','end')
colnames(s2)<-c('chr','start','end')
colnames(re)<-c('chr','start','end')
re$start<-re$start+resolution/2
re$end<-re$end+resolution/2

s1$est<-str_c(s1$chr,'.',s1$start,'.',s1$end)
s2$est<-str_c(s2$chr,'.',s2$start,'.',s2$end)
re$est<-str_c(re$chr,'.',re$start,'.',re$end)

s1$bd<-c(1)
s2$bd<-c(1)

j1<-left_join(re,s1,by="est")
j1<-j1[,-c(4,5,6,7)]
j2<-left_join(re,s2,by="est")
j2<-j2[,-c(4,5,6,7)]

colnames(j1)<-c('chr','start','end','bd')
colnames(j2)<-c('chr','start','end','bd')

j1[is.na(j1)]<-0
j2[is.na(j2)]<-0

j<-cbind(j1,j2)
j<-j[,-c(5,6,7)]
colnames(j)<-c('chr','start','end','bd.s1','bd.s2')
j$bd.comp<-j$bd.s1-j$bd.s2
j$mul<-j$bd.s1*j$bd.s2
write.table(j,file="try1.txt",sep="\t")

#zero<-j[grep(pattern = "1",j[,7]),]
#zero.line<-row.names(zero)
#zero.line<-as.numeric(zero.line)
#zero.num<-length(zero.line)

#Main
shift=0
fushion=0
seperation=0
fushion.shift=0
seperation.shift=0
total=0

#
shift.mat<-c("chr","start","end")
fushion.mat<-c("chr","start","end")
seperation.mat<-c("chr","start","end")
fushion.shift.mat<-c("chr","start","end")
seperation.shift.mat<-c("chr","start","end")

for (o in 1:num.chromosomes){
  t<-subset(j,chr==sprintf("chr%s",o))
  zero<-t[grep(pattern = "1",t[,7]),]
  zero.line<-row.names(zero)
  zero.line<-as.numeric(zero.line)
  zero.num<-length(zero.line)
  zero.num<-zero.num-1
  size<-tail(t,n=1L)[,3]
  total<-total+size
  for (i in 1:zero.num){
    n<-i+1
    test.mat<-j[zero.line[i]:zero.line[n],]
    nrow.t<-nrow(test.mat)
    i.line<-t[grep(pattern = "1",test.mat[,6]),]
    i.line.num<-as.numeric(row.names(i.line))
    i.m.line<-t[grep(pattern = "-1",test.mat[,6]),]
    i.m.line.num<-as.numeric(row.names(i.m.line))
    min.gap<-min(abs(i.m.line.num-i.line.num))
    start=j[zero.line[i],2]+resolution/2
    end=j[zero.line[n],3]-resolution/2
    dis<-end-start
    sum<-sum(test.mat[,6])
    if (sum<0){
      if (1 %in% test.mat[,6]){
        if (min.gap>=gap.threshold){
          fushion.shift<-fushion.shift+dis
          fushion.shift.mat<-rbind(fushion.shift.mat,c(sprintf("chr%s",o),start,end))
        }
        else {
          fushion<-fushion+dis
          fushion.mat<-rbind(fushion.mat,c(sprintf("chr%s",o),start,end))
        }
      }
      else{ 
        fushion<-fushion+dis
        fushion.mat<-rbind(fushion.mat,c(sprintf("chr%s",o),start,end))
      }
    }
    else if (sum>0){
      if (-1 %in% test.mat[,6]){
        if (min.gap>=gap.threshold){
          seperation.shift<-seperation.shift+dis
          seperation.shift.mat<-rbind(seperation.shift.mat,c(sprintf("chr%s",o),start,end))
        }
        else{
          seperation<-seperation+dis
          seperation.mat<-rbind(seperation.mat,c(sprintf("chr%s",o),start,end))
        }
      }
      else{
        seperation<-seperation+dis
        seperation.mat<-rbind(seperation.mat,c(sprintf("chr%s",o),start,end))
      }
    }
    else if (sum==0){
      if (-1 %in% test.mat[,6]){
        if (min.gap>=gap.threshold){
          shift<-shift+dis
          shift.mat<-rbind(shift.mat,c(sprintf("chr%s",o),start,end))
        }
      }
    }
  }
}
t<-subset(j,chr=="chrX")
zero<-t[grep(pattern = "1",t[,7]),]
#rownames(t) <- seq(1,nrow(t),1)
zero.line<-row.names(zero)
zero.line<-as.numeric(zero.line)
zero.num<-length(zero.line)
zero.num<-zero.num-1
size<-tail(t,n=1L)[,3]
total<-total+size
for (i in 1:zero.num){
  n<-i+1
  test.mat<-j[zero.line[i]:zero.line[n],]
  nrow.t<-nrow(test.mat)
  i.line<-t[grep(pattern = "1",test.mat[,6]),]
  i.line.num<-as.numeric(row.names(i.line))
  i.m.line<-t[grep(pattern = "-1",test.mat[,6]),]
  i.m.line.num<-as.numeric(row.names(i.m.line))
  min.gap<-min(abs(i.m.line.num-i.line.num))
  start=j[zero.line[i],2]+resolution/2
  end=j[zero.line[n],3]-resolution/2
  dis<-end-start
  sum<-sum(test.mat[,6])
  if (sum<0){
    if (1 %in% test.mat[,6]){
      if (min.gap>=gap.threshold){
        fushion.shift<-fushion.shift+dis
        fushion.shift.mat<-rbind(fushion.shift.mat,c("chrX",start,end))
      }
      else {
        fushion<-fushion+dis
        fushion.mat<-rbind(fushion.mat,c("chrX",start,end))
      }
    }
    else{ 
      fushion<-fushion+dis
      fushion.mat<-rbind(fushion.mat,c("chrX",start,end))
    }
  }
  else if (sum>0){
    if (-1 %in% test.mat[,6]){
      if (min.gap>=gap.threshold){
        seperation.shift<-seperation.shift+dis
        seperation.shift.mat<-rbind(seperation.shift.mat,c("chrX",start,end))
      }
      else{
        seperation<-seperation+dis
        seperation.mat<-rbind(seperation.mat,c("chrX",start,end))
      }
    }
    else{
      seperation<-seperation+dis
      seperation.mat<-rbind(seperation.mat,c("chrX",start,end))
    }
  }
  else if (sum==0){
    if (-1 %in% test.mat[,6]){
      if (min.gap>=gap.threshold){
        shift<-shift+dis
        shift.mat<-rbind(shift.mat,c("chrX",start,end))
      }
    }
  }
}

conserved<-total-shift-fushion-seperation-fushion.shift-seperation.shift

cat (sprintf("Conserved region takes %s",conserved/total))
cat("\n")
cat (sprintf("Shifted region takes %s",shift/total))
cat("\n")
cat (sprintf("Fushioned region takes %s",fushion/total))
cat("\n")
cat (sprintf("Seperated region takes %s",seperation/total))
cat("\n")
cat (sprintf("Fushioned&shifted region takes %s",fushion.shift/total))
cat("\n")
cat (sprintf("Seperated$shifted region takes %s",seperation.shift/total))
cat("\n")

write.table(shift.mat,file="shifted domains.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(seperation.mat,file="seperated domains.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(fushion.mat,file="fused domains.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(seperation.shift.mat,file="seperated and shifted domains.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(fushion.shift.mat,file="fused and shifted domains.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

                                             
