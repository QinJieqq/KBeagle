`KBeagle` <- function(
  GD=NULL,
  #population_sample=NULL,
  cluster_sample=NULL,
  file.output=FALSE
  ){
    #Object: Genotype file with missing values (format 0,1,2)  
    #Input: cluster_sample - The number of samples the user wants to include in each cluster  
    #Input: GD - Genotype file to be imputed  
    #Output: imputed file (without missing values)
    #Authors: QinJie
    #Start  date: November 17, 2022	



#Fill the missing genotype into "1" to get the cluster file
#supports package
#baegle.jar
print("beagle is being downloaded …………")
download.file("https://raw.githubusercontent.com/liu-xinrui/data/main/beagle.jar",destfile = "beagle.jar" , mode='wb')
print(paste("beagle supports java package downloads to ",getwd(),sep=""))
#data.table
library(data.table)
library(factoextra)
library(ggplot2)


#main
if(colnames(GD)[2]=="IID"){
	GD=GD[,-c(1,3:6)]
}
X=GD[,-1]
n=nrow(X)
m=ncol(X)
index.m=as.matrix(X)
index0=X==0
index1=X==1
index2=X==2
indexna=is.na(X)
X2=as.matrix(X)
X2[index0]="0"
X2[index1]="1"
X2[index2]="2"
X2[indexna]="1"
X3=matrix(as.numeric(X2),n,m)
X4=data.frame(GD[,1],X3)
GM=X4


NA_data=GD[,-1]
Cluster_data=GM[,-1]
#Kmeans cluster
c1=scale(Cluster_data)
aa=fviz_nbclust(c1,kmeans,method="wss",k.max=round(nrow(Cluster_data)/cluster_sample))
a0=aa$data$y
a1=a0[-1]
a2=a0[-length(a0)]
dv=a2-a1
dv[dv<=0]=abs(dv[dv<=0])+sum(dv[dv>=0])
kn=which(dv== min(dv), arr.ind = TRUE)
KN=kn+1

#Extract the individual ID number in the population(IID)
data=data.frame(IID=GD[,1],NA_data)
data1=data.frame(IID=GD[,1],Cluster_data)  


#According to the results of K-Means clustering, individuals in the same cluster are divided into a data set, and each data set is imputed with Beagle software
library(cluster)
data.kmeans<-kmeans(Cluster_data,KN)
t=table(data$IID,data.kmeans$cluster)
colnames(t)=paste("type",1:KN,sep="")
t1=cbind(data.frame(GD[,1]),apply(t,2,as.numeric))
colnames(t1)[1]="IID"
Z.bgl2=NULL
for(i in 1:KN)
{
  d=t1[grep("1",t1[,i+1]),]
  type=merge(d,data,all=FALSE,by="IID")
  type_result=type[,-c(2:(KN+1))]
  Y=type_result[,-1]
  index0=Y==0
  index1=Y==1
  index2=Y==2
  indexna=is.na(Y)
  Y[index0]="A\tA"
  Y[index1]="A\tB"
  Y[index2]="B\tB"
  Y[indexna]="?\t?"
  myGD=cbind("M",type_result[,1],Y)
  write.table(myGD,file="yak.bgl",quote=F,sep="\t",col.name=F,row.name=F)
  system("java -Xmx50g -jar beagle.jar unphased=yak.bgl missing=? out=test" )
  genotype.full <- read.delim("test.yak.bgl.phased.gz",sep=" ",head=T)
  genotype.c=as.matrix(genotype.full[,-(1:2)])
  index.A=genotype.c=="A"
  index.B=genotype.c=="B"
  nr=nrow(genotype.c)
  nc=ncol(genotype.c)
  genotype.n=matrix(0,nr,nc)
  genotype.n[index.A]=0
  genotype.n[index.B]=1
  n2=ncol(genotype.n)
  odd=seq(1,n2-1,2)
  even=seq(2,n2,2)
  g0=genotype.n[,odd]
  g1=genotype.n[,even]
  Z.bgl=g0+g1
  Z.bgl=cbind(type_result[,1],Z.bgl)
  colnames(Z.bgl)=as.character(c("IID",colnames(type_result)[-1]))
  Z.bgl2=rbind(Z.bgl2,Z.bgl)
}
X5=data.frame(Z.bgl2)
loc=match(data$IID,X5$IID)
W=Z.bgl2[loc,]
imputed_file=W
if(file.output==TRUE){}
print("The output file is imputed_file")
print("Deleting a Cached File")
file.remove("beagle.jar")
file.remove(dir()[grep("test*",dir())])
return(list(imputed_file.data=W))
}








