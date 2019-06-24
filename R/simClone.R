# artificial data
simClone <- function(nm, nc, nol=0.05) {

# nc - number of clones
# nm - number of mutations

nmutc=round(nm/nc)
inv=0.5/nc
inv2=0

aVAF=c()
aVAF=abs(rnorm(round(nmutc), mean=0.5, sd=runif(1,0.05,0.1)))
aVAF[aVAF>1]<-1

for(i in 2:nc)
{
aVAF=c(aVAF, rnorm(round(nmutc), mean=runif(1,inv2,inv2+inv), sd=runif(1,0.01,0.1)))
inv2=inv2+inv
}

#add one clone after 0.5
simVAF=abs(aVAF)

##noise to the data
if(nol==0)
{0
  simVAF=abs(aVAF)
}else{

  simVAF=abs(aVAF)
  simVAF[sample(1:nm, nm*nol)]<-rkumar(round(nm*nol), 0.5, 0.5)
}

names<-paste0(rep("mut", length(simVAF)), 1:length(simVAF))
res=data.frame(names,simVAF)

return(res)
}
