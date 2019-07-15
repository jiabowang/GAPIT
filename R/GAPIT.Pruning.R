`GAPIT.Pruning` <-
function(values,DPP=5000){
#Object: To get index of subset that evenly distribute
#Output: Index
#Authors: Zhiwu Zhang
# Last update: May 28, 2011 
##############################################################################################
#No change if below the requirement
if(length(values)<=DPP)return(c(1:length(values)))
  
#values= log.P.values
values=sqrt(values)  #This shift the weight a little bit to the low building.

#Handler of bias plot
rv=runif(length(values))
values=values+rv
values=values[order(values,decreasing = T)]

theMin=min(values)
theMax=max(values)
range=theMax-theMin
interval=range/DPP

ladder=round(values/interval)
ladder2=c(ladder[-1],0)
keep=ladder-ladder2
index=which(keep>0)


return(index)
}#end of GAPIT.Pruning 
#=============================================================================================
