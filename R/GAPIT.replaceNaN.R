`GAPIT.replaceNaN` <-
function(LL) {
#handler of grids with NaN log
#Authors: Zhiwu Zhang
# Last update: may 12, 2011 
##############################################################################################

#handler of grids with NaN log 
index=(LL=="NaN")
if(length(index)>0) theMin=min(LL[!index])
if(length(index)<1) theMin="NaN"
LL[index]=theMin
return(LL)    
}
#=============================================================================================
