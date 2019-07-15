`GAPIT.RemoveDuplicate` <-
function(Y){
#Object: NA
#Output: NA
#Authors: Zhiwu Zhang 
# Last update: Augus 30, 2011 
##############################################################################################
return (Y[match(unique(Y[,1]), Y[,1], nomatch = 0), ] )
}
#=============================================================================================

