`GAPIT.Timmer` <-
function(Timmer=NULL,Infor){
#Object: To report current time
#Output: Timmer
#Authors: Zhiwu Zhang
# Last update: may 8, 2011 
##############################################################################################
Time<- Sys.time()
if(is.null(Timmer)) {
Elapsed=0
Timmer=cbind(Infor,Time,Elapsed)
}else{
Elapsed=0
Timmer.current=cbind(Infor,Time,Elapsed)
Timmer=rbind(Timmer,Timmer.current)
Timmer[nrow(Timmer),3]=as.numeric(as.matrix(Timmer[nrow(Timmer),2]))-as.numeric(as.matrix(Timmer[nrow(Timmer)-1,2]))
}

#print(paste('Time used: ', Timmer[nrow(Timmer),3], ' seconds for ',Infor,sep="" )) 
return (Timmer)
}#end of GAPIT.Timmer function
#=============================================================================================

