`GAPIT.Memory` <-
function(Memory =NULL,Infor){
#Object: To report memory usage
#Output: Memory 
#Authors: Zhiwu Zhang
# Last update: June 6, 2011 
##############################################################################################
gc()
size <- memory.size()
#print(paste("Memory usage: ",size," for", Infor))
if(is.null(Memory)) {
Increased=0
Memory =cbind(Infor,size ,Increased)
}else{
Increased=0
Memory.current=cbind(Infor,size ,Increased)
Memory=rbind(Memory,Memory.current)
Memory[nrow(Memory),3]=as.numeric(as.matrix(Memory[nrow(Memory),2]))-as.numeric(as.matrix(Memory[nrow(Memory)-1,2]))
}

return (Memory)
}#end of GAPIT.Memory function
#=============================================================================================

