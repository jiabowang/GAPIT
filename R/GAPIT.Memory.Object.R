`GAPIT.Memory.Object` <-
function(name.of.trait="Trait"){
# Object: To report memoery usage
# Authors: Heuristic Andrew
# http://heuristically.wordpress.com/2010/01/04/r-memory-usage-statistics-variable/
# Modified by Zhiwu Zhang
# Last update: may 29, 2011 
############################################################################################## 
# print aggregate memory usage statistics 
print(paste('R is using', memory.size(), 'MB out of limit', memory.limit(), 'MB')) 
  
# create function to return matrix of memory consumption 
object.sizes <- function() 
{ 
    return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
        object.size(get(object.name)))))) 
} 

# export file in table format 
memory=object.sizes() 
file=paste("GAPIT.", name.of.trait,".Memory.Object.csv" ,sep = "")
write.table(memory, file, quote = FALSE, sep = ",", row.names = TRUE,col.names = TRUE)


# export file in PDF format 
pdf(paste("GAPIT.", name.of.trait,".Memory.Object.pdf" ,sep = ""))
# draw bar plot 
barplot(object.sizes(), 
    main="Memory usage by object", ylab="Bytes", xlab="Variable name", 
    col=heat.colors(length(object.sizes()))) 
# draw dot chart 
dotchart(object.sizes(), main="Memory usage by object", xlab="Bytes") 
# draw pie chart 
pie(object.sizes(), main="Memory usage by object")
dev.off()  
}
#=============================================================================================

