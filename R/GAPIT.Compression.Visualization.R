#'
#' GAPIT.Compression.Visualization
#' 
#' @description 
#' Compression visualization
#' 
#' @param Compression = Compression,
#' @param name.of.trait = name.of.trait
#' @param file.output = TRUE, should output be automatically written to file.
#'
#' @return 
#' An invisible NULL.
#'
#' @author Alex Lipka and Zhiwu Zhang
#'
#' @export
`GAPIT.Compression.Visualization` <-
  function(Compression = Compression,
           name.of.trait = name.of.trait, 
           file.output = FALSE){
  #Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
  #Output: Three pdfs: One of the log likelihood function, one of the genetic and error variance component,
  #                    and one of the heritabilities
  #Authors: Alex Lipka and Zhiwu Zhang 
  # Last update: May 10, 2011 
  ##############################################################################################
  #Graph the optimum compression 

  print("GAPIT.Compression.Visualization")
  #print(Compression)

  if(length(Compression)<=6) Compression=t(as.matrix(Compression[which(Compression[,4]!="NULL" | Compression[,4]!="NaN"),]))
  if(length(Compression)==6) Compression=matrix(Compression,1,6) 
#print("Compression matrix")
#print(Compression)
#print(length(Compression) )

  if(length(Compression)>6) Compression=Compression[which(Compression[,4]!="NULL" | Compression[,4]!="NaN"),]
  if(length(Compression)<1) return() #no result

#Pie chart for the optimum setting
#-------------------------------------------------------------------------------
  print("Pie chart")
  LL=as.numeric(Compression[,4])
  Compression.best=Compression[1,] 
  variance=as.numeric(Compression.best[5:6])
#colors <- c("grey50","grey70")
  colors <- c("#990000","dimgray")
  varp=variance/sum(variance)
  h2.opt= varp[1]

  labels0 <- round(varp * 100, 1)
  labels <- paste(labels0, "%", sep="")

  legend0=c("Genetic: ","Residual: ")
  legend <- paste(legend0, round(variance*100)/100, sep="")

  LL.best0=as.numeric(Compression.best[4]  )
  LL.best=paste("-2LL: ",floor(LL.best0*100)/100,sep="")
  label.comp=paste(c("Cluster method: ","Group method: ","Group number: "), Compression.best[c(1:3)], sep="")
  theOptimum=c(label.comp,LL.best) 
  #print(variance)
  if( file.output == TRUE ){
    grDevices::pdf(paste("GAPIT.", name.of.trait,".Optimum.pdf", sep = ""), width = 14)
    graphics::par(mfrow = c(1,1), mar = c(1,1,5,5), lab = c(5,5,7))
    graphics::pie(variance,  col=colors, labels=labels,angle=45,border=NA)
    graphics::legend(1.0, 0.5, legend, cex=1.5, bty="n",
                     fill=colors)

    #Display the optimum compression
    graphics::text(1.5,.0, "The optimum compression", col= "gray10")
    for(i in 1:4){
      graphics::text(1.5,-.1*i, theOptimum[i], col= "gray10")
    }
    grDevices::dev.off()
  }

#sort Compression by group number for plot order
Compression=Compression[order(as.numeric(Compression[,3])),]

#Graph compression with multiple groups
#print("Graph compression with multiple groups")


if(length(Compression)==6) return() #For to exit if only one row


#print("It should not go here")

if(length(unique(Compression[,3]))>1)
{
#Create a vector of colors
#print("Setting colors")
color.vector.basic <- c("red","blue","black", "blueviolet","indianred","cadetblue","orange")
color.vector.addition <- setdiff(c(colors()[grep("red",colors())], colors()[grep("blue",colors())]),color.vector.basic )
color.vector.addition.mixed <- sample(color.vector.addition,max(0,((length(unique(Compression[,1])) * length(unique(Compression[,2])))-length(color.vector.basic))))  
color.vector <- c(color.vector.basic,color.vector.addition.mixed )


#Create a vector of numbers for the line dot types
line.vector <-  rep(1:(length(unique(Compression[,1])) * length(unique(Compression[,2]))))

#We want to have a total of three plots, one displaying the likelihood function, one displaying the variance components, and one displaying the
# heritability 
if( file.output == TRUE ){
  grDevices::pdf(paste("GAPIT.", name.of.trait,".Compression.multiple.group", ".pdf", sep = ""), width = 14)
  graphics::par(mfrow = c(2,3), mar = c(5,5,1,1), lab = c(5,5,7))

  # Make the likelihood function plot
  #print("Likelihood")
  k <- 1
  for(i in 1:length(unique(Compression[,1]))){
    for(j in 1:length(unique(Compression[,2]))){
      if((i == 1)&(j == 1)) {
        Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
        x <- as.numeric(Compression.subset[,3])
        y <- as.numeric(Compression.subset[,4])  
        plot(y~x,type="l", pch = 30, lty = line.vector[i], ylim=c(min(as.numeric(Compression[,4])),max(as.numeric(Compression[,4]))), xlim = c(min(as.numeric(Compression[,3])),max(as.numeric(Compression[,3]))),
        col = color.vector[j], xlab = "Number of Groups", ylab = "-2Log Likelihoood",lwd=1 )
        label = paste(c(as.character(unique(Compression[,1]))[k]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
      if((i != 1)|(j != 1)) {
        k <- k+1   
        Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
        x <- as.numeric(Compression.subset[,3])
        y <- as.numeric(Compression.subset[,4])  
        graphics::lines(y~x,type="l", pch = 30, lty = line.vector[i], col = color.vector[j])
        label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }
    }
  }
  #Make a legend
  #legend("topright",  label, fill = color.vector)
  legend.col= 1+floor(length(unique(Compression[,1])) * length(unique(Compression[,2]))/20)
  line.style=rep(1:length(unique(Compression[,1])), each = length(unique(Compression[,2])))      
  line.color=rep(1:length(unique(Compression[,2])), length(unique(Compression[,1])))

  legend("topright",  label, col = color.vector[line.color], lty = line.style, ncol=legend.col,horiz=FALSE,bty="n")
 
  # Make the genetic variance component plots
  #print("genetic variance")
  k <- 1
  for(i in 1:length(unique(Compression[,1]))){
    for(j in 1:length(unique(Compression[,2]))){
      if((i == 1)&(j == 1)) {
        Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
        x <- as.numeric(Compression.subset[,3])
        y <- as.numeric(Compression.subset[,5])  
        plot(y~x,type="l", pch = 17,  lty = line.vector[i], ylim=c(min(as.numeric(Compression[,5])),max(as.numeric(Compression[,5]))), xlim = c(min(as.numeric(Compression[,3])),max(as.numeric(Compression[,3]))),
        col = color.vector[j], xlab = "Number of Groups", ylab = "Genetic Variance", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
      if((i != 1)|(j != 1)) {
        k <- k+1   
        Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
        x <- as.numeric(Compression.subset[,3])
        y <- as.numeric(Compression.subset[,5])  
        graphics::lines(y~x,type="l", pch = 17, lty = line.vector[i], col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }
    }
  }
 #Make a legend
  #legend("topleft",  label, fill = color.vector) 

# Make the residual variance component plots
k <- 1
for(i in 1:length(unique(Compression[,1]))){
  for(j in 1:length(unique(Compression[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,6])  
      plot(y~x,type="l", pch = 17,  ylim=c(min(as.numeric(Compression[,6])),max(as.numeric(Compression[,6]))), xlim = c(min(as.numeric(Compression[,3])),max(as.numeric(Compression[,3]))),
      col = color.vector[j], xlab = "Number of Groups", ylab = "Residual Variance", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression[which( (Compression[,1] == as.character(unique(Compression[,1])[i])) & (Compression[,2] == as.character(unique(Compression[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,6])  
      graphics::lines(y~x,type="l", pch = 17, lty = line.vector[i], col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }  
   }
 }
 #Make a legend
  #legend("topright",  label, fill = color.vector) 


#calculate total variance and h2
#print("h2")
heritablilty.vector <- as.numeric(Compression[,5])/(as.numeric(Compression[,5]) + as.numeric(Compression[,6]))
totalVariance.vector <- as.numeric(as.numeric(Compression[,5]) + as.numeric(Compression[,6]))
Compression.h2 <- cbind(Compression, heritablilty.vector,totalVariance.vector)

# Make the total variance component plots
#print("Total variance")
k <- 1
for(i in 1:length(unique(Compression.h2[,1]))){
  for(j in 1:length(unique(Compression.h2[,2]))){

     if((i == 1)&(j == 1)) {
      Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,8])  
      plot(y~x,type="l", pch = 17,  lty = line.vector[k], ylim=c(min(as.numeric(Compression.h2[,8])),max(as.numeric(Compression.h2[,8]))), xlim = c(min(as.numeric(Compression.h2[,3])),max(as.numeric(Compression.h2[,3]))),
      col = color.vector[1], xlab = "Number of Groups", ylab = "Total Variance", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
    if((i != 1)|(j != 1)) {
      k <- k+1   
      Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
      x <- as.numeric(Compression.subset[,3])
      y <- as.numeric(Compression.subset[,8]) 
      graphics::lines(y~x,type="l", pch = 17, lty = line.vector[i], col = color.vector[j])
      #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }  
    }
  }
  #Make a legend
  #legend("topright",  label, fill = color.vector) 
  

  # Make the heritability plots 
  #print("h2 plot")
  k <- 1
  for(i in 1:length(unique(Compression[,1]))){
    for(j in 1:length(unique(Compression[,2]))){
      if((i == 1)&(j == 1)) {
        Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
        x <- as.numeric(Compression.subset[,3])
        y <- as.numeric(Compression.subset[,7]) 
        plot(y ~ x, type="l", pch = 17,  lty = line.vector[k], 
             ylim=c(min(as.numeric(Compression.h2[,7])), max(as.numeric(Compression.h2[,7]))),
             xlim = c(min(as.numeric(Compression.h2[,3])),max(as.numeric(Compression.h2[,3]))),
             col = color.vector[1], 
             xlab = "Number of Groups", 
             ylab = "Heritability", )
      #label = paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = "")
      }
  
      if((i != 1)|(j != 1)) {
        k <- k+1   
        Compression.subset <- Compression.h2[which( (Compression.h2[,1] == as.character(unique(Compression.h2[,1])[i])) & (Compression.h2[,2] == as.character(unique(Compression.h2[,2])[j]))  ),              ]
        x <- as.numeric(Compression.subset[,3])
        y <- as.numeric(Compression.subset[,7])  
        graphics::lines(y~x,type="l", lty = line.vector[i], pch = 17, col = color.vector[j])
        #label = c(label, paste(c(as.character(unique(Compression[,1]))[i]," ",as.character(unique(Compression[,2]))[j]), collapse = ""))
      }
    }
  }
 
 #Make a legend
  #legend("topleft",  label, fill = color.vector) 
  legend.col= 1+floor(length(unique(Compression[,1])) * length(unique(Compression[,2]))/20)
  line.style=rep(1:length(unique(Compression[,1])), each = length(unique(Compression[,2])))      
  line.color=rep(1:length(unique(Compression[,2])), length(unique(Compression[,1])))

  # Make labels
  plot(0~0,axes=FALSE, type="l", ylab = "", xlab = "", frame.plot=FALSE)
  legend("topleft",  label, col = color.vector[line.color], lty = line.style,
         ncol=legend.col, horiz=FALSE) 
  
  grDevices::dev.off()
}

}#end of Graph compression with multiple groups

#Graph compression with single groups
#print("Graph compression with single groups")
if(length(unique(Compression[,3]))==1& length(unique(Compression[,1]))*length(unique(Compression[,2]))>1)
{

#Graph the compression with only one group
if( file.output == TRUE ){
  grDevices::pdf(paste("GAPIT.Compression.single.group.", name.of.trait, ".pdf", sep = ""),
                 width = 14)
  graphics::par(mfrow = c(2,2), mar = c(5,5,1,1), lab = c(5,5,7))

  nkt=length(unique(Compression[,1]))
  nca=length(unique(Compression[,2]))
  kvr=rep(c(1:nkt),nca)
  kvc0=rep(c(1:nca),nkt)
  kvc=as.numeric(t(matrix(kvc0,nca,nkt)))
  kt.name=Compression[1:nkt,1]

  ca.index=((1:nca)-1)*nkt+1
  ca.name=Compression[ca.index,2]

  KG<- t(tapply(as.numeric(Compression[,4]), list(kvr, kvc), mean))
  colnames(KG)=kt.name
  graphics::barplot(as.matrix(KG),  ylab= "-2 Log Likelihood",beside=TRUE, col = grDevices::rainbow(length(unique(Compression[,2]))))

  KG<- t(tapply(as.numeric(Compression[,5]), list(kvr, kvc), mean))
  colnames(KG)=kt.name
  graphics::barplot(as.matrix(KG),  ylab= "Genetic varaince", beside=TRUE,
                    col = grDevices::rainbow(length(unique(Compression[,2]))))

  KG<- t(tapply(as.numeric(Compression[,6]), list(kvr, kvc), mean))
  colnames(KG)=kt.name
  graphics::barplot(as.matrix(KG),
                    ylab= "Residual varaince", 
                    beside=TRUE, 
                    col=grDevices::rainbow(length(unique(Compression[,2])))
                    )

  KG<- t(tapply(as.numeric(Compression[,5])/(as.numeric(Compression[,5])+as.numeric(Compression[,6])), list(kvr, kvc), mean))
  colnames(KG)=kt.name
  graphics::barplot(as.matrix(KG),  ylab= "Heritability", beside=TRUE, col=grDevices::rainbow(length(unique(Compression[,2]))),ylim=c(0,1))

  graphics::legend("topleft", paste(t(ca.name)), cex=0.8,bty="n", fill=grDevices::rainbow(length(unique(Compression[,2]))),horiz=TRUE)
  grDevices::dev.off()
  }
} #end of Graph compression with single groups

print("GAPIT.Compression.Visualization accomplished successfully!")

#return(list(compression=Compression.h2,h2=h2.opt))
  return(invisible(NULL))

}#GAPIT.Compression.Plots ends here
#=============================================================================================

