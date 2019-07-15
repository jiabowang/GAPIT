`GAPIT.Create.Indicator` <-
function(xs, SNP.impute = "Major" ){
#Object: To esimate variance component by using EMMA algorithm and perform GWAS with P3D/EMMAx
#Output: ps, REMLs, stats, dfs, vgs, ves, BLUP,  BLUP_Plus_Mean, PEV
#Authors: Alex Lipka and Zhiwu Zhang
# Last update: April 30, 2012
##############################################################################################
#Determine the number of bits of the genotype

bit=nchar(as.character(xs[1]))


#Identify the SNPs classified as missing

if(bit==1)  {
xss[xss=="xs"]="N"
xs[xs=="-"]="N"
xs[xs=="+"]="N"
xs[xs=="/"]="N"
xs[xs=="K"]="Z" #K (for GT genotype)is is replaced by Z to ensure heterozygose has the largest value
}

if(bit==2)  {
xs[xs=="xsxs"]="N"
xs[xs=="--"]="N"
xs[xs=="++"]="N"
xs[xs=="//"]="N"
xs[xs=="NN"]="N"
}

#Create the indicators

#Sort the SNPs by genotype frequency
xs.temp <- xs[-which(xs == "N")]

frequ<- NULL
for(i in 1:length(unique(xs.temp))) frequ <- c(frequ, length(which(xs == unique(xs)[i])))

unique.sorted <- cbind(unique(xs.temp), frequ)

print("unique.sorted is")
print(unique.sorted)

unique.sorted <- unique.sorted[order(unique.sorted[,2]),]
unique.sorted <- unique.sorted[,-2]


#Impute based on the major and minor allele frequencies
if(SNP.impute == "Major") xs[which(is.na(xs))] = unique.sorted[1]
if(SNP.impute == "Minor") xs[which(is.na(xs))] = unique.sorted[length(unique.sorted)]
if(SNP.impute == "Middle") xs[which(is.na(xs))] = unique.sorted[2]
x.ind <- NULL
for(i in unique.sorted){
 x.col <- rep(NA, length(xs))
 x.col[which(xs==i)] <- 1
 x.col[which(xs!=i)] <- 0
 x.ind <- cbind(x.ind,x.col)                         
}



return(x.ind)

print("GAPIT.Create.Indicator accomplished successfully!")
}#end of GAPIT.Create.Indicator function
#=============================================================================================

