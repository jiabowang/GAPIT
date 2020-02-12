`GAPIT.RandomModel` <-
function(GWAS,Y,CV=NULL,X,cutOff=0.01,GT=NULL,n_ran=30){
    #Object: To calculate the genetics variance ratio of Candidate Genes
    #Output: The genetics variance raito between CG and total
    #Authors: Jiabo Wang and Zhiwu Zhang
    # Last update: Nov 6, 2019
    ##############################################################################################
    if(!require(lme4))  install.packages("lme4")
    library("lme4")
    #GWAS=myGAPIT_GLM$GWAS
    #CV=myGAPIT_GLM$PCA
    #cut.set=0.01
    #return(list(GVs=NULL))
    print("GAPIT.RandomModel beginning...")
    if(is.null(GT))GT=as.character(Y[,1])
    name.of.trait=colnames(Y)[2]
    cutoff=cutOff/nrow(GWAS)
    P.value=as.numeric(GWAS[,4])
    P.value[is.na(P.value)]=1
    index=P.value<cutoff
    #index[c(1:2)]=TRUE
    
    geneGD=X[,index]
    geneGWAS=GWAS[index,]
    # print(dim(Y))
    # print(dim(CV))
    # print(dim(geneGD))
    if(length(unique(index))==1)
    {
    	print("There is no significant marker for VE !!")
    	return(list(GVs=NULL))

    }
    index_T=as.matrix(table(index))
    in_True=index_T[rownames(index_T)=="TRUE"]
    if(in_True!=1)
    {
    
    	colnames(geneGD)=paste("gene_",1:in_True,sep="")

    }
    colnames(Y)=c("taxa","trait")
    if(is.null(CV))
    {
        if(in_True>n_ran)
        {
    	print("The candidate markers are more than threshold value !")
    	return(list(GVs=NULL))
    	}     	
    	taxa_Y=as.character(Y[,1])
        geneGD=geneGD[GT%in%taxa_Y,]
     # if(!is.null(PC))PC=PC[taxa_GD%in%taxa_Y,]
        Y=Y[taxa_Y%in%GT,]
        tree2=cbind(Y,geneGD)
    	# CV[,2]=1
    }else{
    	if(ncol(CV)==1)
    	{
    		if(in_True+1>n_ran)
            {
    	    print("The candidate markers are more than threshold value !")
    	    return(list(GVs=NULL))
    	    }  
    	taxa_Y=as.character(Y[,1])
        geneGD=geneGD[GT%in%taxa_Y,]
     # if(!is.null(PC))PC=PC[taxa_GD%in%taxa_Y,]
        Y=Y[taxa_Y%in%GT,]
        tree2=cbind(Y,geneGD)
    	}else{
    		if(in_True+ncol(CV)-1>n_ran)
            {
    	    print("The candidate markers are more than threshold value !")
    	    return(list(GVs=NULL))
    	    }
    	colnames(CV)=c("Taxa",paste("CV",1:(ncol(CV)-1),sep=""))
    	taxa_Y=as.character(Y[,1])
    	taxa_CV=as.character(CV[,1])
        geneGD=geneGD[GT%in%taxa_Y,]
     # if(!is.null(PC))PC=PC[taxa_GD%in%taxa_Y,]
        Y=Y[taxa_Y%in%GT,]
        CV=CV[taxa_CV%in%GT,]
    	tree2=cbind(Y,CV[,-1],geneGD)
        }
    }
    if(in_True==1)colnames(tree2)[ncol(tree2)]=paste("gene_",1,sep="")

    	# print(head(tree2))
     #    print(dim(CV))
     #    print(dim(geneGD))
    n_cv=ncol(CV)-1
        # print(n_cv)
    n_gd=in_True
    n_id=nrow(Y)
    # print(n_gd)
    # print(n_cv)
    # print(n_id)
    # if((n_gd+n_cv)>n_ran)
    # {
    # 	print("The candidate markers are more than threshold value !")
    # 	return(list(GVs=NULL))
    # 	} 
if(!is.null(CV))
{
#ff <- paste("trait~1+PC1+PC2+PC3+(1|gene_1)+(1|gene_2)+(1|gene_3)+(1|gene_4)+(1|gene_5)+(1|gene_6)"
#dflme <- lmer(ff, data=tree2)
    if(ncol(CV)==1)
    {
    command0=paste("trait~1",sep="")
    command1=command0
    
    command2=command1
    for(j in 1:n_gd)
{
	command2=paste(command2,"+(1|gene_",j,")",sep="")
}

    }else{
    command0=paste("trait~1",sep="")
    command1=command0
    for(i in 1:n_cv)
{	
	command1=paste(command1,"+CV",i,sep="")
}
    command2=command1
    for(j in 1:n_gd)
{
	command2=paste(command2,"+(1|gene_",j,")",sep="")
}
    }
}else{

    command0=paste("trait~1",sep="")
    command1=command0
    
    command2=command1
    for(j in 1:n_gd)
{
	command2=paste(command2,"+(1|gene_",j,")",sep="")
}

}
#command3=paste(command2,"+(1|gene_",j,")",sep="")
    dflme <- lmer(command2, data=tree2,control=lmerControl(check.nobs.vs.nlev = "ignore",
     check.nobs.vs.rankZ = "ignore",
     check.nobs.vs.nRE="ignore"))

    carcor_matrix=as.data.frame(summary(dflme)$varcor)
    var_gene=as.numeric(carcor_matrix[1:(nrow(carcor_matrix)-1),4])
    var_res=carcor_matrix[nrow(carcor_matrix),4]

    print(paste("Candidate Genes could explain genetics variance :",sep=""))
    print(var_gene/sum(var_gene+var_res))
    v_rat=var_gene/sum(var_gene+var_res)
    gene_list=cbind(geneGWAS,v_rat)
    colnames(gene_list)[ncol(gene_list)]="Variance_Explained"

    write.csv(gene_list,paste("GAPIT.", name.of.trait,".Phenotype_Variance_Explained_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
#gene_list=read.csv("GAPIT.Weight.GrowthIntercept.Phenotype_Variance_Explained_by_Association_Markers.csv",head=T)
if(!is.na(sum(gene_list[1,c(4:8)])))
{
        pdf(paste("GAPIT.", name.of.trait,".Effect_VP.pdf" ,sep = ""), width = 7,height=5.75)
        par(mar=c(4,5,4,4),cex=0.8)

        gene_list=gene_list[order(gene_list$effect),]
        plot(gene_list$effect,gene_list$Variance_Explained,
        	xlab="Estimated Effect",
        	ylab="Variance Explained of Phenotype"
        	)
        dev.off()

        pdf(paste("GAPIT.", name.of.trait,".MAF_VP.pdf" ,sep = ""), width = 7,height=5.75)
        par(mar=c(4,5,4,4),cex=0.8)
        gene_list=gene_list[order(gene_list$maf),]
        plot(gene_list$maf,gene_list$Variance_Explained,xlab="MAF",ylab="Variance Explained of Phenotype")
        dev.off()

    if(n_gd>=10)
        {
        pdf(paste("GAPIT.", name.of.trait,".MAF_Effect_VP.pdf" ,sep = ""), width = 9,height=5.75)
        
        n=10
        layout(matrix(c(1,1,2,1,1,1,1,1,1),3,3,byrow=TRUE), c(2,1), c(1,1), TRUE)
        do_color=colorRampPalette(c("green", "red"))(n)

            par(mar=c(4,5,2,8),cex=0.8)
            y=gene_list$maf
            x=gene_list$effect
            x.lim=max(x)+max(x)/10
            y.lim=max(y)+max(y)/10
            z=gene_list$Variance_Explained
            quantile_cut=quantile(z)
            r2_color=rep("black",n_gd)
        for(i in 1:(n/2))
        {
        	r2_color[z<=quantile_cut[i+1]&z>=quantile_cut[i]]=do_color[2*i]
        }
            plot(y~x,type="p", ylim=c(0,y.lim), xlim = c(min(x), max(x)),col = r2_color, xlab = "",ylab = "", cex.lab=1.2,pch=21,bg=r2_color)
            mtext("Estimated Effect",side=1,line=2.5)
            mtext("MAF",side=2,line=2.5)


            par(mar=c(2,13,5,4),cex=0.5)
            
            barplot(matrix(rep(0.4,times=n),n,1),beside=T,col=do_color,border=do_color,axes=FALSE,horiz =T)
        #legend(x=10,y=2,legend=expression(R^"2"),,lty=0,cex=1.3,bty="n",bg=par("bg"))
            step=length(seq(0,round(max(z),3),by=0.01))
            small_bar=round(seq(0,round(max(z),3),by=(max(z)-min(z))/10),2)
            #main()
            mtext("Variance Explained",side=2,line=0.4,col="black",cex=0.5)

            axis(4,c(1,6,11),c(min(small_bar),median(small_bar),max(small_bar)),las=2,cex.axis = 0.9,tick=F,line=0)
        
        dev.off()
        }
}






return(list(GVs=var_gene/sum(var_gene+var_res)))

}#end of GAPIT.RandomModel function
#=============================================================================================
          



