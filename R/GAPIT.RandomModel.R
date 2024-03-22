`GAPIT.RandomModel` <-
function(GWAS,Y,CV=NULL,X,cutOff=0.01,GT=NULL,name.of.trait=NULL,N.sig=NULL,n_ran=500,ld.cut=FALSE){
    #Object: To calculate the genetics variance ratio of Candidate Genes
    #Output: The genetics variance raito between CG and total
    #Authors: Jiabo Wang and Zhiwu Zhang
    # Last update: Nov 6, 2019
    ##############################################################################################
    if(!require(lme4))  install.packages("lme4")
    library("lme4")
    print("GAPIT.RandomModel beginning...")
    if(is.null(GT))GT=as.character(Y[,1])
    # name.of.trait=colnames(Y)[2]
    # GWAS=GWAS[order(GWAS[,3]),]
    # GWAS=GWAS[order(GWAS[,2]),]
    P.value=as.numeric(GWAS[,4])
    P.value[is.na(P.value)]=1
    if(is.null(N.sig))
    {
    cutoff=cutOff/nrow(GWAS)
    index=P.value<cutoff    
    }else{
    sort.p=sort(P.value)
    print("GAPIT setup Number of significant markers into Random model:")
    print(N.sig)
    cutoff=max(sort.p[1:N.sig])
    index=P.value<=cutoff 
    }
    
    if(length(unique(index))==1)
    {
    	print("There is no significant marker for VE !!")
    	return(list(GVs=NULL))
    }
    geneGD=X[,index,drop=FALSE]
    geneGWAS=GWAS[index,,drop=FALSE]
    var.gd=diag(var(geneGD))
    var.index=var.gd>0.0001
    geneGD=geneGD[,var.index,drop=FALSE]
    geneGWAS=geneGWAS[var.index,,drop=FALSE]
    if(ld.cut)
    {
        gene.licols=GAPIT.Licols(X=geneGD)
        geneGD=gene.licols$Xsub
        geneGWAS=geneGWAS[gene.licols$idx,]
    }
    index_T=as.matrix(table(index))
    # print(index_T)
    in_True=ncol(geneGD)

    print(in_True==1)
    if(sum(var.index)==0)
    {
        print("There is no significant marker for VE !!")
        return(list(GVs=NULL))
    }
    if(!is.null(geneGD))
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
        Y=Y[taxa_Y%in%GT,]
        CV=CV[taxa_CV%in%GT,]
    	tree2=cbind(Y,CV[,-1,drop=FALSE],geneGD) # thanks jeremyde 2022.9.14
        }
    }
    if(in_True==1)colnames(tree2)[ncol(tree2)]=paste("gene_",1,sep="")
    n_cv=ncol(CV)-1
    n_gd=in_True
    n_id=nrow(Y)
    
if(!is.null(CV))
{
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

    dflme <- lme4::lmer(command2, data=tree2, control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
     check.nobs.vs.rankZ = "ignore",
     check.nobs.vs.nRE="ignore"))
    gene_names=paste("gene_",1:n_gd,sep="")
    carcor_matrix=as.data.frame(summary(dflme)$varcor)
    carcor_matrix=carcor_matrix[-nrow(carcor_matrix),]
    carcor_matrix=carcor_matrix[match(gene_names,as.character(carcor_matrix[,1])),]
    var_gene=as.numeric(carcor_matrix[,4])
    var_res=as.data.frame(summary(dflme)$varcor)[nrow(as.data.frame(summary(dflme)$varcor)),4]

    print(paste("Candidate Genes could Phenotype_Variance_Explained(%) :",sep=""))
    print(100*var_gene/sum(var_gene,var_res))
    v_rat=100*var_gene/sum(var_gene,var_res)
    # print(dim(geneGWAS))
    # print(length(v_rat))
    gene_list=cbind(geneGWAS,v_rat)
    # print("!!!!")
    # print(gene_list)
    colnames(gene_list)[ncol(gene_list)]="Phenotype_Variance_Explained(%)"
    utils::write.csv(var_gene,paste("GAPIT.Association.Vairance_markers.", name.of.trait,".csv",sep=""),quote = FALSE,  row.names = FALSE)
    utils::write.csv(gene_list,paste("GAPIT.Association.PVE.", name.of.trait,".csv",sep=""),quote = FALSE,  row.names = FALSE)
    colnames(gene_list)[ncol(gene_list)]="Variance_Explained"
    colnames(gene_list)[which(colnames(gene_list)%in%c("maf","MAF"))]="MAF"

if(sum(is.na(gene_list[1,c(4:8)]))==0)
{
     
        gene_list=gene_list[order(as.numeric(gene_list$effect)),]

    if(n_gd>=5)
        {
        n=10
        do_color = grDevices::colorRampPalette(c("green", "red"))(n)
            # graphics::par(mar=c(4,5,4,4),cex=1)
            x=as.numeric(gene_list$MAF)
            if(min(x)<0)
            {
                print("The MAF present negative values!!!")
                print("GAPIT will not output PVE against MAF plots!!!")
                return(list(GVs=var_gene/sum(var_gene+var_res),PVEs=gene_list))
            } 
            y=as.numeric(gene_list$effect)
            x.lim=max(x)+max(x)/10
            y.lim=max(y)+max(y)/10
            z=gene_list$Variance_Explained
            quantile_cut = stats::quantile(z)
            r2_color=rep("black",n_gd)
        for(i in 1:(n/2))
        {
            r2_color[z<=quantile_cut[i+1]&z>=quantile_cut[i]]=do_color[2*i]
        }
            
            print("Creating marker p-value, MAF, estimated effect, PVE 3 plot...")

            grDevices::pdf(paste("GAPIT.Association.Significant_SNPs.", name.of.trait,".pdf" ,sep = ""), width =10, height = 3.5)      
            layout.matrix <- matrix(c(1,2,3), nrow = 1, ncol = 3)
            layout(mat = layout.matrix,
                   heights = c(100,80,120), # Heights of the two rows
                   widths = c(2, 2,2)) # Widths of the two columns
            par(mar = c(5, 5, 2, 1))
            # print(head(gene_list))
            # print(length(gene_list$maf))
            # print(length(gene_list$P.value))
            plot(gene_list$MAF,-log10(gene_list$P.value),xlab="MAF",las=1,
            cex=1.2,xlim =c(0,x.lim) ,main="a",
            ylab=expression(-log[10](italic(p))))
            # par(mar = c(5, 5, 2, 1))
            # print(min(y))
            # print(max(y))
            plot(gene_list$MAF,gene_list$effect,cex=1.2,main="b",
            xlab="MAF",ylim=c(min(y), max(y)), xlim =c(0,x.lim) ,las=1,
            ylab="Estimated Effect")
            # par(mar = c(5, 5, 2, 1))
            plot(gene_list$MAF,gene_list$Variance_Explained,cex=1.2,las=1,
            xlab="MAF",xlim =c(0,x.lim) ,main="c",
            ylab="Phenotypic Variance Explained (%)")
            grDevices::dev.off()


        }
}
return(list(GVs=var_gene/sum(var_gene+var_res),PVEs=gene_list))
}#end of GAPIT.RandomModel function
#=============================================================================================
          



