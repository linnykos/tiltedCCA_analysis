# from nancy
# Assumes that seurat.obj has an assay "RNA" for gene expression,
# and an assay "ATAC" for ATAC.
getPeaksForGenes=function(seurat.obj, gene.names, get.nearby.genes=FALSE,
                          include.cell.features=NULL, ext.upstream=50000, ext.downstream=50000){
  DefaultAssay(seurat.obj) <- "RNA"
  Y <- Seurat::FetchData(seurat.obj, vars=gene.names, slot="data")
  peaks.gr <- GenomicRanges::granges(seurat.obj@assays$ATAC)
  
  n.genes <- length(gene.names)
  Seurat::DefaultAssay(seurat.obj)="ATAC"
  X <- vector("list", n.genes)
  cat("Getting peaks for genes ...\n")
  annotations <- seurat.obj@assays$ATAC@annotation
  genes.near <- rep("",0)
  
  TSS.all <- ifelse(BiocGenerics::strand(annotations)=="+", BiocGenerics::start(annotations), BiocGenerics::end(annotations))
  TSS.all <- GenomicRanges::GRanges(seqnames = seqnames(annotations), 
                                    ranges = IRanges::IRanges(start=TSS.all, width=1), 
                                    strand = BiocGenerics::strand(annotations), 
                                    gene_name = annotations$gene_name)
  for(g in 1:n.genes){
    cat("Gene ",g," of ", n.genes,"... \n" )
    temp=Signac::LookupGeneCoords(seurat.obj, gene.names[g])
    if(is.null(temp)) next
    GenomeInfoDb::seqlevelsStyle(temp)="UCSC"  # just so that the seqlevelsStyle matches up with peaks.gr
    temp2=annotations[which(annotations$gene_name==gene.names[g] & annotations$type=="cds")]
    if(length(temp2)==0) next
    gene.strand=BiocGenerics::strand(temp2)[1]
    TSS.position <- ifelse(gene.strand == "+", BiocGenerics::start(temp), BiocGenerics::end(temp))
    TSS <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(temp),
                   ranges = IRanges::IRanges(start = TSS.position, width = 1),
                   strand = BiocGenerics::strand(temp),
                   gene_name = gene.names[g])
    upstream.gr = getRegion(c(ext.upstream, 0), TSS)
    downstream.gr = getRegion(c(0, ext.downstream), TSS)
    
    # populate the peaks counts array.
    peaks.upstream = GenomicRanges::countOverlaps(peaks.gr, upstream.gr)>0
    peaks.downstream = GenomicRanges::countOverlaps(peaks.gr, downstream.gr)>0
    
    counts1 <- t(as.matrix(seurat.obj@assays$ATAC@counts[peaks.upstream,,drop=F]))# counts within these peaks.
    counts2 <- t(as.matrix(seurat.obj@assays$ATAC@counts[peaks.downstream,,drop=F])) # counts within these peaks.
    
    X[[g]]$counts = cbind(counts1,counts2)
    X[[g]]$peaks.gr = c(peaks.gr[peaks.upstream], peaks.gr[peaks.downstream])
    TSS.subset = TSS.all[GenomicRanges::countOverlaps(TSS.all, getRegion(c(2*ext.upstream, 2*ext.downstream), TSS))>0]
    X[[g]]$genes.near.peaks = getGenesNearPeaks(gr=X[[g]]$peaks.gr, TSS.gr=TSS.subset, exclude.gene=gene.names[g], ext.downstream=ext.downstream, ext.upstream=ext.upstream)
    X[[g]]$isDownstream = c(rep(FALSE,sum(peaks.upstream)), rep(TRUE, sum(peaks.downstream)))
    genes.near = c(genes.near, unique(unlist(X[[g]]$genes.near.peaks)))
  }
  
  
  Seurat::DefaultAssay(seurat.obj)="RNA"
  Y.near.cell=NULL
  if(get.nearby.genes) Y.near.cell = Seurat::FetchData(seurat.obj, vars=genes.near)
  
  Seurat::DefaultAssay(seurat.obj)="RNA"
  if(all(is.null(include.cell.features))){
    include.cell.features=c("nCount_RNA", "nCount_ATAC")
  }
  features.cell=Seurat::FetchData(seurat.obj, vars=include.cell.features)
  obj=list(X.cell=X, Y.cell=Y, Y.near.cell=Y.near.cell, features.cell=features.cell, gene.names=gene.names, 
           ext.downstream=ext.downstream, ext.upstream=ext.upstream)
  obj
}

getRegion <- function(vicinity.size, TSS){
  upstream.size <- vicinity.size[1]
  downstream.size <- vicinity.size[2]
  if (upstream.size == 0) {
    upstream <- TSS
  } else if (upstream.size > 0) {
    upstream <- GenomicRanges::flank(TSS, width=upstream.size, start=T) # upstream excluding the TSS
  }
  if (downstream.size == 0) {
    downstream <- TSS
  } else if (downstream.size > 0) {
    downstream <- GenomicRanges::flank(TSS, width=downstream.size, start=F) # downstream excluding the TSS
  }
  vicinity <- GenomicRanges::punion(upstream, downstream, fill.gap=T)
  BiocGenerics::start(vicinity) <- pmax(0, BiocGenerics::start(vicinity))
  # end(vicinity) <- pmin(end(vicinity), choromosome.size))
  vicinity$gene_name <- TSS$gene_name
  return(vicinity)
}

getGenesNearPeaks = function(gr, TSS.gr, exclude.gene=NULL, 
                             ext.downstream=5e5, ext.upstream=5e5){
  cat("Getting genes near ", length(gr), " peaks.\n", sep="")
  BiocGenerics::start(gr)= BiocGenerics::start(gr)-ext.upstream
  BiocGenerics::end(gr) = BiocGenerics::end(gr)+ext.downstream
  genes.near.peaks = vector("list", length(gr))
  for(i in 1:length(gr)){
    genes.near.peaks[[i]]=setdiff(unique(TSS.gr$gene_name[GenomicRanges::countOverlaps(TSS.gr,gr[i])>0]), exclude.gene)  
  }
  genes.near.peaks
}


#####################3

aggregateCells = function(obj, aggregate.over){
  cat("Aggregating X and Y matrices over group: ", aggregate.over," ... \n")
  
  aggr.ind=which(colnames(obj$features.cell)==aggregate.over)
  group.names=unique(obj$features.cell[,aggr.ind])
  group=obj$features.cell[,aggr.ind]
  X.aggr=initializeX(obj$X.cell, group.names)
  Y.aggr=matrix(nrow=length(group.names), ncol=ncol(obj$Y.cell))
  Y.near.aggr=NULL
  if(!is.null(obj$Y.near.cell)) Y.near.aggr=matrix(nrow=length(group.names), ncol=ncol(obj$Y.near.cell))
  
  metacell.barcodes=vector("list", length(group.names))
  cell.barcodes=row.names(obj$features.cell)
  for(i in 1:length(group.names)){
    cat(i, " ")
    which.cells=which(group==group.names[i])
    metacell.barcodes[[i]]=cell.barcodes[which.cells]
    Y.aggr[i,]=apply(obj$Y.cell[which.cells,,drop=F],2,mean)
    if(!is.null(obj$Y.near.cell)) Y.near.aggr[i,]=apply(obj$Y.near.cell[which.cells,,drop=F],2,mean)
    
    for(j in 1:length(obj$X.cell)){
      if(length(obj$X.cell[[j]])>0){
        X.aggr[[j]]$counts[i,]=apply(obj$X.cell[[j]]$counts[which.cells,, drop=F], 2, mean)
      }
    }
  }
  colnames(Y.near.aggr)=colnames(obj$Y.near.cell)
  
  cat("\n")
  obj$metacell.barcodes=metacell.barcodes
  obj$X.aggr=X.aggr
  obj$Y.aggr=Y.aggr
  obj$Y.near.aggr=Y.near.aggr
  obj$aggregate.over=aggregate.over
  obj
}

initializeX = function(X, group.names){
  newX=vector("list", length=length(X))
  for(i in 1:length(X)){
    if(length(X[[i]])>0){
      newX[[i]]$counts=matrix(nrow=length(group.names), ncol=ncol(X[[i]]$counts), data=NA)
      rownames(newX[[i]]$counts)=group.names
      colnames(newX[[i]]$counts)=colnames(X[[i]]$counts)
      newX[[i]]$peaks.gr=X[[i]]$peaks.gr
      newX[[i]]$isDownstream = X[[i]]$isDownstream
      newX[[i]]$genes.near.peaks = X[[i]]$genes.near.peaks
    }
  }
  newX
}

#####################3

# stores the pvalues in the granges object "peaks.gr" under the appropriate X list
# if using aggregate, then X=X.aggr, otherwise doing per-cell correlations, then X=X.cell.
findLinks <- function(obj, use.aggregate=TRUE, method="spearman"){
  cat("Finding links, use.aggregate=", use.aggregate,", method=", method,"\n")
  if(use.aggregate){
    X=obj$X.aggr
    Y=obj$Y.aggr
  } else {
    X=obj$X.cell
    Y=obj$Y.cell
  }
  ngenes=length(X)
  
  for(g in 1:ngenes){
    cat("Gene",g,"... ")
    y=Y[,g]
    if(length(X[[g]])>0){
      X[[g]]$peaks.gr$pval.spearman=rep(NA,length(X[[g]]$peaks.gr))      
      for(j in 1:ncol(X[[g]]$counts)){
        x = X[[g]]$counts[,j]
        # can also add covariates here "u"
        
        if(method=="spearman"){
          res = stats::cor.test(x,y, method="spearman")
          X[[g]]$peaks.gr$pval.spearman[j]=res$p.value
          X[[g]]$peaks.gr$estimate[j]=res$estimate
          
        }
      }
    }
  }
  cat("\n")
  if(use.aggregate){
    obj$X.aggr = X
  } else {
    obj$X.cell = X
  }
  obj
}

###################

# simply aggregates all of the peaks.gr Granges object across genes into one
# big Granges object.  Does not do any filtering.
collapseLinks <- function(obj, use.aggregate=TRUE){
  if(use.aggregate){
    X=obj$X.aggr
    Y=obj$Y.aggr
  } else {
    X=obj$X.cell
    Y=obj$Y.cell
  }
  
  peaks.all=NULL
  for(i in 1:length(X)){
    temp=X[[i]]$peaks.gr
    if(length(temp)>0){
      temp$target.gene=obj$gene.names[i]
      if(is.null(peaks.all)){
        peaks.all=temp    
      } else {
        peaks.all=c(peaks.all,temp)
      }
    }
  }
  peaks.all
}

