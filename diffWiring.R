### Function
#Rewritten differential wiring function, exponentially faster than the previous 
#for each implementation used by Dan.
#
#Takes as input two correlation matrices, like a power=1 adjacency matrix,
#such as [adjacency(X,power = 1)].  Examines all gene pairs with a change in
#correlation above deltaThresh, estimates significance based on the number of
#samples in each grouping and carries forward with any with p_fdr < pThresh.
#
#The number of significantly different pairings for each gene AKA the number of
#differential edges is computed, then significant is again computed using a
#binomial test where p is the overall proportion of detected significantly
#differential expressed edges.  A p_fdr test is again computed and a tibble
#is output containing each gene, the number of sig. DE, the genes corresponding
#to those edges and the p-values.
#
#General needs:
#corrA/B need to be named matrices/dataframes with geneinfo in the names, 
#and they need to be sorted the same as the geneInfo
diffWiring <- function(corrA,
                       corrB,
                       nA,
                       nB,
                       pThresh = 0.05,
                       deltaThresh = 0.5,
                       geneInfo = Region_Info) #how different are the 2 correlations - threshold of what is a significant change
                       {
  #Tests
  stopifnot(require(plyr))
  stopifnot(require(psych))
  stopifnot(require(tidyverse))
  stopifnot(require(dplyr))
  stopifnot(nrow(corrA) == nrow(corrB))
  stopifnot(ncol(corrA) == ncol(corrB))
  stopifnot(nrow(corrA) == nrow(geneInfo))
  stopifnot(all(rownames(corrA) == rownames(corrB)))
  stopifnot(all(colnames(corrA) == colnames(corrB)))
  stopifnot(nA > 0)
  stopifnot(nB > 0)
  stopifnot(pThresh < 1)
  stopifnot(deltaThresh < 1)
 
  #try to figure out which column of the geneInfo is associated with corrA/B
  geneID_A <- sapply(geneInfo,
                     function(x, y) sum(y %in% x), y = rownames(corrA)) %>%
    which.max %>% names
  geneID_B <- sapply(geneInfo,
                     function(x, y) sum(y %in% x), y = rownames(corrB)) %>%
    which.max %>% names
  
  #Make sure A and B are the same ones, and they are in order
  stopifnot(geneID_A == geneID_B)
  stopifnot(all(rownames(corrA) == geneInfo %>% pull({{geneID_A}})))
  stopifnot(all(colnames(corrA) == geneInfo %>% pull({{geneID_A}})))
  
  #Find gene-pairs with changes in correlation greater than deltaThresh between
  #the two correlation matrices. Note this is a symmetric matrix, so we have
  #redundant information/pairs in the upper/lower triangular matricies. We will
  #just keep the upper (sparse matrix)
  changedGenePairs <- abs(corrA - corrB) > deltaThresh
  geneIndices      <- which(changedGenePairs, arr.ind = T, useNames = T) #Only pulls gene/regions where this is true into a matrix
  geneIndices      <- geneIndices[geneIndices[,1] < geneIndices[,2], ] 
  
  
  #This is a "sparse" matrix with row/column info and the correlation from
  #corrA and corrB.  We only keep the upper triangular matrix here to minimize
  #the number of statistical tests we run
  diffEdges <- data.frame(row = geneIndices[, 1],
                          col = geneIndices[, 2],
                          A = corrA[geneIndices],
                          B = corrB[geneIndices]) 
  
  #Erase unused big things
  rm(geneIndices, changedGenePairs)
  
  #Wrapper for r.test to use map2
  corrSigTest <- function(A, B, n, n2) {
    return(r.test(n = n, r12 = A, r34 = B, n2 = n2)$p)} 
  
  diffEdges$p <- 
    map2_dbl(diffEdges$A, diffEdges$B, 
                    corrSigTest,
                    n = nA,
                    n2 = nB)
  
  #FDR correction to edge detection and filtering. Note: depending on our
  #nA/nB and deltaThresh this may not actually do anything.
  diffEdges$p_adj <- p.adjust(diffEdges$p, method="fdr")
  diffEdges <- diffEdges %>% filter(p_adj < pThresh)
  
  #Relevant stats for binomial tests
  affectedEdges   <- nrow(diffEdges)
  totalGenes      <- nrow(geneInfo)
  totalEdges      <- totalGenes^2
  edgeChangeRate  <- affectedEdges / totalEdges
  
  stopifnot(affectedEdges > 0)
  
  #Add back in the lower matrix, just stack rows after swapping row/col indices
  #We do this because we want to do stats on each gene, not on each edge. even
  #though edges are not directional here.
  diffEdges <-
    bind_rows(
      diffEdges,
      diffEdges %>% 
        rename(col = row,
               row = col)
    )
  
  #Count the number of edges significantly changed for each gene, this is the
  #equivalent of a rowsum or colsum on our sparse matrix.  This is why we added
  #back in the lower triangular matrix, because this is no longer symmetric
  DW <- diffEdges %>%
    group_by(row) %>%
    dplyr::summarize(n = n(),
              .groups="rowwise") %>%
    ungroup()
  
  #This applies a simple binomial test to test if the differential edging is 
  #significant based on its frequency. A wrapper for binom.test for map. x is number of differential edges.
  #n - total number of edges. 
  DWtest <- function(x, n , p, alternative) {
    return(binom.test(x, n, p , alternative)$p.value)}
  
  DW$p <- map_dbl(DW$n,
                         DWtest,
                         n = totalGenes,
                         p = edgeChangeRate,
                         alternative = "g")
  #FDR correction
  DW$p.adj <- p.adjust(DW$p, method = "fdr")
  
  #Finds all the column numbers matching a given row, then pulls the external 
  #gene names corresponding to them as well as the total correlation in each
  #correlation matrix so we can find a logFC. 
  
  edgeGenes <- unstack(diffEdges[,c(1,2)]) %>%
    enframe(name = "row",
           value = "genes") %>%
    mutate(row = as.numeric(row)) %>%
    rowwise() %>%
    mutate(
      corrA = sum(abs(corrA[row, genes])),
      corrB = sum(abs(corrB[row, genes])),
      logFC = log2((corrA+0.1)/(corrB=0.1)), #AC added 0.1 to values to prevent infinities
      regions = sort(geneInfo[genes, "name"]) %>%
        paste(collapse =", ")) %>%
    ungroup()
  
  DW <- DW %>%
    left_join(edgeGenes, by = "row")
  
  #Erase unused big things
  rm(corrA, corrB, edgeGenes, diffEdges)
  
  #Change rows to ensembl_gene_id
  DW <- DW %>%
    mutate(name = geneInfo[row, "name"])
  
  #Convert to full gene list, fill in NA
  
  DW <- DW %>% 
    replace_na(list(n = 0, p = 1, p.adj = 1, regions = "",
                    corrA = 0, corrB = 0, logFC = 0))
  
  return(DW)

}
