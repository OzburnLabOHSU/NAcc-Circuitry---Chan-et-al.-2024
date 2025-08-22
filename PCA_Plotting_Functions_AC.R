#Inputs: 
# gene count matrix (rows = genes, columns = samples)
# factor matrix (rows = samples, columns = factors), expects sample names in
# rownames.
# factor of interest (columns to plot counts by, defaults to libraryprepBatch)
#Ouputs:
#stats matrix
#

plot_stats <- function(geneCounts, 
                       sampleFactors, 
                       factorofInterest = "LibraryPrepBatch",
                       ...
                       ) {
  #Tests
  stopifnot(require(plyr))
  stopifnot(require(tidyverse))
  stopifnot(require(ggpubr))
  stopifnot(require(ggstatsplot))
  stopifnot(is.data.frame(geneCounts))
  stopifnot(is.data.frame(sampleFactors))
  stopifnot(ncol(geneCounts) == nrow(sampleFactors))
  stopifnot(all(colnames(geneCounts) == rownames(sampleFactors)))
  stopifnot(factorofInterest %in% colnames(sampleFactors))

  #Label options for Amy
  extra <- list(...)
  options <- c("xlab", "ylab", "titleL", "titleR")
  if (length(names(extra))) {
    stopifnot(all(names(extra) %in% options))
  }
  if(!length(extra[["xlab"]])) {
    extra[["xlab"]] <- paste("Total counts by",factorofInterest)
  }
  if(!length(extra[["ylab"]])) {
    extra[["ylab"]] <- "Total Gene Counts (millions)"
  }
  if(!length(extra[["titleL"]])) {
    extra[["titleL"]] <- paste("Count by", factorofInterest)
  }
  if(!length(extra[["titleR"]])) {
    extra[["titleR"]] <- "Count summary stats:"
  }  
  
  sampleStats <- data.frame(t(do.call(cbind,
                                      lapply(geneCounts, summary))))
  sampleStats <- sampleStats %>%
    mutate(zeros = colSums(geneCounts == 0),
           percentZeros = zeros / nrow(geneCounts),
           sums = colSums(abs(geneCounts)),
           joinKey = rownames(.)) %>%
    left_join(sampleFactors %>% mutate(joinKey = rownames(.)), by = "joinKey")
  
  rownames(sampleStats) <- sampleStats$joinKey
  
  plotBatch <- 
    ggbetweenstats(
      data = sampleStats %>% mutate(sums = sums / 10^6),
      x = {{factorofInterest}},
      y = sums,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      #var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic()) +
    xlab(extra[["xlab"]]) + ylab(extra[["ylab"]]) 
  
  plotMean <-
    ggbetweenstats(
      data = sampleStats %>%
        mutate(x = "Mean"),
      x = x,
      y = Mean,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic(),
      xlab="",
      point.args = list(
        position = ggplot2::position_jitter(width = 0.15),
        alpha = 0.4, size = 3, stroke = 0, color = 1))
  
  plotMedian <- 
    ggbetweenstats(
      data = sampleStats %>%
        mutate(x = "Median"),
      x = x,
      y = Median,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic(),
      xlab="",
      point.args = list(
        position = ggplot2::position_jitter(width = 0.15),
        alpha = 0.4, size = 3, stroke = 0, color = 2))
  
  plotMinimum <-
    ggbetweenstats(
      data = sampleStats %>%
        mutate(x = "Minimum"),
      x = x,
      y = Min.,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic(),
      xlab="",
      point.args = list(
        position = ggplot2::position_jitter(width = 0.15),
        alpha = 0.4, size = 3, stroke = 0, color = 3))
  
  plotMaximum <-
    ggbetweenstats(
      data = sampleStats %>%
        mutate(x = "Maximum"),
      x = x,
      y = Max.,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic(),
      xlab="",
      point.args = list(
        position = ggplot2::position_jitter(width = 0.15),
        alpha = 0.4, size = 3, stroke = 0, color = 4))
  
  plotZeros <-
    ggbetweenstats(
      data = sampleStats %>%
        mutate(x = "Percent Zeros"),
      x = x,
      y = percentZeros,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic(),
      xlab="",
      point.args = list(
        position = ggplot2::position_jitter(width = 0.15),
        alpha = 0.4, size = 3, stroke = 0, color = 5))
  
  plotSums <-      
    ggbetweenstats(
      data = sampleStats %>%
        mutate(x = "Total"),
      x = x,
      y = sums,
      results.subtitle = FALSE,
      centrality.plotting = FALSE,
      pairwise.comparisons = FALSE,
      outlier.tagging = TRUE,
      outlier.label = joinKey,
      outlier.label.args = list(color = "red"),
      outlier.coef = 2.2,
      var.equal = TRUE,
      plot.type = "violin",
      ggtheme = theme_classic(),
      xlab="",
      point.args = list(
        position = ggplot2::position_jitter(width = 0.15),
        alpha = 0.4, size = 3, stroke = 0, color = 6))
  
  tmp_plot <- wrap_plots( 
    wrap_elements(plotBatch + 
                    plot_annotation(title = extra[["titleL"]])),
    wrap_elements(wrap_plots(plotZeros,
                             plotSums,
                             plotMean,
                             plotMedian,
                             plotMinimum,
                             plotMaximum) +
                    plot_annotation(title = extra[["titleR"]]) &
                    ylab("")), 
    widths = c(1/3,2/3))
  
  return(list(stats = sampleStats,
              plot = tmp_plot))
} #plot_stats

#Inputs: 
# gene count matrix (rows = genes, columns = samples)
# factor matrix (rows = samples, columns = factors), expects sample names in
# rownames.
# factors of interest (columns to test)
# Number of PCs (defaults to all PCs with explained variance above 1%)
#Ouputs:
#stats matrix
plot_PCA <- function(geneCounts,
                     sampleFactors, 
                     factorsofInterest,
                     nPC = 0,
                     plot = FALSE,
                     fdr = "none",
                     pval = 0.05) {

  #Tests
  stopifnot(require(plyr))
  stopifnot(require(tidyverse))
  stopifnot(require(ggpubr))
  stopifnot(require(ggstatsplot))
  stopifnot(require(patchwork))
  stopifnot(require(WRS2))
  stopifnot(require(ComplexHeatmap))
  stopifnot(is.data.frame(geneCounts))
  stopifnot(is.data.frame(sampleFactors))
  stopifnot(ncol(geneCounts) == nrow(sampleFactors))
  stopifnot(all(colnames(geneCounts) == rownames(sampleFactors)))
  stopifnot(length(factorsofInterest)>0)
  stopifnot(all(factorsofInterest %in% colnames(sampleFactors)))

  #Make sure factors of interest ARE factors
  sampleFactors  <- sampleFactors %>%
    mutate(joinKey = rownames(.),
           across({{factorsofInterest}}, factor))
  
  #Make sure we have enough samples in each 
  keepFactors <- c()
  for (tmp_factor in factorsofInterest) {
    if ((sampleFactors %>% pull(!!sym(tmp_factor)) %>% levels %>% length > 1) &
        (sum(sampleFactors %>% 
             count(!!sym(tmp_factor)) %>%      
             pull(n) < 2) == 0)) {
      keepFactors <- c(keepFactors, tmp_factor)
    } else {
      warning(paste0("plot_PCA: Too few samples in factor of interest ",
                     tmp_factor, ".  Excluding from further analysis.\n"),
              call. = FALSE)
    }
  } #end for
  factorsofInterest <- keepFactors
  stopifnot(length(factorsofInterest)>0)

  if(length(factorsofInterest)>2) {
    warning(paste0("plot_PCA: More than 2 factors of interest have been included, 
           note that ONLY two-way interactions will be included in analysis."),
           call. = F)
  }
  
  tmp_pcaRes <- prcomp(t(geneCounts))
  tmp_importance <- summary(tmp_pcaRes)$importance[2, ]
  
  if (nPC == 0) {nPC <- sum(tmp_importance > 0.01)}
  
  #This is exploratory function, we don't want to keep all the data
  tmp_pcaRes$sdev <- tmp_pcaRes$sdev[1:nPC]
  tmp_pcaRes$rotation <- tmp_pcaRes$rotation[, 1:nPC]
  tmp_pcaRes$center <- tmp_pcaRes$center[1:nPC]
  tmp_pcaRes$x <- tmp_pcaRes$x[, 1:nPC]
  
  tmp_plottingData_PCA <- tmp_pcaRes$x[, 1:nPC] %>%
    data.frame %>%
    mutate(joinKey = rownames(.)) %>%
    left_join(sampleFactors %>%
                select(joinKey, all_of(factorsofInterest)), by = "joinKey")
  
  rownames(tmp_plottingData_PCA) <- tmp_plottingData_PCA$joinKey
  
  tmp_plottingData_PCA_long <- tmp_plottingData_PCA %>%
    pivot_longer(!c("joinKey", all_of(factorsofInterest)),
                 names_to="PC",
                 values_to="Score",
                 names_prefix="PC") %>%
    mutate(
      var = tmp_importance[as.numeric(PC)],
      PC = as.factor(as.numeric(PC)),
      x = str_c("PC", PC, "(", round(var*100,digits=2), "%)"),
      x = factor(x,unique(x)))
  
  #We really want to use heteroskedastic anova here, because there's absolutely no reason
  #to think that PCs are going to be normally distributed OR have equal variances for random
  #factors of interest. We could consider sorting PCs by the effect size instead, but that would
  #involve some more work
  
  #One-way's first
  tmp_factor_anova <- sapply(1:nPC, function(pc_i) {
    sapply(factorsofInterest, function(factor_j) {
      tmp <- WRS2::t1way(formula(paste0("PC",pc_i,"~",factor_j)), data = tmp_plottingData_PCA)
      return(data.frame(
        PC = pc_i,
        effect = factor_j,
        p.value = tmp$p.value)
      )
      }, simplify = F) %>%
        bind_rows()
    }, simplify = F) %>% bind_rows()
  
  #Two-way's if we got them
  if(length(factorsofInterest) > 1) {
    tmp_factor_anova_twoway <- 
      sapply(1:nPC, function(pc_i) {
        apply(combn(factorsofInterest,2),2, function(effect) {
          tmp <- WRS2::t2way(formula(paste0("PC",pc_i,"~",paste0(effect, collapse = "+"))),
                      data = tmp_plottingData_PCA)
          return(data.frame(
            PC = pc_i,
            effect = paste0(effect, collapse = ":"),
            p.value = tmp$AB.p.value))
        }) %>% bind_rows()
      }, simplify = F) %>% bind_rows()
    
    tmp_factor_anova <- bind_rows(
      tmp_factor_anova,
      tmp_factor_anova_twoway
    )
    
    #Add combined factors
    tmp_effects <- combn(factorsofInterest,2)
    for(j in 1:ncol(tmp_effects)) {
      tmp_plottingData_PCA <- tmp_plottingData_PCA %>%
        mutate(!!sym(paste0(tmp_effects[,j], collapse = ":")) := paste0(!!sym(tmp_effects[1,j]), "x", !!sym(tmp_effects[2,j])))
      tmp_plottingData_PCA_long <- tmp_plottingData_PCA_long %>%
        mutate(!!sym(paste0(tmp_effects[,j], collapse = ":")) := paste0(!!sym(tmp_effects[1,j]), "x", !!sym(tmp_effects[2,j])))
    }
  }

  #Multiple comparison correction, if called for, and widen
  if( fdr != "none" ) {
  tmp_factor_anova <- tmp_factor_anova %>%
    mutate(p.value = p.adjust(p.value, method = fdr)) 
  }
  
  tmp_factor_anova <- tmp_factor_anova %>%
    pivot_wider(names_from = effect,
                values_from = p.value)
  
  tmp_plotList <- list()
  tmp_plotList2 <- list()
  for (tmp_factor in colnames(tmp_factor_anova)[2:ncol(tmp_factor_anova)]) {

    #Are there any PCs which are significantly associated with this factor?
    tmp <- which(tmp_factor_anova[,tmp_factor] < pval)

    #"Most associated PCs" - this should be changed to effect size someday
    tmp_sigPC <- tmp_factor_anova %>%
      slice_min(order_by = !!sym(tmp_factor), n = 2)
    tmp_sigPC <- tmp_sigPC[1:2,]
    
    tmp_plot <- 
      ggscatterhist(tmp_plottingData_PCA,
                    x = paste0("PC", tmp_sigPC$PC[1]),
                    y = paste0("PC", tmp_sigPC$PC[2]),
                    label = tmp_factor,
                    xlab = paste0("PC",
                                  tmp_sigPC$PC[1],
                                  "(",
                                  round(tmp_importance[tmp_sigPC$PC[1]] * 100,
                                        digits = 2),
                                  "%), p_val: ",
                                  format(tmp_factor_anova[[tmp_sigPC$PC[1], tmp_factor]],
                                         digits = 2)
                    ),
                    ylab = paste0("PC",
                                  tmp_sigPC$PC[2],
                                  "(",
                                  round(tmp_importance[tmp_sigPC$PC[2]] * 100,
                                        digits = 2),
                                  "%), p_val: ",
                                  format(tmp_factor_anova[[tmp_sigPC$PC[2], tmp_factor]],
                                         digits = 2)
                    ),
                    title = paste0(tmp_factor,":"),
                    color = tmp_factor,
                    margin.params = list(fill = tmp_factor,
                                         color = "black", 
                                         size = 0.2),
                    margin.plot = "density",
                    palette = hcl.colors(
                      n = n_distinct(tmp_plottingData_PCA[, tmp_factor]),
                      "Dark 2"),
                    repel = TRUE,
                    print = FALSE) 
    
    tmp_plot <- ggpar(tmp_plot, tickslab = FALSE, legend = "none")
    
  
    tmp_caption <- "";
    if (length(tmp)) {
      tmp_caption <- paste0(tmp_factor," is significantly associated with PCs ",
                            paste(tmp, collapse = ", "), 
                            ".\n These corresponding to ",
                            sum(tmp_importance[tmp])*100,
                            "% of total explained variance.\n")
      
      if (sum(tmp_importance[tmp] > 0.01) > 0) {
        tmp_plotList2[[paste0(tmp_factor)]] <-
          grouped_ggbetweenstats(
            data = tmp_plottingData_PCA_long %>%
              filter(PC %in% tmp, var > 0.01) %>%
              mutate(p_value = tmp_factor_anova[as.numeric(PC), tmp_factor][[1]],
                     x = str_c(x, " p=", format(p_value, digits = 2)),
                     x = factor(x)),
            x = {{tmp_factor}},
            y = Score,
            grouping.var = x,
            plotgrid.args = list(ncol = if_else(length(tmp) > 1, 2, 1)),
            annotation.args = list(
              title = paste("PCs associated with", tmp_factor),
              subtitle = "(One way Anova p-value < 0.05, explained var > 0.01)"
            ),
            outlier.tagging = TRUE,
            outlier.label = joinKey,
            outlier.label.args = list(color = "red"),
            outlier.coef = 2.2,
            pairwise.comparisons = T,
            type = "r",
            pairwise.display = "s",
            p.adjust.method = "fdr",
            results.subtitle = T,
            centrality.plotting = T,
            plot.type = "violin",
            ggplot.component = list(
              scale_color_manual(values = hcl.colors(
                n = n_distinct(tmp_plottingData_PCA[, tmp_factor]),
                "Dark 2")),
              scale_fill_manual(values = hcl.colors(
                n = n_distinct(tmp_plottingData_PCA[, tmp_factor]),
                "Dark 2")),
              theme(axis.text.x = element_text(angle = 45),
                    axis.text.y = element_blank())
            ),
            ggtheme = theme_classic())
      }#end if
    } #end if
    
    tmp_layout <- c(area(t = 2, l = 1, b = 5, r = 4),
                    area(t = 1, l = 1, b = 1, r = 4),
                    area(t = 2, l = 5, b = 5, r = 5),
                    area(t = 1, l = 5, b = 1, r = 5))
    tmp_plot <- wrap_elements(wrap_plots(tmp_plot$sp,
                                         tmp_plot$xplot,
                                         tmp_plot$yplot,
                                         plot_spacer()) +
                                plot_layout(design = tmp_layout) +
                                plot_annotation(#title = paste0("Top 2 PCs associated with ", tmp_factor),
                                  #subtitle = "(by smallest one way Anova p-value)",
                                  caption = tmp_caption) &
                                theme(legend.position = "none") &
                                guides(fill = "none"))
    
    tmp_plotList[[tmp_factor]] <- tmp_plot;
  } #end for
  
  
  tmpHM <- tmp_factor_anova %>%
    column_to_rownames(var = "PC") %>%
    mutate(across(everything(), function(x) {
      x <- -log10(x)
      x[x>5] <- 5
      return(x)
    })) 
  
 tmpHM_matrix <- as.matrix(tmpHM) %>% t()  
  tmpHM <- Heatmap(tmpHM_matrix,
                   row_names_side = "left",
                   cluster_rows = F,
                   cluster_columns = F,
                   heatmap_legend_param = list(title = "-log10(p-val)",
                                               title_position = "leftcenter-rot",
                                               at = c(0, 1.3, 2.6, 3.9, 5)),
                   rect_gp = gpar(col = "white", lwd = 2),
                   cell_fun = function(j, i, x, y, width, height, fill) 
                   {
                     if(tmpHM_matrix[i,j] >= -log10(0.05)) {
                       grid.text("*", x, y, gp = gpar(fontsize = 20, col = "red"))
                     }
                     else grid.text("")
                   },
                   col = circlize::colorRamp2(c(0, 1.3,
                                                2.225, 3.150, 4.075, 5.000), c("black", "white",
                                                                               "#00B785","#53CC67","#B2DC3C","#FDE333")),
                   top_annotation = HeatmapAnnotation("%Variance" = anno_barplot(round(100*tmp_importance[1:nPC], digits = 1),
                                                                                 add_numbers = T,
                                                                                 numbers_rot = 0,
                                                                                 border = F)),
                   right_annotation = rowAnnotation(test = 
                                                      anno_text(paste0(round(digits = 1,
                                                                             100*colSums(tmp_importance[1:nPC]*(tmpHM>-log10(pval)), na.rm=TRUE)),
                                                                       "%"))))
  
  #Make sure we include all significant factors, or the core ones at least
  tmp_sigFactors <- (tmp_factor_anova %>%
    column_to_rownames(var = "PC") %>%
    mutate(across(everything(), function(x) x<pval)) %>%
      colSums() > 0)
  tmp_sigFactors <- unique(c(names(tmp_sigFactors)[tmp_sigFactors], factorsofInterest))
  
  
    tmp_plot <- wrap_plots(tmp_plotList[tmp_sigFactors], 
                         ncol = if_else(length(tmp_sigFactors) > 1, 2, 1)) +
    plot_annotation(title = paste0("Top 2 PCs associated with factors of interest"),
                    subtitle = "(by heteroskedastic Anova p-value)")
  
  #rename the secondary plots
  if (plot) {
    draw(tmpHM)
    print(tmp_plot)

    
    
     for (tmp_figure in names(tmp_plotList)) {
       if(!is.null(tmp_plotList2[[tmp_figure]]))
             print(wrap_plots(
               tmp_plotList[[tmp_figure]],
               tmp_plotList2[[tmp_figure]],
               nrow = 1))
       else {print(tmp_plotList[[tmp_figure]])}
    }#end for 
  } #end if
  
  tmp_factor_anova$var = tmp_importance[1:nPC]
    
  return(list(plottingData = tmp_plottingData_PCA,
              plottingData_long = tmp_plottingData_PCA_long,
              PCA_data =  tmp_pcaRes,
              heatmap = tmpHM,
              summary_plot = tmp_plot,
              all_plots = tmp_plotList,
              factor_plots = tmp_plotList2,
              anova = tmp_factor_anova))
} #plot_PCA

#Quiet plot_PCA
quiet_plot_PCA <- function(geneCounts,
                           sampleFactors, 
                           factorsofInterest,
                           nPC = 0,
                           fdr = "none",
                           pval = 0.05) {
  stopifnot(require(plyr))
  stopifnot(require(tidyverse))
  
  return( quietly(plot_PCA)(geneCounts,
                            sampleFactors,
                            factorsofInterest,
                            nPC,
                            plot = FALSE,
                            fdr = fdr,
                            pval = pval)$result)
}#end quiet_plot_PCA
