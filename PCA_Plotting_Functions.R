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
                     plot = TRUE) {
  #Tests
  stopifnot(require(plyr))
  stopifnot(require(tidyverse))
  stopifnot(require(ggpubr))
  stopifnot(require(ggstatsplot))
  stopifnot(is.data.frame(geneCounts))
  stopifnot(is.data.frame(sampleFactors))
  stopifnot(ncol(geneCounts) == nrow(sampleFactors))
  stopifnot(all(colnames(geneCounts) == rownames(sampleFactors)))
  stopifnot(length(factorsofInterest)>0)
  stopifnot(all(factorsofInterest %in% colnames(sampleFactors)))

  #Make sure we have enough samples in each 
  # keepFactors <- c()
  # for (tmp_factor in factorsofInterest) {
  #   if ((sampleFactors %>% pull(!!sym(tmp_factor)) %>% levels %>% length > 1) &
  #       (sum(sampleFactors %>% 
  #           count(!!sym(tmp_factor)) %>%      
  #           pull(n) < 2) == 0)) {
  #     keepFactors <- c(keepFactors, tmp_factor)
  #   } else {
  #     warning(paste0("plot_PCA: Too few samples in factor of interest ",
  #                    tmp_factor, ".  Excluding from further analysis.\n"),
  # #             call. = FALSE)
  # #   }
  # # } #end for
  # factorsofInterest <- keepFactors
  # stopifnot(length(factorsofInterest)>0)
  
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
                mutate(across(all_of(factorsofInterest), as_factor)) %>%
                mutate(joinKey = rownames(.)) %>%
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
  
  #Use a one-way Anova to see if there are any differences based on factors of
  #interest, and estimate the total explained variance associated with them
  tmp_factor_anova = list()
  for (tmp_factor in factorsofInterest) {
    tmp <- tmp_plottingData_PCA_long %>%
      group_by(PC) %>%
      group_modify(~statsExpressions::oneway_anova(
        data = .x,
        x = {{tmp_factor}},
        y = Score), .keep = TRUE) %>%
      arrange(PC)
    
    #Store the p.values
    tmp_factor_anova[[tmp_factor]] <- tmp$p.value
  } #end for
  #convert to data frame, add PCs
  tmp_factor_anova <- tmp_factor_anova %>% 
    data.frame %>%
    mutate(PC = row_number(),
           var = tmp_importance[as.numeric(PC)])
  
  # We want to associate each PC with a single factor
  tmp_factor_anova <- tmp_factor_anova %>%
    rowwise %>%
    mutate(prinFactor = which.min(c_across({{factorsofInterest}})),
           prinFactor = {{factorsofInterest}}[prinFactor]) %>%
    ungroup()
  
  tmp_plotList <- list()
  for (tmp_factor in factorsofInterest) {
    #We want this function to be as informative as possible, so we will
    #allow for multiple factors being associated with 
    tmp <- which(tmp_factor_anova[, tmp_factor] < 0.05)

    tmp_sigPC <- tmp_factor_anova %>%
      slice_min(order_by = !!sym(tmp_factor), n = 2)
    
    tmp_plot <- 
      ggscatterhist(tmp_plottingData_PCA,
                    x = paste0("PC", tmp_sigPC$PC[1]),
                    y = paste0("PC", tmp_sigPC$PC[2]),
                    label = tmp_factor,
                    font.label = c(8, "plain"),
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
                    show.legend = FALSE,
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
        tmp_plotList[[paste0(tmp_factor,"-2")]] <-
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
            pairwise.comparisons = FALSE,
            pairwise.display = "s",
            p.adjust.method = "fdr",
            results.subtitle = FALSE,
            centrality.plotting = FALSE,
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
                                         guide_area()) +
                                plot_layout(design = tmp_layout,
                                            guides = "collect") +
                                plot_annotation(#title = paste0("Top 2 PCs associated with ", tmp_factor),
                                  #subtitle = "(by smallest one way Anova p-value)",
                                  caption = tmp_caption) &
                                theme(legend.position = "none") &
                                guides(fill = "none"))
    
    tmp_plotList[[tmp_factor]] <- tmp_plot;
  } #end for
  
  tmp_plot <- wrap_plots(tmp_plotList[factorsofInterest], 
                         ncol = if_else(length(factorsofInterest) > 1, 2, 1)) +
    plot_annotation(title = paste0("Top 2 PCs associated with factors of interest"),
                    subtitle = "(by one way Anova p-value)")
  
  #rename the secondary plots
  tmp_plotList2 <- list()
  for (tmp_factor in factorsofInterest) {
    if (!is.null(paste0(tmp_factor,"-2")))
      tmp_plotList2[[tmp_factor]] <- tmp_plotList[[paste0(tmp_factor,"-2")]]
  }
  if (plot) {
    print(tmp_plot)

     for (tmp_factor in factorsofInterest) {
      if (!is.null(paste0(tmp_factor,"-2")))
        print(tmp_plotList[[paste0(tmp_factor,"-2")]])
    }#end for 
  } #end if
  
  return(list(plottingData = tmp_plottingData_PCA,
              plottingData_long = tmp_plottingData_PCA_long,
              PCA_data =  tmp_pcaRes,
              summary_plot = tmp_plot,
              factor_plots = tmp_plotList2,
              anova = tmp_factor_anova))
} #plot_PCA
