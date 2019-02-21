# Analysis of CEDAR's Value Recommender results

library(ggplot2)
library(gridExtra)
library(ggpubr)
library("RColorBrewer")

### CONSTANTS ###
workspace = "../../../data/results/main_experiment"
color1 <- "#fdb863" # light orange
color2 <- "#e66101" # dark orange
color3 <- "#b2abd2" # light violet
color4 <- "#5e3c99" # dark violet
  
# Results data (free text)
file_NCBItoNCBI_text <- paste(workspace, "/text/results_trainNCBI_testNCBI_text.csv", sep="")
file_NCBItoEBI_text <- paste(workspace, "/text/results_trainNCBI_testEBI_text.csv", sep="")
file_EBItoEBI_text <- paste(workspace, "/text/results_trainEBI_testEBI_text.csv", sep="")
file_EBItoNCBI_text <- paste(workspace, "/text/results_trainEBI_testNCBI_text.csv", sep="")

# Results data (annotated, using the same ontologies for NCBI and EBI datasets and using mappings)
file_NCBItoNCBI_annotated <- paste(workspace, "/annotated/results_trainNCBI_testNCBI_annotated.csv", sep="")
file_EBItoEBI_annotated <- paste(workspace, "/annotated/results_trainEBI_testEBI_annotated.csv", sep="")
file_NCBItoEBI_annotated <- paste(workspace, "/annotated/results_trainNCBI_testEBI_annotated.csv", sep="")
file_EBItoNCBI_annotated <- paste(workspace, "/annotated/results_trainEBI_testNCBI_annotated.csv", sep="")

### FUNCTION DEFINITIONS ###

standardFieldName <- function(field_name) {
  if (field_name == 'cell_line' | field_name == 'cellLine') {
    return('cell line')
  }
  else if (field_name == 'cell_type' | field_name == 'cellType') {
    return('cell type')
  }
  else if (field_name == 'diseaseState') {
    return('disease')
  }
  else if (field_name == 'organismPart') {
    return('tissue')
  }
  else {
    return(field_name)
  }
}

# Aggregation by no_populated_fields
aggregate_data_1 <- function(df, reciprocal_rank_vr_column, reciprocal_rank_baseline_column) {
  # aggregation for the 'recommender' method
  agg1 <- aggregate(list(mrr=df[[reciprocal_rank_vr_column]]), by=list(no_populated_fields = df$populated_fields_size), FUN=mean)
  agg1$method <- "recommender"
  # aggregation for the 'baseline' method
  agg2 <- aggregate(list(mrr=df[[reciprocal_rank_baseline_column]]), by=list(no_populated_fields = df$populated_fields_size), FUN=mean)
  agg2$method <- "baseline"
  # final aggregation
  agg_final <- rbind(agg1, agg2)
  # Limit it to no_populated_fields <5
  agg_final <- agg_final[agg_final$no_populated_fields < 5,]
  agg_final$experiment <- df$experiment[1]
  return(agg_final)
}

# Aggregation by no_populated_fields
aggregate_data_1_2 <- function(df, reciprocal_rank_vr_column, reciprocal_rank_baseline_column, recom_method_name="recommender", baseline_method_name="baseline") {
  # aggregation for the 'recommender' method
  agg1 <- aggregate(list(mrr=df[[reciprocal_rank_vr_column]]), by=list(no_populated_fields = df$populated_fields_size), FUN=mean)
  agg1$method <- recom_method_name
  if (!is.null(reciprocal_rank_baseline_column)) {
    # aggregation for the 'baseline' method
    agg2 <- aggregate(list(mrr=df[[reciprocal_rank_baseline_column]]), by=list(no_populated_fields = df$populated_fields_size), FUN=mean)
    agg2$method <- baseline_method_name
    # final aggregation
    agg_final <- rbind(agg1, agg2)
  }
  else {
    agg_final <- agg1
  }
  # Limit it to no_populated_fields <5
  agg_final <- agg_final[agg_final$no_populated_fields < 5,]
  agg_final$experiment <- df$experiment[1]
  return(agg_final)
}

# Aggregation by target_field and no_populated_fields
aggregate_data_2 <- function(df, reciprocal_rank_vr_column, reciprocal_rank_baseline_column, recom_method_name="recommender", baseline_method_name="baseline") {
  # aggregation for the 'recommender' method
  agg1 <- aggregate(list(mrr=df[[reciprocal_rank_vr_column]]), by=list(field = df$target_field, no_populated_fields = df$populated_fields_size), FUN=mean)
  agg1$method <- recom_method_name
  # aggregation for the 'baseline' method
  agg2 <- aggregate(list(mrr=df[[reciprocal_rank_baseline_column]]), by=list(field = df$target_field, no_populated_fields = df$populated_fields_size), FUN=mean)
  agg2$method <- baseline_method_name
  # final aggregation
  agg_final <- rbind(agg1, agg2)
  # convert from factors to characters to be able to modify them
  agg_final$field <- as.character(agg_final$field)
  
  # Limit it to no_populated_fields <5
  agg_final <- agg_final[agg_final$no_populated_fields < 5,]
  agg_final$experiment <- df$experiment[1]
  
  # Generate the right field name
  for(i in 1:length(agg_final$field)){
    s <- agg_final$field[i]
    if (startsWith(s, '[')) {
      separator <- '\\]\\('
      pos_separator <- regexpr(pattern = separator, s)
      var_name <- substr(s, pos_separator + 2, nchar(s) - 1)
    }
    else {
      var_name <- s
    }
    var_name <- standardFieldName(var_name)
    agg_final$field[i] <- var_name
  }
  
  agg_final$field <- as.factor(agg_final$field)
    
  return(agg_final)
}

generate_plot_2 <- function(df, title="title") {
  # Custom order for the factors
  df$method <- factor(df$method, c("baseline (text)", "baseline (ontology terms)","recommender (text)", "recommender (ontology terms)"))
  
  plot <- ggplot(data=df, aes(x=no_populated_fields, y=mrr, group=method, colour=method)) +
    theme_bw(base_size = 9)  +
    geom_line(aes(linetype=method), size=0.5) +
    scale_linetype_manual(values=c("dotted", "dotted", "solid", "solid")) +
    scale_color_manual(values=c(color2, color4, color2, color4)) +
    geom_point() +
    ylim(0,1) + ggtitle(title) + xlab("No. populated fields") + ylab("Mean Reciprocal Rank")
  return(plot)
}

generate_plot_field <- function(df, title="title") {
  # Custom order for the factors
  df$method <- factor(df$method, c("baseline (text)", "baseline (ontology terms)","recommender (text)", "recommender (ontology terms)"))
  plot <- ggplot(data=df, aes(x=field, y=mrr, fill=method)) + geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c(color1, color3, color2, color4)) +
    ylim(0,1) + ggtitle(title) + xlab("Field") + ylab("Mean Reciprocal Rank") +
    theme_bw() + theme(text = element_text(size=8))
  return(plot)
}

setClass("EvaluationSet", representation(datasets = "vector", description = "character"))

generate_all_plots <- function(evaluation_set, reciprocal_rank_vr_column, reciprocal_rank_baseline_column) {
  
  data_NCBItoNCBI <- read.csv(evaluation_set@datasets[1])
  data_NCBItoEBI <- read.csv(evaluation_set@datasets[2])
  data_EBItoEBI <- read.csv(evaluation_set@datasets[3])
  data_EBItoNCBI <- read.csv(evaluation_set@datasets[4])
  description <- evaluation_set@description
  
  # remove the 'treatment' field from the analysis
  data_NCBItoNCBI <- data_NCBItoNCBI[data_NCBItoNCBI$target_field!="treatment",]
  
  data_NCBItoNCBI$experiment <- "NCBItoNCBI"
  data_NCBItoEBI$experiment <- "NCBItoEBI"
  data_EBItoEBI$experiment <- "EBItoEBI"
  data_EBItoNCBI$experiment <- "EBItoNCBI"

  # 1) Recommender vs Baseline 2x2 plots
  p1 <- generate_plot(aggregate_data_1(data_NCBItoNCBI, reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: NCBI; Testing: NCBI")
  p2 <- generate_plot(aggregate_data_1(data_NCBItoEBI, reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: NCBI; Testing: EBI")
  p3 <- generate_plot(aggregate_data_1(data_EBItoEBI, reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: EBI; Testing: EBI")
  p4 <- generate_plot(aggregate_data_1(data_EBItoNCBI, reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: EBI; Testing: NCBI")
  fig1 <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  desc_text <- paste("Metadata Recommender vs Baseline (", description, ")", sep = "")
  fig1_annotated <- annotate_figure(fig1, top = text_grob(label=desc_text, color = "black", face = "bold", size = 11))
  print(fig1_annotated)
  
  # Export plot
  dev.copy(pdf, paste("plot1_", gsub(" ", "_", description), "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".pdf", sep=""))
  dev.off()
  
  # 2) Recommender vs Baseline per target field
  p1 <- generate_plot_field(aggregate_data_2(data_NCBItoNCBI,reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: NCBI; Testing: NCBI")
  p2 <- generate_plot_field(aggregate_data_2(data_NCBItoEBI,reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: NCBI; Testing: EBI")
  p3 <- generate_plot_field(aggregate_data_2(data_EBItoEBI,reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: EBI; Testing: EBI")
  p4 <- generate_plot_field(aggregate_data_2(data_EBItoNCBI,reciprocal_rank_vr_column, reciprocal_rank_baseline_column), "Training: EBI; Testing: NCBI")
  fig2 <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  desc_text <- paste("Metadata Recommender vs Baseline by field (", description, ")", sep = "")
  fig2_annotated <- annotate_figure(fig2, top = text_grob(label=desc_text, color = "black", face = "bold", size = 11))
  print(fig2_annotated)
  
  # Export plot
  dev.copy(pdf, paste("plot2_", gsub(" ", "_", description), "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".pdf", sep=""))
  dev.off()
}

generate_all_plots_overlapped <- function(evaluation_set1, evaluation_set2, reciprocal_rank_vr_column, reciprocal_rank_baseline_column) {
  
  data_NCBItoNCBI1 <- read.csv(evaluation_set1@datasets[1])
  data_NCBItoEBI1 <- read.csv(evaluation_set1@datasets[2])
  data_EBItoEBI1 <- read.csv(evaluation_set1@datasets[3])
  data_EBItoNCBI1 <- read.csv(evaluation_set1@datasets[4])
  description1 <- evaluation_set1@description
  
  data_NCBItoNCBI2 <- read.csv(evaluation_set2@datasets[1])
  data_NCBItoEBI2 <- read.csv(evaluation_set2@datasets[2])
  data_EBItoEBI2 <- read.csv(evaluation_set2@datasets[3])
  data_EBItoNCBI2 <- read.csv(evaluation_set2@datasets[4])
  description2 <- evaluation_set2@description
  
  data_NCBItoNCBI1$experiment <- "NCBItoNCBI"
  data_NCBItoEBI1$experiment <- "NCBItoEBI"
  data_EBItoEBI1$experiment <- "EBItoEBI"
  data_EBItoNCBI1$experiment <- "EBItoNCBI"
  
  data_NCBItoNCBI2$experiment <- "NCBItoNCBI"
  data_NCBItoEBI2$experiment <- "NCBItoEBI"
  data_EBItoEBI2$experiment <- "EBItoEBI"
  data_EBItoNCBI2$experiment <- "EBItoNCBI"
  
  # Used to include or exclude the baseline from the evaluation
  baseline_col = reciprocal_rank_baseline_column
  #baseline_col = NULL
  
  # 1) Recommender vs Baseline 2x2 plots
  
  data_p1_1 <- aggregate_data_1_2(data_NCBItoNCBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p1_2 <- aggregate_data_1_2(data_NCBItoNCBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p1 <- rbind(data_p1_1, data_p1_2)
  
  data_p2_1 <- aggregate_data_1_2(data_NCBItoEBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p2_2 <- aggregate_data_1_2(data_NCBItoEBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p2 <- rbind(data_p2_1, data_p2_2)
  
  data_p3_1 <- aggregate_data_1_2(data_EBItoEBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p3_2 <- aggregate_data_1_2(data_EBItoEBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p3 <- rbind(data_p3_1, data_p3_2)
  
  data_p4_1 <- aggregate_data_1_2(data_EBItoNCBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p4_2 <- aggregate_data_1_2(data_EBItoNCBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p4 <- rbind(data_p4_1, data_p4_2)

  p1 <- generate_plot_2(data_p1, "Training: NCBI; Testing: NCBI")
  p2 <- generate_plot_2(data_p2, "Training: NCBI; Testing: EBI")
  p3 <- generate_plot_2(data_p3, "Training: EBI; Testing: EBI")
  p4 <- generate_plot_2(data_p4, "Training: EBI; Testing: NCBI")
  fig1 <- ggarrange(p1, p3, p2, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  description = paste(description1, " vs ", description2, sep = "")
  #desc_text <- paste("Results", sep = "")
  desc_text <- ""
  fig1_annotated <- annotate_figure(fig1, top = text_grob(label=desc_text, color = "black", face = "bold", size = 11))
  print(fig1_annotated)

  # Export plot
  dev.copy(pdf, paste("plot_MRR", "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".pdf", sep=""))
  dev.off()
  
  # 2) Recommender vs Baseline (per target field)
  data_p5_1 <- aggregate_data_2(data_NCBItoNCBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p5_2 <- aggregate_data_2(data_NCBItoNCBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p5 <- rbind(data_p5_1, data_p5_2)
  
  data_p6_1 <- aggregate_data_2(data_NCBItoEBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p6_2 <- aggregate_data_2(data_NCBItoEBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p6 <- rbind(data_p6_1, data_p6_2)
  
  data_p7_1 <- aggregate_data_2(data_EBItoEBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p7_2 <- aggregate_data_2(data_EBItoEBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p7 <- rbind(data_p7_1, data_p7_2)
  
  data_p8_1 <- aggregate_data_2(data_EBItoNCBI1, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (text)", baseline_method_name="baseline (text)")
  data_p8_2 <- aggregate_data_2(data_EBItoNCBI2, reciprocal_rank_vr_column, baseline_col, recom_method_name="recommender (ontology terms)", baseline_method_name="baseline (ontology terms)")
  data_p8 <- rbind(data_p8_1, data_p8_2)
  
  p5 <- generate_plot_field(data_p5, "Training: NCBI; Testing: NCBI")
  p6 <- generate_plot_field(data_p6, "Training: NCBI; Testing: EBI")
  p7 <- generate_plot_field(data_p7, "Training: EBI; Testing: EBI")
  p8 <- generate_plot_field(data_p8, "Training: EBI; Testing: NCBI")
  
  fig2 <- ggarrange(p5, p7, p6, p8, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  desc_text <- ""
  fig2_annotated <- annotate_figure(fig2, top = text_grob(label=desc_text, color = "black", face = "bold", size = 11))
  print(fig2_annotated)
  
  
  # Export plot
  dev.copy(pdf, paste("plot_MRR_per_field", "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".pdf", sep=""))
  dev.off()
  
}

### MAIN BODY ###
setwd(workspace)
evaluation_set_1 <- new("EvaluationSet", datasets=c(file_NCBItoNCBI_text, file_NCBItoEBI_text, file_EBItoEBI_text, file_EBItoNCBI_text), description="free text")
evaluation_set_2 <- new("EvaluationSet", datasets=c(file_NCBItoNCBI_annotated, file_NCBItoEBI_annotated, file_EBItoEBI_annotated, file_EBItoNCBI_annotated), description="annotated; diff ontologies; with mappings")
generate_all_plots_overlapped(evaluation_set_1, evaluation_set_2, 'RR_top5_vr', 'RR_top5_baseline')

