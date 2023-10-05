source("summary_funcs.r")
#read in results
results <- read.csv("results/all_results.csv", row.names=1)
#produce tables as before for just resNMTF
results_nmtf <- subset(results, Method == "ResNMTF")
colnames(results_nmtf) <- c("Method", "Factor", "Choice", "Measure", "Mean", "S.d.")
#make and saves plot
number <- length(row_cl_dims)
n_repeats <- 10
subhead <- TeX(paste0(4, " views, ",
                 "$\\phi$", " = ", 200, ", ",
                  n_repeats, " repetitions."))
x_title <- "Number of biclusters"
plot_title <- "The effect of increasing the number of biclusters"
path_to_save <-  "results/tester_increasing_bicl_plot.pdf"
measures <- c("CSR", "Relevance", "Recovery", "F score", "BiS")
measures <- c("BiS", "CSR", "F score", "Recovery", "Relevance")
suppressMessages(make_plot_line(results_nmtf, "Overall", plot_title,
                subhead, x_title, y_title,
                path_to_save, measures))

colnames(results) <- c("Method", "Factor", "Choice", "Measure", "Mean", "S.d.")
suppressMessages(make_plot_measure(results, "Overall", plot_title,
                subhead, x_title, "Measure",
                path_to_save, "Rel", c("GFA", "iSSVD", "NMTF", "ResNMTF")))
p_rec <- suppressMessages(make_plot_measure(results, "Overall", plot_title,
                subhead, x_title, "Recovery",
                path_to_save, "Rec", c("GFA", "iSSVD","NMTF", "ResNMTF")))

p_rel <- suppressMessages(make_plot_measure(results, "Overall", plot_title,
                subhead, x_title, "Relevance",
                path_to_save, "Rel", c("GFA", "iSSVD", "NMTF", "ResNMTF")))     
     
library("ggpubr")

measures <- c("Relevance", "Recovery", "F score", "BiS", "CSR")
n_measures <- length(measures)
measure_vec <-  c("Rel", "Rec", "F_score", "BiS", "CSR")
plots <- vector("list", length = n_measures)
for(i in 1:n_measures){
  path_to_save <-   paste0("results/tester_",measure_vec[i], "_plot.pdf")
  plots[[i]] <- make_plot_measure(results, "Overall", plot_title,
                subhead, x_title, measures[i],
                path_to_save, measure_vec[i], c("GFA", "iSSVD","NMTF", "ResNMTF"))
}
  ggarrange(
    plots[[1]], plots[[2]],plots[[3]], plots[[4]], plots[[5]],  labels = c("A", "B", "C", "D", "E"), ncol=2, nrow=3,
    common.legend = TRUE, legend = "bottom"
    )

  ggarrange(
    plots, labels = c("A", "B", "C", "D", "E"), ncol=2, nrow=3,
    common.legend = TRUE, legend = "bottom"
    )
plot(1:10)
TeX(paste0(4, " views, ",
                 "$\\phi$", " = ", 200, ", ",
                  n_repeats, " repetitions."), output="expression")
paste0(4, " views, ",
                  "$\\phi$", " = ", 200, ", ",
                 n_repeats, " repetitions.")
#produce correlation table
corr_table <- function(results, subheading, filename){
  measures <- c("F_score", "Rel", "Rec", "BiS", "CSR")
  corr_tab <- matrix(0, nrow=3, ncol=2)
  for(i in 1:3){
    for(j in 1:2){
       vec1 <- results$Mean[results$Measure == measures[i]]
      vec2 <- results$Mean[results$Measure == measures[3 + j]]
      corr_tab[i, j] <- cor(vec1, vec2)
    }
  }
  rownames(corr_tab) <- c("F score", "Relevance", "Recovery")
  colnames(corr_tab) <- c("BiS", "CSR")
  latex_tab <- kbl(corr_tab, booktabs=T, "latex", escape=FALSE)
  latex_tab <- kable_styling(add_footnote(latex_tab, subheading))
   #save
  sink(paste0(filename, "corr_tab.txt"))
  print(latex_tab)
  sink()
  return(latex_tab)
}


cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
cor_res <-corr_table(results, "hello", "results/")



df_overall <- kable_styling(add_footnote(df_overall2, subheading))

col_names <- c("Mes", "CSR", "Relevance", "Recovery", "F score", "BiS") 
path_to_save <- "results/resNMTF_table_test"
measure_labels <- c("CSR", "Rel", "Rec", "F_score", "BiS")
measure_labels <- c("F_score",  "Rel", "Rec", "CSR" ,  "BiS")

col_names <- c("", "F score", "Relevance", "Recovery", "CSR", "BiS") 
all_table2(results_nmtf, measure_labels, path_to_save, 10, col_names, subhead, feat = "Factor")



method_labels <- c("ResNMTF", "NMTF", "GFA", "iSSVD")
col_names <- c("", method_labels)
n_methods <- length(method_labels)
for(i in 1:1){
  path_to_save <-   paste0("results/tester_",measure_vec[i], "_table")
  results_measure <- subset(results, Measure == measure_vec[i])
  all_table2(results_measure, method_labels, path_to_save, 10, col_names, subhead, feat = "Factor")
}
feat <- "Factor"
ndig <-4



rel <- measure_table(results,"Rel", "results/test", 10, 2:5, "hello", "subscript of table", ndig = 4)

method_labels <- c("ResNMTF", "NMTF", "GFA", "iSSVD")
col_names <- c("", method_labels)
n_methods <- length(method_labels)
tables <- vector("list", length = length(measure_vec))
for(i in 1:length(measure_vec)){
  path_to_save <-   paste0("results/tester_",measure_vec[i], "_table")
  tables[[i]] <- measure_table(results,measure_vec[i], 
    path_to_save, 10, 2:5, "hello", "subscript of table", ndig = 4)
}
measures <- c("Relevance", "Recovery", "F score", "BiS", "CSR")
rbind.data.frame(tables[[1]], tables[[2]], tables[[3]], tables[[4]])
combine_tables(tables, "results/test", measures, "hello", "hi", as.character(2:5))

vec1 <- c(" " = 1, "table_head" = 4)
vec1[["a"]] <- 2
name <- "hi"
vec1[[name]] <- 5
vec <- c(" " = 1)
vec[["table_head"]] <- 5
