suppressWarnings(suppressPackageStartupMessages(library("rio")))
library(ggplot2)
library("tidyr")
library("dplyr")
#import as list of matrices instead of as a list
import_matrix <- function(filename){
    return(lapply(import_list(filename), function(x) as.matrix(x)))
}
calc_averages <- function(mat_res){return(c(rowMeans(mat_res), mean(mat_res)))}

get_averages <- function(list_avgs, data_name, filename, repeat_no, measure){
    results <- import_matrix(paste0(paste0(data_name,"/"), paste0(filename,".xlsx")))
    for (j in 1:length(measure)){
            #list_avgs[[j]][repeat_no, ] <- calc_averages(results[[j]])
            list_avgs[[j]][repeat_no, ] <- mean(results[[measure[j]]])      
    }
    return(list_avgs)
    }

avg_over_repeats <- function(avgs, i, measure, name, cluster){
    rep_avgs <- colMeans(avgs[[i]])
    rep_sd <- apply(avgs[[i]], 2, sd)
    return(cbind(name, cluster, measure, rep_avgs, rep_sd))
}

avg_over_repeats2 <- function(avgs, i, measure, name1, name2, cluster){
    rep_avgs <- colMeans(avgs[[i]])
    rep_sd <- apply(avgs[[i]], 2, sd)
    return(cbind(name1, name2, cluster, measure, rep_avgs, rep_sd))
}

make_plot <- function(data_frame, cluster, plot_title, sub_title, filename){
    p <- ggplot(subset(data_frame, Choice==cluster), aes(x=Measure, y=Mean, fill=Method)) + 
    geom_bar(stat="identity", position=position_dodge()) +
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(.9))+
     expand_limits(y=c(0,1))+
     scale_y_continuous(breaks = seq(0, 1, length=6))+
     labs(
     title = plot_title,
     subtitle = sub_title
   )
    p<-p+  scale_fill_brewer() + theme(text = element_text(size = 16),legend.position = "bottom")    
    ggsave(filename,plot=p, compress=FALSE,device="pdf")
}

make_plot_params <- function(data_frame, measure, plot_title, sub_title, filename){
    p <- ggplot(subset(data_frame, Measure==measure), aes(x=Method, y=Mean, fill=Choice)) + 
    geom_bar(stat="identity", position=position_dodge()) +
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(.9))+
     expand_limits(y=c(0,1))+
     scale_y_continuous(breaks = seq(0, 1, length=6))+
     labs(
     title = plot_title,
     subtitle = sub_title
   )
    p<-p+  scale_fill_brewer() + theme(text = element_text(size = 16),legend.position = c(0.8,0.2))    
    ggsave(filename,plot=p, compress=FALSE,device="pdf")
}
#position_dodge(0.05)
make_plot_line <- function(data_frame, cluster, plot_title, sub_title, x_title, y_title, filename, measures){
    if(min(data_frame$Mean)<0){
        y_min <- min(data_frame$Mean)
        y_lims <- c(y_min,1)
        if(y_min< -0.2){
            y_breaks <- c(sort(seq(-0.2, y_min, by = -0.2)),
                     seq(0, 1, length = 6))
        }else{
            y_breaks <- c(-0.2, seq(0, 1, length = 6) )
        } 
    }else{
        y_lims <- c(0,1)
        y_breaks <- seq(0, 1, length = 6)
    }
    p <- ggplot(subset(data_frame, Choice==cluster), aes(x=Factor, y=Mean, group=Measure, color=Measure)) + 
    geom_line(stat="identity") +
    geom_point() +
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(0.01))+
     expand_limits(y=y_lims)+
     scale_color_hue(labels = measures) +
     scale_y_continuous(breaks = y_breaks)+
     labs(
     title = plot_title,
     subtitle = sub_title,
     x = x_title, 
     y = y_title
   )
    p <- p +
        theme(text = element_text(size = 11), legend.position = "bottom")
    suppressMessages(ggsave(filename, plot = p, compress = FALSE, device="pdf"))
    return(p)
}

make_plot_measure <- function(data_frame, cluster, plot_title, sub_title, x_title, y_title, filename, measure, methods){
    data_frame <- subset(data_frame, Measure==measure)
    if(min(data_frame$Mean)<0){
        y_min <- min(data_frame$Mean)
        y_lims <- c(y_min,1)
        if(y_min< -0.2){
            y_breaks <- c(sort(seq(-0.2, y_min, by = -0.2)),
                     seq(0, 1, length = 6))
        }else{
            y_breaks <- c(-0.2, seq(0, 1, length = 6) )
        } 
    }else{
        y_lims <- c(0,1)
        y_breaks <- seq(0, 1, length = 6)
    }
    p <- ggplot(subset(data_frame, Choice==cluster), aes(x=Factor, y=Mean, group=Method, color=Method)) + 
    geom_line(stat="identity") +
    geom_point() +
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(0.01))+
     expand_limits(y=y_lims)+
     scale_color_hue(labels = methods) +
     scale_y_continuous(breaks = y_breaks)
    p_save <- p + labs(x = x_title, y = y_title)
    p <- p + labs(
     title = plot_title,
     subtitle = sub_title,
     x = x_title, 
     y = y_title
   )
    p <- p +
        theme(text = element_text(size = 11), legend.position = "bottom")
    suppressMessages(ggsave(filename, plot = p, compress = FALSE, device="pdf"))
    return(list("a" = p_save, "b" = p))
}


make_plot_line_phi <- function(data_frame, cluster, plot_title, sub_title, x_title, filename){
    p <- ggplot(subset(data_frame, Choice==cluster), aes(x=Method, y=Mean, group=Measure, color=Measure)) + 
    geom_line(stat="identity") +
    geom_point() +
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(0.01))+
     expand_limits(y=c(0,1))+
     scale_color_hue(labels = c("CSR", "F score", "BiS", "Restrictions")) +
     scale_y_continuous(breaks = seq(0, 1, length = 6))+
     labs(
     title = plot_title,
     subtitle = sub_title,
     x = x_title
   )
    p <- p +
        theme(text = element_text(size = 11), legend.position = c(0.8, 0.2))
    suppressMessages(ggsave(filename, plot = p, compress = FALSE, device="pdf"))
}

all_plots <- function(plots, filename, n_col_plots, n_row_plots){
  all_plots <- ggarrange(
        plots[[1]]$a, plots[[2]]$a,plots[[3]]$a, plots[[5]]$a,
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 8),
           ncol = n_col_plots, nrow = n_row_plots,
        common.legend = TRUE, legend = "bottom")
  suppressMessages(ggsave(filename, plot = all_plots, compress = FALSE, device="pdf"))
  return(all_plots)
}

sd_brackets <- function(mean,sd,n, measure="other"){
    if(measure == "CSR"){
        return(paste0(paste0(n * as.numeric(mean), " ("),
                paste0(n * as.numeric(sd), ")")))
    }
    return(paste0(paste0(mean, " ("), paste0(sd,")")))
}
sd_brackets <- function(mean,sd,n, measure="other"){
    return(paste0(paste0(mean, " ("), paste0(sd,")")))
}
collapse_rows_df <- function(df, variable){
  group_var <- enquo(variable)
  df %>%
    group_by(!! group_var) %>%
    mutate(groupRow = 1:n()) %>%
    ungroup() %>%
    mutate(!!quo_name(group_var) := ifelse(groupRow == 1, as.character(!! group_var), "")) %>%
    select(-c(groupRow))
}

max_index <- function(col) {
    return(which.max(as.numeric(sub(" .*", "", as.vector(col)))))
}
min_index <- function(col) {
    return(which.min(as.numeric(sub(" .*", "", as.vector(col)))))
}
#updated all_table
all_table <- function(df, measures, filename, repeats, ndig = 4, feat = "Method"){

    df <- mutate_at(df, c(feat, "Cluster", "Measure"), as.factor)
    df <- pivot_wider(df, names_from = Measure, values_from = c("Mean","S.d."))
    #round to 4 digits
    df <- mutate_if(df, is.numeric,
            function(x) format(round(x, ndig), nsmall = ndig))
    #
    #add info in brackets
    mean_col_names <- paste0("Mean_",measures)
    sd_col_names <- paste0("S.d._",measures)
    #name of mean column
    for (i in 1:length(measures)){
        df[[measures[i]]] <- sd_brackets(df[[mean_col_names[i]]],
                                df[[sd_col_names[i]]], repeats, measures[i])
    }
    df_all <- subset(df, select = c("Cluster", measures))
    measures[["CSR"]] <- "CS"
    colnames(df_all) <- c(feat, measures)


    #collapse method rows
    group_tab <- kable_styling(kbl(df_all, booktabs = T, "latex"))
    #select groups
    feats <- unique(df[[feat]])
    for(i in 1:length(feats)){
        group_tab <- pack_rows(group_tab,feats[i], 3*i-2, 3*i)
    }
    #subset only overall results and remove extras
    df_overall <- subset(df,Cluster=="Overall",select=c(feat,measures))
    #make the largest values down a column bold
    for(i in 1:length(measures)){
        max_id <- max_index(df_overall[[measures[i]]])
        cond <- (1:length(df_overall[,1][[1]]))==max_id
        df_overall[[measures[i]]] <- text_spec(df_overall[[measures[i]]] , "latex", 
            bold = cond)
    }

    df_overall <- kbl(df_overall,booktabs=T,"latex",escape=FALSE)
    df_overall <- kable_styling(add_footnote(df_overall, c("Best performance in bold.")))
    #save
    sink(paste0(filename, "_overall.txt"))
    print(df_overall)
    sink()
    sink(paste0(filename, "_all.txt"))
    print(group_tab)
    sink()
}


all_table2 <- function(df,measures, filename, repeats, col_names_list, subheading, table_head, ndig = 4,feat = "Method", group=FALSE){
    df <- mutate_at(df,c(feat, "Measure"), as.factor)
    df <- pivot_wider(df, names_from = Measure, values_from = c("Mean", "S.d."))
    #round to 4 digits
    df <- mutate_if(df, is.numeric, function(x) format(round(x, ndig), nsmall = ndig))
    #add info in brackets
    mean_col_names <- paste0("Mean_",measures)
    sd_col_names <- paste0("S.d._",measures)
    #name of mean column
    for (i in 1:length(measures)){
        df[[measures[i]]] <- sd_brackets(df[[mean_col_names[i]]],
                                df[[sd_col_names[i]]], repeats, measures[i])
    }
    
    df_all <- subset(df, select = c("Method", "Choice", measures))
    df_overall <- subset(df, select = c(feat, measures))
    #measures <- ifelse(measures == "CSR", "CS", measures)
    colnames(df_all) <- c(feat, "Choice", measures)
    #colnames(df_overall) <- c(feat, measures)
    colnames(df_overall) <- col_names_list

    #subset only overall results and remove extras
    #make the largest values down a column bold
    # for(i in 1:length(measures)){
    #     max_id <- max_index(df_overall[[measures[i]]])
    #     if(measures[i]=="Error"){
    #         max_id <- min_index(df_overall[[measures[i]]])
    #     }
    #     cond <- (1:length(df_overall[,1][[1]]))==max_id
    #     df_overall[[measures[i]]] <- text_spec(df_overall[[measures[i]]] , "latex", 
    #         bold = cond)
    # }
    
    df_overall <- as.data.frame(t(df_overall))
    df_overall <- df_overall[-1,]
    colnames(df_overall) <- unique(df[[feat]])
    n <- ncol(df_overall) 
    vec <- c(" " = 1)
    vec[[table_head]] <- n 
    df_overall2 <- kable(df_overall)
    df_overall2 <- kbl(df_overall,booktabs=T,"latex",escape=FALSE)
    df_overall2 <- add_header_above(df_overall2, header = vec)
    df_overall2 <- kable_styling(add_footnote(df_overall2, subheading, escape = FALSE))
    #save
    sink(paste0(filename, "_overall.txt"))
    print(df_overall2)
    sink()

    #collapse method rows
    if (group){
        group_tab <- kable_styling(kbl(df_overall, booktabs = T,"latex"))
    #select groups
        feats <- unique(df_overall[[feat]])
        for(i in 1:length(feats)){
            group_tab <- pack_rows(group_tab,feats[i], 3*i-2, 3*i)
        }
    sink(paste0(filename, "_all.txt"))
    print(group_tab)
    sink()
    }
    return(df_overall2)
}



measure_table <- function(df, measure_val, filename, repeats, col_names, subheading, table_head, ndig = 4){
    factors <- unique(df$Factor)
    df <- subset(df,Measure == measure_val)
    df <- mutate_at(df,c( "Measure"), as.factor)
    df <- pivot_wider(df, names_from = Factor, values_from = c("Mean", "S.d."))
        #round to 4 digits
    df <- mutate_if(df, is.numeric, function(x) format(round(x, ndig), nsmall = ndig))
        #add info in brackets
    mean_col_names <- paste0("Mean_",factors)
    sd_col_names <- paste0("S.d._",factors)
        #name of mean column
    for (i in 1:length(factors)){
        df[[as.character(factors[i])]] <- sd_brackets(df[[mean_col_names[i]]],
                                    df[[sd_col_names[i]]], repeats, factors[i])
    }
    c_factors <- as.character(factors)
    df_all <- subset(df, select = c("Method", c_factors))
    df_overall <- subset(df, select = c("Method", c_factors))
    colnames(df_all) <- c("Method", col_names)
    colnames(df_overall) <- c("", col_names)

    #subset only overall results and remove extras
    #make the largest values down a column bold
    for(i in 1:length(factors)){
            max_id <- max_index(df_overall[[c_factors[i]]])
            cond <- (1:length(df_overall[,1][[1]]))==max_id
            df_overall[[c_factors[i]]] <- text_spec(df_overall[[c_factors[i]]] , "latex", 
                bold = cond)
    }
    n <- ncol(df_overall) - 1
    vec <- c(" " = 1)
    vec[[table_head]] <- n 
    df_overall2 <- kable(df_overall)
    df_overall2 <- kbl(df_overall,booktabs=T, "latex",escape = FALSE)
    df_overall2 <- add_header_above(df_overall2, header = vec)
    df_overall2 <- kable_styling(add_footnote(df_overall2, subheading, escape = FALSE))
        #save
    sink(paste0(filename, "_overall.txt"))
    print(df_overall2)
    sink()
    return(list("table" = df_all, "latex" = df_overall2))
}

combine_tables <- function(tables, filename, measures, table_head, subheading, col_names_tables){
    full_table <- rbind.data.frame(tables[[1]]$table, tables[[2]]$table,
         tables[[3]]$table, tables[[4]]$table, tables[[5]]$table)
    #collapse method rows
    n <- ncol(full_table) - 1
    n2 <- length(unique(full_table[["Method"]]))
    colnames(full_table) <- c("", col_names_tables)
    group_tab <- kable_styling(kbl(full_table, booktabs = T,"latex"))
    vec <- c(" " = 1)
    vec[[table_head]] <- n 
    group_tab <- add_header_above(group_tab, header = vec)
    group_tab <- kable_styling(add_footnote(group_tab, subheading, escape = FALSE))
    #select groups
    for(i in 1:length(measures)){
        group_tab <- pack_rows(group_tab, measures[i], n2*i-3, n2*i)
    }

    sink(paste0(filename, "_all.txt"))
    print(group_tab)
    sink()
    return(group_tab)
}

corr_table <- function(results, filename){
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
  #latex_tab <- kable_styling(add_footnote(latex_tab, subheading))
   #save
  sink(paste0(filename, "corr_tab.txt"))
  print(latex_tab)
  sink()
  return(latex_tab)
}