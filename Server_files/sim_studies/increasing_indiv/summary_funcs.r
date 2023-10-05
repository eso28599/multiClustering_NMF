library("rio")
library(ggplot2)
library("tidyr")
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
make_plot_line <- function(data_frame, cluster, plot_title, sub_title, xtitle, filename){
    p <- ggplot(subset(data_frame, Choice==cluster), aes(x=Method, y=Mean, group=Measure, color=Measure)) + 
    geom_line(stat="identity") +
    geom_point() +
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(0.01))+
     expand_limits(y=c(0,1))+
     scale_color_hue(labels = c("CSR", "F score", "BiS")) +
     scale_y_continuous(breaks = seq(0, 1, length=6))+
     labs(
     title = plot_title,
     subtitle = sub_title,
     x = xtitle
   )
    p <- p +  
        theme(text = element_text(size = 13), legend.position = c(0.8, 0.2))
    ggsave(filename, plot=p, compress=FALSE, device="pdf")
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


all_table2 <- function(df,measures, filename, repeats, ndig = 4,feat = "Method", group=FALSE){
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
    colnames(df_overall) <- c(feat, measures)

    #subset only overall results and remove extras
    #make the largest values down a column bold
    for(i in 1:length(measures)){
        max_id <- max_index(df_overall[[measures[i]]])
        if(measures[i]=="Error"){
            max_id <- min_index(df_overall[[measures[i]]])
        }
        cond <- (1:length(df_overall[,1][[1]]))==max_id
        df_overall[[measures[i]]] <- text_spec(df_overall[[measures[i]]] , "latex", 
            bold = cond)
    }
    df_overall2 <- kbl(df_overall,booktabs=T,"latex",escape=FALSE)
    df_overall2 <- kable_styling(add_footnote(df_overall2, c("Best performance in bold.")))
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
}