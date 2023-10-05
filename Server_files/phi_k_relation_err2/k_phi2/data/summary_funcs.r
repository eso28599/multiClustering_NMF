library("rio")
library(ggplot2)
library("tidyr")
#import as list of matrices instead of as a list
import_matrix <- function(filename){
    return(lapply(import_list(filename), function(x) as.matrix(x)))
}
calc_averages <- function(mat_res){return(c(rowMeans(mat_res), mean(mat_res)))}

get_averages <- function(list_avgs, data_name, filename, measures, repeat_no){
    results <- import_matrix(paste0(paste0(data_name,"/"), paste0(filename,".xlsx")))
    for (j in 1:length(measures)){
        #list_avgs[[j]][repeat_no, ] <- calc_averages(results[[j]])
        #list_avgs[[j]][repeat_no, ] <- mean(results[[j]])
        list_avgs[[j]][repeat_no, ] <- mean(results[[measures[j]]])
    }
    return(list_avgs)
    }

avg_over_repeats <- function(avgs, i, measure, name,name2){
    rep_avgs <- colMeans(avgs[[i]])
    rep_sd <- apply(avgs[[i]],2,sd)
    #return(cbind(rep(name,3), c("Row", "Column", "Overall"), rep(measure,3), rep_avgs, rep_sd))
    return(cbind(name, name2, measure, rep_avgs, rep_sd))
}

make_plot <- function(data_frame, cluster, plot_title, sub_title, filename){
    p <- ggplot(subset(data_frame, Cluster==cluster), aes(x=Measure, y=Mean, fill=Method)) + 
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
make_plot_k <- function(data_frame, cluster, plot_title, sub_title, filename){
    p <- ggplot(subset(data_frame, Cluster==cluster), aes(x=k, y=Mean, group=Measure, color=Measure)) + 
    geom_line(stat="identity", position=position_dodge(0.05)) +
    geom_point()+
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(0.05))+
     expand_limits(y=c(0,1))+
     scale_y_continuous(breaks = seq(0, 1, length=6))+
     labs(
     title = plot_title,
     subtitle = sub_title
   )
    p<-p+  scale_fill_brewer() + theme(text = element_text(size = 16),legend.position = c(0.8,0.2))    
    ggsave(filename,plot=p, compress=FALSE,device="pdf")
}

make_plot_phi <- function(data_frame, plot_title, sub_title, filename){
    p <- ggplot(data_frame, aes(x=Phi, y=Mean, group=k, color=k)) + 
    geom_line(stat="identity", position=position_dodge(0.05)) +
    geom_point()+
     geom_errorbar(aes(ymin=Mean-S.d., ymax=Mean+S.d.), width=.2,
                 position=position_dodge(0.05))+
     expand_limits(y=c(0,1))+
     scale_y_continuous(breaks = seq(0, 1, length=6))+
     labs(
     title = plot_title,
     subtitle = sub_title
   )
    p<-p+  scale_fill_brewer() + theme(text = element_text(size = 16),legend.position = c(0.8,0.2))    
    ggsave(filename,plot=p, compress=FALSE,device="pdf")
}

sd_brackets <- function(mean,sd){
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

#updated all_table
all_table <- function(df,phi,filename, ndig=4,feat="Method", group=TRUE){

    df <- mutate_at(df,c(feat, "Measure"), as.factor)
    df <- pivot_wider(df, names_from = Phi, values_from = c("Mean", "S.d."))
    #round to 4 digits
    df <- mutate_if(df, is.numeric, function(x) format(round(x, ndig), nsmall = ndig))
    #
    #add info in brackets
    mean_col_names <- paste0("Mean_",phi)
    sd_col_names <- paste0("S.d._",phi)
    #name of mean column
    for (i in 1:length(phi)){
        df[[as.character(phi[i])]] <-sd_brackets(df[[mean_col_names[i]]],df[[sd_col_names[i]]])
    }
    df_all <- subset(df,select=c("k", as.character(phi)))
    colnames(df_all) <- c(feat,as.character(phi))

    #collapse method rows
    if(group){
        group_tab <- kable_styling(kbl(df_all, booktabs = T,"latex"))
    #select groups
        feats <- unique(df[[feat]])
        for(i in 1:length(feats)){
            group_tab <- pack_rows(group_tab,feats[i], 3*i-2, 3*i)
        }
    sink(paste0(filename, "_all.txt"))
    print(group_tab)
    sink()
    }
    #subset only overall results and remove extras
    df_overall <- subset(df,select=c(feat,as.character(phi)))
    #make the largest values down a column bold
    for(i in 1:length(phi)){
        max_id <- max_index(df_overall[[as.character(phi)]])
        cond <- (1:length(df_overall[,1][[1]]))==max_id
        df_overall[[as.character(phi)]] <- text_spec(df_overall[[as.character(phi)]] , "latex", 
            bold = cond)
    }
    df_overall <- kbl(df_overall,booktabs=T,"latex",escape=FALSE)
    df_overall <- kable_styling(add_footnote(df_overall, c("Best performance in bold.")))
    #save
    sink(paste0(filename, "_overall.txt"))
    print(df_overall)
    sink()
    
}