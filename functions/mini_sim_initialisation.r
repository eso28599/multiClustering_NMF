source('/Users/ellaorme/GitHub/Data_integration/multiClustering_NMF/working_method_new_updates.r')
source('/Users/ellaorme/GitHub/Data_integration/multiClustering_NMF/multi_CoNMTF.r')
# dependence of initial start - svd with random s
# generate and save data
n_views <- 3
#each dataset with 200 samples, same clusters across views
rowClusters1<-list(c(75, 75, 50),c(75, 75, 50),c(75, 75, 50))
#100, 50,250 features respectively
colClusters1<-list(c(30,30,40),c(10, 20, 20),c(100,50,100))

#can change the noise parameter here for level of noise in views
our_data <- multi_view(rowClusters1,colClusters1,noise=5)
data_views <- our_data$data_views
#save data as a file in current directory 
#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(data_views, file = 'data_views.xlsx') 
phi <- matrix(0, nrow = length(data_views), ncol = length(data_views))
phi[1,2:3]<-0.5
nIter <-1000

#initialise storage of results
#number of simulations 
n <- 10
iter_time <- rep(0,n)
iter_error <- rep(0,n)
iter_acc <- rep(0,n)
conv_time <- rep(0,n)
conv_error <- rep(0,n)
conv_acc <- rep(0,n)

for(i in 1:n){
    print(i)
    startTimeIter <- Sys.time()
    data_views_NMTF_iter <- restMultiNMTF_run(Xinput = data_views,
                                        KK = c(3,3,3),
                                        LL = c(3,3,3),
                                        phi = phi,
                                        nIter = nIter)
    endTimeIter <- Sys.time()
    eval_dv2<- evaluate_simulation(X_nmtf = data_views_NMTF_iter,
                                            true_row_clustering = our_data$truth_rows,
                                            true_col_clustering = our_data$truth_cols)
    iter_time[i]<- endTimeIter-startTimeIter
    iter_error[i] <- mean(data_views_NMTF_iter$Error[990:1000])
    iter_acc[i] <- mean(eval_dv$accuracy)

    startTimeConv <- Sys.time()
    data_views_NMTF_conv <- restMultiNMTF_run(Xinput = data_views,
                                        KK = c(3,3,3),
                                        LL = c(3,3,3),
                                        phi = phi)
    endTimeConv <- Sys.time()
    eval_dv<- evaluate_simulation(X_nmtf = data_views_NMTF_conv,
                                            true_row_clustering = our_data$truth_rows,
                                            true_col_clustering = our_data$truth_cols)
    conv_time[i]<- endTimeConv - startTimeConv
    conv_error[i] <- mean(data_views_NMTF_conv$Error[length(data_views_NMTF_conv$Error)])
    conv_acc[i] <- mean(eval_dv$accuracy)
}

#sim study for multi_co against restrictive NMTF
n <- 4
res_error <- rep(0,n)
res_acc <- rep(0,n)
co_error <- rep(0,n)
co_acc <- rep(0,n)
for(i in 1:n){
    print(i)
    data_views_NMTF<- restMultiNMTF_run(Xinput = data_views,
                                        KK = c(3,3,3),
                                        LL = c(3,3,3),
                                        phi = phi,
                                        nIter = nIter)
    eval_dv<- evaluate_simulation(X_nmtf = data_views_NMTF,
                                            true_row_clustering = our_data$truth_rows,
                                            true_col_clustering = our_data$truth_cols)
    res_error[i] <- mean(data_views_NMTF$Error[990:1000])
    res_acc[i] <- mean(eval_dv$accuracy)
    data_views_coNMTF <- mulCoNMTF_run(Xinput = data_views,
                                        KK = c(3,3,3),
                                        LL = c(3,3,3),
                                        phi = phi, nIter = nIter)
    eval_dv<- evaluate_simulation(X_nmtf = data_views_coNMTF,
                                            true_row_clustering = our_data$truth_rows,
                                            true_col_clustering = our_data$truth_cols)
    co_error[i] <- mean(data_views_coNMTF$Error[990:1000])
    co_acc[i] <- mean(eval_dv$accuracy)
}




conv_error
# better to use non-random initialisation 

test_phi <- c(0, 0.01, 0.05,0.1, 0.2, 0.5, 0.65, 0.8, 1.0, 5.0, 10.0, 20.0, 50.0)
tuningAccuracy <- vector("list", length = length(test_phi))
tuningAdjRandValue <- vector("list", length = length(test_phi))
tuningNMIvalue <- vector("list", length = length(test_phi))
tuningError <- vector("list", length = length(test_phi))


for (k in 1:length(test_phi)){
  phi[1,2] <- test_phi[k]
  dataExample_NMTF <- restMultiNMTF_run(Xinput = data_views,
                                            Finput = Finit,
                                            Sinput = Sinit,
                                            Ginput = Ginit,
                                            phi = phi,
                                            xi = xi,
                                            psi = psi,
                                            nIter = nIter)
  eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                 true_row_clustering = true_row_clusterings,
                                                 true_col_clustering = true_col_clusterings)
  tuningAccuracy[[k]] <- eval_measures_dataExample_NMTF$accuracy
  tuningAdjRandValue[[k]] <- eval_measures_dataExample_NMTF$ARI
  tuningNMIvalue[[k]] <- eval_measures_dataExample_NMTF$NMI
  tuningError[[k]] <- dataExample_NMTF$Error
  if (k == 1){
    finalAccuracy_single <- tuningAccuracy[[k]]
    finalAdjRandValue_single <- tuningAdjRandValue[[k]]
    finalNMIvalue_single <- tuningNMIvalue[[k]]
    finalError_single <- tuningError[[k]]
  }
  if (test_phi[k] == 1){
    finalAccuracy_allOnes <- tuningAccuracy[[k]]
    finalAdjRandValue_allOnes <- tuningAdjRandValue[[k]]
    finalNMIvalue_allOnes <- tuningNMIvalue[[k]]
    finalError_allOnes <- tuningError[[k]]
  }
}

make_plot(test1, "Overall", paste0("Overall", " performance with different initialisations"),
                paste0(10, " repetitions"),
                path_to_save)


all_table(test1,measures,"test",4)
#make tables
test1 <- read.csv("/Users/ellaorme/GitHub/Data_integration/multiClustering_NMF/Server_files/NMTF_initialisation/data_1/all_results.csv",row.names=1)
kable(test1,format="html")
test1 <-mutate_at(test1,c("Method", "Cluster", "Measure"), as.factor)
test_tab <-kbl(filter(test1,),format="html")
#collapse_rows(test_tab,columns = 1:2, valign = "middle")

path_to_sim_folder <- "NMTF_initialisation"
cluster <- "Overall"
path_to_save <- paste0(paste0(path_to_sim_folder,"/"), paste0(paste0(cluster, "_diff_init"),".pdf"))
test_tab2 <- pivot_wider(test1,names_from=Measure, values_from= c("Mean","S.d."))
test_tab2 <-mutate_at(test_tab2,c("Method","Cluster"), as.factor)
test_tab2 <- mutate_if(test_tab2, is.numeric, function(x) format(round(x, 4), nsmall = 4))
#new columns 
test_tab2$Accuracy <-sd_brackets(test_tab2$Mean_accuracy,test_tab2$S.d._accuracy)
test_tab2$NMI <-sd_brackets(test_tab2$Mean_NMI,test_tab2$S.d._NMI)
test_tab2$ARI <-sd_brackets(test_tab2$Mean_ARI,test_tab2$S.d._ARI)
test_tab3 <- subset(test_tab2,select=c("Method","Cluster","Accuracy","ARI","NMI"))
kable_styling(collapse_rows(kable(test_tab2),columns = c(1,2),target=1, valign = "middle"))


table1<-kable_styling(kable(collapse_rows_df(test_tab3,Method),"latex"))
table1<-kable(collapse_rows_df(test_tab3,Method),"latex")
grouped <-kbl(subset(test_tab3,select=c("Cluster",measures)), caption = "Group Rows", booktabs = T) %>%
  kable_styling() %>%
  pack_rows("NMTF", 1, 3) %>%
  pack_rows("NMTF - k means", 4, 6)

test_tab4 <-pivot_wider(test_tab3,names_from=Cluster, values_from= c("Accuracy","ARI","NMI"))
colnames(test_tab4) <- c("Method",rep(averages,3))
test5 <-kbl(test_tab4, booktabs = T) %>%
  kable_styling() %>%
  add_header_above(c(" " = 1, "Accuracy" = 3, "ARI" = 3, "NMI" = 3))

max_index <-which.max(as.numeric(sub(" .*", "", as.vector(test_tab4[,2][[1]]))))
sub(".*:", "", string)
test_tab4 <- kable(test_tab4,"latex")
test_tab4[,2] = lapply(as.vector(test_tab4[,2])[[1]], function(x) text_spec(x, "latex", escape=FALSE,
            bold = ifelse(x == max_index,TRUE, FALSE)))
test_tab4[,2][[1]] = text_spec(test_tab4[,2][[1]] , "latex", 
            bold = ifelse(1:2 == max_index,TRUE, FALSE))

check <- kable_styling(kable(test_tab4,"latex"))
mutate(test_tab4, ~cell_spec(.x, 'latex', bold = ifelse(.x == max(.x), TRUE,  FALSE)))
test_tab4 %>% 
  adorn_totals('row') %>% 
  colwise() %>% 
  mutate(a,  ~cell_spec(.x, 'latex', bold = ifelse(.x == max(c_across(var.1:var.4)), TRUE,  FALSE))))

