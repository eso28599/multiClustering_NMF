#load in data
breast <- vector("list")
breast[["Gene"]] <- read.table("tcga_data/Breast/BREAST_Gene_Expression.txt")
breast[["Methy"]] <- read.table("tcga_data/Breast/BREAST_Methy_Expression.txt")
breast[["Mirna"]] <- read.table("tcga_data/Breast/BREAST_Mirna_Expression.txt")
#load in data
lung <- vector("list")
lung[["Gene"]] <- read.table("tcga_data/Lung/LUNG_Gene_Expression.txt")
lung[["Methy"]] <- read.table("tcga_data/Lung/LUNG_Methy_Expression.txt")
lung[["Mirna"]] <- read.table("tcga_data/Lung/LUNG_Mirna_Expression.txt")
#kidney
kidney <- vector("list")
kidney[["Gene"]] <- read.table("tcga_data/Kidney/KIDNEY_Gene_Expression.txt")
kidney[["Methy"]] <- read.table("tcga_data/Kidney/KIDNEY_Methy_Expression.txt", col.names = 1:123)
kidney[["Mirna"]] <- read.table("tcga_data/Kidney/KIDNEY_Mirna_Expression.txt")
library("stringr")
kidney[["Methy"]] <- kidney[["Methy"]][, -1]
col_names <- str_replace_all(read.table("tcga_data/Kidney/col_names.txt")[1,], "-", ".")
#replace - with .
colnames(kidney[["Methy"]]) <- c( col_names)
#kidney[["Methy"]] <- kidney[["Methy"]][, -1]
new_names <- sub("LLL.*", "", rownames(kidney[["Gene"]]))
not_unique <- c()
for(i in 1:length(new_names)){
    tot <- sum(new_names == new_names[i])
    if(tot > 1){
        not_unique <- c(not_unique, new_names[new_names==new_names[i]])
    }
}

kidney[["Gene"]] <- kidney[["Gene"]][!(new_names %in% not_unique), ]
rownames(kidney[["Gene"]]) <- sub("LLL.*", "", rownames(kidney[["Gene"]]))

#ensure for each cancer individuals are in the correct order
colnames(breast[["Gene"]]) <- substr(colnames(breast[["Gene"]]), 1, 12)
colnames(breast[["Methy"]]) <- substr(colnames(breast[["Methy"]]), 1, 12)
colnames(breast[["Mirna"]]) <- substr(colnames(breast[["Mirna"]]), 1, 12)
cols <- colnames(breast[["Gene"]])
breast[["Methy"]] <- (breast[["Methy"]])[, cols]
breast[["Mirna"]] <- breast[["Mirna"]][, cols]

#ensure for each cancer individuals are in the correct order
colnames(kidney[["Gene"]]) <- substr(colnames(kidney[["Gene"]]), 1, 12)
colnames(kidney[["Methy"]]) <- substr(colnames(kidney[["Methy"]]), 1, 12)
colnames(kidney[["Mirna"]]) <- substr(colnames(kidney[["Mirna"]]), 1, 12)
cols <- colnames(kidney[["Gene"]])
kidney[["Methy"]] <- (kidney[["Methy"]])[, cols]
kidney[["Mirna"]] <- kidney[["Mirna"]][, cols]

#ensure for each cancer individuals are in the correct order
colnames(lung[["Gene"]]) <- substr(colnames(lung[["Gene"]]), 1, 12)
colnames(lung[["Methy"]]) <- substr(colnames(lung[["Methy"]]), 1, 12)
colnames(lung[["Mirna"]]) <- substr(colnames(lung[["Mirna"]]), 1, 12)
cols <- colnames(lung[["Gene"]])
lung[["Methy"]] <- (lung[["Methy"]])[, cols]
lung[["Mirna"]] <- lung[["Mirna"]][, cols]



data_process <- function(cancer){
    omics <- c("Gene", "Methy", "Mirna")
    data <- vector("list")  
    for(omic in omics){
        file_path <- paste0("tcga_data/", cancer, "/", toupper(cancer), "_", omic, "_Expression.txt")
        if((cancer=="Kidney")&(omic == "Methy")){
            data[[omic]] <- read.table(file_path, col.names = 1:123)
        }else{
            data[[omic]] <- read.table(file_path)
            colnames(data[[omic]]) <- substr(colnames(data[[omic]]), 1, 12)
        }  
    }
    if(cancer == "Kidney"){
        data[["Methy"]] <- data[["Methy"]][, -1]
        col_names <- str_replace_all(read.table("tcga_data/Kidney/col_names.txt")[1,], "-", ".")
        col_names <- substr(col_names, 1, 12)
        #replace - with .
        colnames(data[["Methy"]]) <- c(col_names)
        #kidney[["Methy"]] <- kidney[["Methy"]][, -1]
        new_names <- sub("LLL.*", "", rownames(data[["Gene"]]))
        not_unique <- c()
        for(i in 1:length(new_names)){
            tot <- sum(new_names == new_names[i])
            if(tot > 1){
                not_unique <- c(not_unique, new_names[new_names==new_names[i]])
            }
        }
    data[["Gene"]] <- data[["Gene"]][!(new_names %in% not_unique), ]
    rownames(data[["Gene"]]) <- sub("LLL.*", "", rownames(data[["Gene"]]))
    }
    cols <- colnames(data[["Gene"]])
    data[["Methy"]] <- (data[["Methy"]])[, cols]
    data[["Mirna"]] <- data[["Mirna"]][, cols]

    #check portion of NA's
    return(data)
}

library(stringr)

breast <- data_process("Breast")
kidney <- data_process("Kidney")
lung <- data_process("Lung")
data_combine <- function(cancer1, cancer2, cancer3){

}
#select common variables 
#dataset
dataset <- vector("list")
true_rows <- vector("list")
#breast
for(omic in c("Gene", "Methy", "Mirna")){
    vars_breast <- rownames(breast[[omic]])
    vars_kidney <- rownames(kidney[[omic]])
    vars_lung <- rownames(lung[[omic]])
    vars_all <- intersect(intersect(vars_breast, vars_kidney), vars_lung)
    data_breast <- t(breast[[omic]][vars_all, ])
    data_kidney <- t(kidney[[omic]][vars_all, ])
    data_lung <- t(lung[[omic]][vars_all, ])
    dataset[[omic]] <- rbind(data_breast, data_kidney, data_lung)
    true_rows[[omic]] <- matrix(0, 
                 nrow(dataset[[omic]]), 3)
    nb <- nrow(data_breast)
    nk <- nrow(data_kidney)
    nl <- nrow(data_lung)
    true_rows[[omic]][1:nb,1] <- 1
    true_rows[[omic]][(nb+1):(nb+nk),2] <- 1
    true_rows[[omic]][(nb+nk+1):(nb+nk+nl),3] <- 1
}

#shuffle
new_ind <- sample(dim(dataset[["Gene"]])[1])
for(omic in c("Gene", "Methy", "Mirna")){
    dataset[[omic]] <- dataset[[omic]][new_ind, ]
    true_rows[[omic]] <- true_rows[[omic]][new_ind, ]
    #add minimum 
    #dataset[[omic]] <- apply(dataset[[omic]], 2, function(x) x+min(x))

    #dataset[[omic]] <- apply(dataset[[omic]], 2, function(x) x/sum(x))
    print(omic)
}

for(omic in c("Gene", "Methy", "Mirna")){
    #add minimum 
    dataset[[omic]] <- apply(dataset[[omic]], 2, function(x) x + abs(min(x)))

    dataset[[omic]] <- apply(dataset[[omic]], 2, function(x) x/sum(x))
    print(omic)
}
test <- apply(dataset[[omic]], 2, function(x)  min(x))
source("resNMTF_funcs.r")
sd_vec <- apply(dataset[["Gene"]], 2, sd)
dataset[["Gene"]] <- dataset[["Gene"]][, sd_vec>0.0015]
#subset genes for methylation
sd_vec <- apply(dataset[["Methy"]], 2, sd)
dataset[["Methy"]] <- dataset[["Methy"]][, sd_vec>0.0035]


phi_mat <- matrix(0, 3, 3)
phi_mat[1, c(2, 3)] <- 1
phi_mat[2, c(3)] <- 1

res200 <- restMultiNMTF_run(Xinput = dataset, k_min=3,
            phi = 200*phi_mat, relErr = "Sil") 

#calculate jaccard score
source("evaluation_funcs.r")
mean(evaluate_real(res200$row_clusters, res200$col_clusters,
            true_rows, dataset)$F_score)
source("resNMTF_funcs.r")
res0 <- restMultiNMTF_run(Xinput = dataset, KK=c(3,3,3), LL = c(3,3,3),
            phi = 0*phi_mat, relErr = "Sil") 
mean(evaluate_real(res0$row_clusters, res0$col_clusters,
            true_rows, dataset)$F_score)


res400 <- restMultiNMTF_run(Xinput = dataset, KK=c(3,3,3), LL = c(3,3,3),
            phi = 400*phi_mat, relErr = "Sil") 
mean(evaluate_real(res400$row_clusters, res400$col_clusters,
            true_rows, dataset)$F_score)
res600 <- restMultiNMTF_run(Xinput = dataset, KK=c(3,3,3), LL = c(3,3,3),
            phi = 600*phi_mat, relErr = "Sil") 
mean(evaluate_real(res600$row_clusters, res600$col_clusters,
            true_rows, dataset)$F_score)

res800 <- restMultiNMTF_run(Xinput = dataset, KK=c(3,3,3), LL = c(3,3,3),
            phi = 800*phi_mat, relErr = "Sil") 
mean(evaluate_real(res800$row_clusters, res800$col_clusters,
            true_rows, dataset)$F_score)

source("gfa_funcs.r")
gfa_res <- gfa_apply(dataset, 3)
mean(evaluate_real(gfa_res$row_clusters, gfa_res$col_clusters,
            true_rows, dataset)$F_score)
jaccard_row(res800$row_clusters[[1]], true_rows[[1]])

jaccard_results(res200$row_clusters[[2]], true_rows[[2]])
jaccard_results(res200$row_clusters[[3]], true_rows[[3]])

evaluate_real(res800$row_clusters, res800$col_clusters,
            true_rows, dataset)
evaluate_real(gfa_res$row_clusters, gfa_res$col_clusters,
            true_rows, dataset)

jaccard_row(gfa_res$row_clusters[[1]], true_rows[[1]])
