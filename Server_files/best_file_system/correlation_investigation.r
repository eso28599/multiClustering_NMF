#correlation between f score and BiS - phi 4
results <- read.csv("increasing_phi4/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)

#correlation between f score and BiS - phi 3
results <- read.csv("increasing_phi3/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)


#correlation between f score and BiS - phi 3
results <- read.csv("increasing_phi3/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
cor(rec, rel)

#correlation between f score and BiS - phi 2
results <- read.csv("increasing_phi2/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
cor(rec, rel)


#correlation between f score and new BiS
results <- read.csv("results/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
