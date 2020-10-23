
all_score <- readRDS("final_scores/all_score.1.RDS")
template_score <- read.table("template.profile", stringsAsFactors=F, header=T)

fileind <- which(colnames(all_score) == "bentham.1.clump")
new_score <- read.table("try_single_score/score.bentham.clump.1.profile", stringsAsFactors=F, header=T)
new_score <- new_score[order(new_score[,1])[rank(template_score[,1])],]
cor(all_score[,fileind], new_score$SCORE)
cor(rank(new_score$SCORE, ties.method="average"), rank(all_score[,fileind], ties.method="average"))
mean_dif <- mean(abs(rank(new_score$SCORE, ties.method="average") - rank(all_score[,fileind], ties.method="average")))
max_dif <- max(abs(rank(new_score$SCORE, ties.method="average") - rank(all_score[,fileind], ties.method="average")))

###

fileind <- which(colnames(all_score) == "christophersen.1.clump")
new_score <- read.table("try_single_score/score.christophersen.clump.1.profile", stringsAsFactors=F, header=T)
new_score <- new_score[order(new_score[,1])[rank(template_score[,1])],]
cor(all_score[,fileind], new_score$SCORE)
cor(rank(new_score$SCORE, ties.method="average"), rank(all_score[,fileind], ties.method="average"))
mean_dif <- mean(abs(rank(new_score$SCORE, ties.method="average") - rank(all_score[,fileind], ties.method="average")))
max_dif <- max(abs(rank(new_score$SCORE, ties.method="average") - rank(all_score[,fileind], ties.method="average")))

###

fileind <- which(colnames(all_score) == "imsgc.1.clump")
new_score <- read.table("try_single_score/score.imsgc.clump.1.profile", stringsAsFactors=F, header=T)
new_score <- new_score[order(new_score[,1])[rank(template_score[,1])],]
cor(all_score[,fileind], new_score$SCORE)
cor(rank(new_score$SCORE, ties.method="average"), rank(all_score[,fileind], ties.method="average"))
mean_dif <- mean(abs(rank(new_score$SCORE, ties.method="average") - rank(all_score[,fileind], ties.method="average")))
max_dif <- max(abs(rank(new_score$SCORE, ties.method="average") - rank(all_score[,fileind], ties.method="average")))
