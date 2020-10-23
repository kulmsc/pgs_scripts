library(vroom)

defs <- read.table("../descript_defs/covar_defs_singlecol", stringsAsFactors=F, header = T)
preg_codes <- strsplit(defs[which(defs[,1] == "pregnant"),2], ",")[[1]]
use_codes <- c(defs$single_col[-which(defs[,1] == "pregnant")], preg_codes)

phen <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb26867.csv.gz", delim = ","))
sub_phen_1 <- phen[,colnames(phen) %in% c(use_codes, "eid")]

phen <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb33822.csv.gz", delim = ","))
sub_phen_2 <- phen[,colnames(phen) %in% c(use_codes, "eid")]

phen <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb42385.csv.gz", delim = ","))
sub_phen_3 <- phen[,colnames(phen) %in% c(use_codes, "eid")]

sub_phen_1 <- sub_phen_1[sub_phen_1$eid %in% sub_phen_2$eid & sub_phen_1$eid %in% sub_phen_3$eid,]
sub_phen_2 <- sub_phen_2[sub_phen_2$eid %in% sub_phen_1$eid & sub_phen_2$eid %in% sub_phen_3$eid,]
sub_phen_3 <- sub_phen_3[sub_phen_3$eid %in% sub_phen_2$eid & sub_phen_3$eid %in% sub_phen_1$eid,]

sub_phen_2 <- sub_phen_2[order(sub_phen_2$eid)[rank(sub_phen_1$eid)],]
sub_phen_3 <- sub_phen_3[order(sub_phen_3$eid)[rank(sub_phen_1$eid)],]

new_phen <- cbind(sub_phen_1, sub_phen_2[,-1, drop = F], sub_phen_3[,-1, drop = F])
preg_ans <- (new_phen[["3140-0.0"]] == 1 & !is.na(new_phen[["3140-0.0"]])) | (!is.na(new_phen[["2754-0.0"]]))
preg_ans <- preg_ans*1

inds <- 1:nrow(defs)[-10]
for(i in inds){
  colnames(new_phen)[colnames(new_phen) == defs[i,2]] <- defs[i,1]
}
new_phen <- new_phen[,!grepl("0.0", colnames(new_phen))]
new_phen$pregnant <- preg_ans

##################### CLEAN BY HAND ##############################

new_phen$alcohol[is.na(new_phen$alcohol) | new_phen$alcohol == -3] <- 4
new_phen$smoking[is.na(new_phen$smoking) | new_phen$smoking == -3] <- 0
new_phen <- new_phen[,-which(colnames(new_phen) == "alcohol.1")]
new_phen <- new_phen[,-which(colnames(new_phen) == "smoking.1")]

new_phen$bmi[is.na(new_phen$bmi)] <- mean(new_phen$bmi, na.rm = T)

new_phen$exercise[is.na(new_phen$exercise) | new_phen$exercise == -3 | new_phen$exercise == -1] <- 2

new_phen$use_of_sun_protection[is.na(new_phen$use_of_sun_protection) | new_phen$use_of_sun_protection == -3 | new_phen$use_of_sun_protection == -1] <- 2

new_phen$age_menarche[is.na(new_phen$age_menarche) | new_phen$age_menarche == -1 | new_phen$age_menarche == -3] <- mean(new_phen$age_menarche[!(is.na(new_phen$age_menarche) | new_phen$age_menarche == -1 | new_phen$age_menarche == -3 | new_phen$sex == 1)])
new_phen$age_menarche[new_phen$sex == 1] <- NA

new_phen$hormone_replacement_therapy[is.na(new_phen$hormone_replacement_therapy) | new_phen$hormone_replacement_therapy == -1 | new_phen$hormone_replacement_therapy == -3] <- 0
new_phen$hormone_replacement_therapy[new_phen$sex == 1] <- NA

new_phen$had_menopause[is.na(new_phen$had_menopause) | new_phen$had_menopause == 2 | new_phen$had_menopause == 3] <- 0
new_phen$had_menopause[new_phen$sex == 1] <- NA
new_phen$age_menopause[new_phen$had_menopause == 0] <- NA
new_phen$age_menopause[new_phen$sex == 0 & is.na(new_phen$age_menopause)] <- mean(new_phen$age_menopause, na.rm = T)
new_phen$age_menopause[new_phen$sex == 1] <- NA

new_phen$ever_used_oral_contraceptive[is.na(new_phen$ever_used_oral_contraceptive) | new_phen$ever_used_oral_contraceptive == -3 | new_phen$ever_used_oral_contraceptive == -1] <- 0
new_phen$ever_used_oral_contraceptive[new_phen$sex == 1] <- NA
new_phen$age_started_oral_contraceptive[new_phen$ever_used_oral_contraceptive == 0 | new_phen$age_started_oral_contraceptive < 10] <- NA
new_phen$age_started_oral_contraceptive[new_phen$sex == 0 & is.na(new_phen$age_started_oral_contraceptive)] <- mean(new_phen$age_started_oral_contraceptive, na.rm = T)
new_phen$age_started_oral_contraceptive[new_phen$sex == 1] <- NA

new_phen$epstein_barr_virus[is.na(new_phen$epstein_barr_virus)] <- 0.5

##########################################################################

write.table(new_phen, "covar_data/single_col_covars", row.names = F, col.names = T, quote = F, sep = "\t")
