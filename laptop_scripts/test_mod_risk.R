library(stringr)
library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(rmda)
theme_set(theme_cowplot())



convert_names <- function(x, the_dict = method_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}

method_names <- read.table("local_info/method_names", stringsAsFactors = F)
disease_names <- read.table("local_info/disease_names", stringsAsFactors = F, sep = "\t")

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
author_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = '\t')

mod_names <- read.table("local_info/mod_splits.txt", stringsAsFactors = F, sep = "\t")

all_authors <- unique(str_split(list.files("test_results/", "class_data"), "_", simplify = T)[,1])
if(any(all_authors == "Xie")){
  all_authors <- all_authors[-which(all_authors == "Xie")]
}

#plot_authors
#plot_vars

keep_factors <- c("days_mod_activity", "smoking_status", "time_tv", "vitamin_d", "cook_veg", "alcohol_intake", "time_summer")
all_mod_df <- list()
all_supp_df <- list()
best_supp_quals <- list()
big_count <- 1
best_count <- 1

for(author in all_authors){

    res <- readRDS(paste0("test_results/", author, ".arr_data.RDS"))
    all_mod_factors <- list()
    single_supp_df <- list()

    #rows are mod risk factor
    #cols are prs risk group
    
    #for(i in which(names(res[[1]]) %in% keep_factors)){
    for(i in 1:length(res[[1]])){
      df <- data.frame(res[["cont_tables"]][[i]])
      sedf <- data.frame(res[["se_cont_tables"]][[i]])
      fishp <- res[["fish_stat"]][res[["fish_stat"]]$mod_factor == names(res[[1]])[i],1]
      fishor <- res[["fish_stat"]][res[["fish_stat"]]$mod_factor == names(res[[1]])[i],2]
  
      single_supp_df[[i]] <- as.data.frame(t(as.data.frame(signif( c(as.numeric(res[[1]][[i]]), fishp, fishor), 3))))
      colnames(single_supp_df[[i]]) <- c("prs lo - mod lo", "prs lo - mod inter", "prs lo - mod hi",
                                    "prs inter - mod lo", "prs inter - mod inter", "prs inter - mod hi",
                                    "prs hi - mod lo", "prs hi - mod inter", "prs hi - mod hi", 
                                    "pval - lo", "pval - inter", "pval - hi", "or - lo", "or - inter", "or - hi")
      single_supp_df[[i]]$mod_factor <- names(res[[1]])[i]
      single_supp_df[[i]]$author <- author
      
      colnames(df) <- c("prs_lo", "prs_inter", "prs_hi")
      df <- melt(df)
      df$mod_group <- rep(as.character(mod_names[which(names(res[[1]])[i] == mod_names[,1]), 3:5]), 3)
      
      colnames(sedf) <- c("prs_lo", "prs_inter", "prs_hi")
      sedf <- melt(sedf)
      df$se <- sedf$value
  
      df$mod_group <- factor(df$mod_group, levels = mod_names[which(names(res[[1]])[i] == mod_names[,1]), 3:5]) 

                       
      
      mod_drop <- df$value[seq(1,9,3)] - df$value[seq(3,9,3)]
      
      the_plot <- ggplot(df, aes(variable, value, fill = mod_group)) +
        geom_bar(position = position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin = value - se, ymax = value + se), position = position_dodge(0.9), width = 0) +
         labs (x = "Polygenic Risk Score Group", y = "Incidence Rate", 
               caption = paste("P-Vals:",  paste(as.character(signif(fishp, 3)), collapse = ", "))) +
        scale_x_discrete(labels = c("Low", "Inter.", "High")) +
        scale_fill_viridis_d(name = gsub("_", "\n", mod_names[which(names(res[[1]])[i] == mod_names[,1]),2]))
  
      #plot(the_plot)
      #ggsave(paste0("mod_factor_plots/", tolower(author), ".", poss_eps[kk], ".", mod_names[i,1], ".png"),
      #       the_plot, "png", height=5, width=6)
      if(all(fishp < 0.05) | sum(fishp < 0.005) == 2 | any(fishp[c(1,3)] < 0.0005)){
        if(abs(mod_drop[3]) > abs(mod_drop[1])){
          if(all(df$value > 0)){
          ggsave(paste0("great_mod_factor_plots/", tolower(author), ".", names(res[[1]])[i], ".png"),
                 the_plot, "png", height=4.5, width=5)
          best_supp_quals[[best_count]] <- c(author, names(res[[1]])[i])
          best_count <- best_count + 1
          }
        }
      }
      
      all_mod_factors[[i]] <- df
      all_mod_factors[[i]]$mod_factor <- names(res[["cont_tables"]])[i]
      all_mod_factors[[i]]$author <- author
  
    }
    
  
    all_mod_df[[big_count]] <- all_mod_factors
    all_supp_df[[big_count]] <- single_supp_df
    big_count <- big_count + 1
  
}



redo_list <- list()
k <- 1
for(i in 1:length(all_supp_df)){
  for(j in 1:length(all_supp_df[[i]])){
    if(!is.null(all_supp_df[[i]][[j]])){
      redo_list[[k]] <- all_supp_df[[i]][[j]]
      k <- k + 1
    }
  }
}

big_df <- do.call("rbind", redo_list)

keep_rows <- do.call("rbind", best_supp_quals)
best_df <- big_df[paste0(big_df$author, "_", big_df$mod_factor) %in% paste0(keep_rows[,1], "_", keep_rows[,2]),]
best_df$disease <- convert_names(best_df$author, author_dict)
best_df$mod_name <- convert_names(best_df$mod_factor, mod_names)

for(i in 1:9){
  best_df[,i] <- signif(best_df[,i], 2)
}

best_df <- best_df[best_df$mod_factor != "cook_veg",]
best_df$mod_name <- gsub("_", " ", best_df$mod_name)
supp1 <- best_df[,c(18, 19, 1:9)]
supp2 <- best_df[,c(18, 19,  10:15)]

supp1$mod_name[supp1$mod_name == "Pieces Fruit Eaten Per Day"] <- "Fruit"
supp1$mod_name[supp1$mod_name == "Hours per Day Watch TV"] <- "TV"
supp1$mod_name[supp1$mod_name == "Processed Meat Intake"] <- "Meat"
supp1$mod_name[supp1$mod_name == "Glasses of Water Per Day"] <- "Water Intake"
supp1$mod_name[supp1$mod_name == "Days Mod. Activity Per Week"] <- "Mod. Activity"
supp1$mod_name[supp1$mod_name == "Hours per Day Driving"] <- "Driving"
supp1$mod_name[supp1$mod_name == "Tbsp Raw Veg. Per Day"] <- "Vegetable"
supp1$mod_name[supp1$mod_name == "Alcohol Consumed"] <- "Alcohol"
supp1$mod_name[supp1$mod_name == "Min. Walked Per Day"] <- "Min. Walked"

supp2$mod_name[supp2$mod_name == "Pieces Fruit Eaten Per Day"] <- "Fruit"
supp2$mod_name[supp2$mod_name == "Hours per Day Watch TV"] <- "TV"
supp2$mod_name[supp2$mod_name == "Processed Meat Intake"] <- "Meat"
supp2$mod_name[supp2$mod_name == "Glasses of Water Per Day"] <- "Water Intake"
supp2$mod_name[supp2$mod_name == "Days Mod. Activity Per Week"] <- "Mod. Activity"
supp2$mod_name[supp2$mod_name == "Hours per Day Driving"] <- "Driving"
supp2$mod_name[supp2$mod_name == "Tbsp Raw Veg. Per Day"] <- "Vegetable"
supp2$mod_name[supp2$mod_name == "Alcohol Consumed"] <- "Alcohol"
supp2$mod_name[supp2$mod_name == "Min. Walked Per Day"] <- "Min. Walked"


write.table(supp1, paste0("supp_tables/mod_factor.vals.txt"), col.names = T, row.names = F, quote = F, sep = "\t") 
write.table(supp2, paste0("supp_tables/mod_factor.sigs.txt"), col.names = T, row.names = F, quote = F, sep = "\t")


# for(uval in unique(big_df$mod_factor)){
#   small_df <- big_df[big_df$mod_factor == uval,]
#   small_df <- small_df[,-which(colnames(small_df) == "mod_factor")]
#   small_df$author <- convert_names(small_df$author, author_dict)
#   small_df <- small_df[,c(10,1:9)]
#   new_df <- as.data.frame(rbind(str_split(colnames(small_df), "-", simplify = T)[,1],
#                                 str_split(colnames(small_df), "-", simplify = T)[,2]), stringsAsFactors = F)
#   colnames(new_df) <- colnames(small_df)
#   small_df <- rbind(new_df, small_df)
#   write.table(small_df, paste0("supp_tables/mod_factor.", uval, ".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
# }




