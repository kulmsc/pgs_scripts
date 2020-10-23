# 5 - sex
# 6 - year of birth
# 40 - month of birth
# 1310 - genetic PCs
# 2862 - date of death
# 213 - date lost to follow up
# 2745 - scotland
# 2728 - england
# 2742 - wales

zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f1,5,6,40,1310-1320 > into_base_covars
zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f1,213,2862 > censor_covars
zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f1,2728,2742,2745 > nation_covars
