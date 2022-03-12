# pgs_scripts

These scripts align with the manuscript entitled: "The Impact of Polygenic Risk Score Calculation on Prediction Accuracy"

The process begins with twenty-three sets of genome wide summary statistics from the GWAS Catalog located within raw_ss.
Quality control procedures are applied to these raw summary statistics, matching the variants to those of the UK Biobank.
Next each set of summary statistics are adjusted within adjust_ss. Many different methods are applied, all of which attempt to select variants and/or adjust effect sizes such that the adjusted summary statistics will create the most accurate polygenic risk score.
To determine how well each adjustment method actually works we form/calculate the polygenic risk scores within do_score.  Specifically, the scores are calculated upon the UK Biobank, nearly 500,000 individuals.
Following, the polygenic risk scores are analyzed within the aptly named analyze_score directory.  Various statistics are computed such as AUC and odds ratio by comparing the score, adjusted by age, sex and other covariates, against electronic health record defined disease labels.
Lastly, the accuracy statistics are visualized with the scripts within the laptop_scripts directory.

All of these scripts are far better described within the manuscript - or within the pgs_book directory.
