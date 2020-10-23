
#self diag from 472 to 573
#date attend assessment center from 41 to 43
#medication from 574 to 717

zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f1 > eid.csv
zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f472-573 > self_report_diag.csv
zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f454-471 > self_report_cancer.csv
zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f41-43 > date_assessed.csv
zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f574-717 > medications.csv
