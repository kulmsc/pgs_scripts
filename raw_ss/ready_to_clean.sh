#Reduce the full list of RSIDs to those with INFO score greater than 0.9
cat ~/athena/ukbiobank/qc/imputed/ukb_mfi_all.txt | awk '$8 > 0.9 {print $0}' > common_files/impute_rsids

#Extract all of the duplicated RSIDs
cat ~/athena/ukbiobank/qc/imputed/ukb_mfi_all.txt | cut -f2 | sort | uniq -d > common_files/dup_ids

#Remove the duplicated RSIDs from the full list of (INFO approved) SNPs
cat common_files/impute_rsids | fgrep -w -v -f common_files/dup_ids > common_files/temp

mv common_files/temp common_files/impute_rsids
