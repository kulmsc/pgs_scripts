rm try_single_score/ss
rm try_single_score/rsids
rm try_single_score/geno_list

Author=IMSGC
author=imsgc

for chr in {1..22};do
  cat ../mod_sets/${Author}/${author}.${chr}.clump.1.ss >> try_single_score/ss
  cat ../mod_sets/${Author}/${author}.${chr}.clump.1.ss | cut -f3 > try_single_score/rsids
  bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chr}.bgen -incl-rsids try_single_score/rsids > try_single_score/temp.bgen
  plink2_new --memory 12000 --threads 12 --bgen try_single_score/temp.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.1.sample --keep-fam temp_files/brit_eid --make-bed --out try_single_score/geno.${chr}
  echo try_single_score/geno.${chr} >> try_single_score/geno_list
done

plink --merge-list try_single_score/geno_list --make-bed --out try_single_score/all
plink --bfile try_single_score/all --keep-allele-order --score try_single_score/ss 3 4 7 no-sum --out try_single_score/score.${author}.clump.1
