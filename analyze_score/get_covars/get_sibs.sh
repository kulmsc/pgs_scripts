# defs from https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-69927-7/MediaObjects/41598_2020_69927_MOESM1_ESM.pdf
cat ~/athena/ukbiobank/qc/ukb47137_rel_s488292.dat | awk '$4 > 0.0012 {print $0}' | awk '$5 > 0.176 {print $0}' | cut -f1,2 -d' ' > sibs

