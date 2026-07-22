#/!/bash

awk 'NR>1 {
  split($3, counts, "/")
  het = counts[2]
  total = counts[1] + counts[2] + counts[3]
  if (total > 0) {
    heterozygosity = het / total
    print $1, $2, heterozygosity
  }
}' Pmart_site_hardy-het.hwe > snp_heterozygosity.hwe