# compare overlap of som mut calls

# cd /n/irizarryfs01_backed_up/kkorthauer/MCC/PREPROCESS/DNA

# compare mutect-results/DFCI-5367-T-01_muts.vcf (matched norm) and 
# mutect-results-test2/DFCI-5367-T-01_muts.vcf (unmatched norm)

####################################################################
## overlap between the two
####################################################################
awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' \
  mutect-results/DFCI-5367-T-01_muts.vcf \
  mutect-results-test2/DFCI-5367-T-01_muts.vcf > \
  mutect-results-test2/DFCI-5367-T-01_overlap_with_matchednorm.vcf


####################################################################
## only detected with matched norm
####################################################################
awk 'NR==FNR{a[$1$2];next}!($1$2 in a)' \
  mutect-results-test2/DFCI-5367-T-01_muts.vcf \
  mutect-results/DFCI-5367-T-01_muts.vcf > \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm.vcf


####################################################################
## only detected with altnorm
####################################################################
awk 'NR==FNR{a[$1$2];next}!($1$2 in a)' \
  mutect-results/DFCI-5367-T-01_muts.vcf \
  mutect-results-test2/DFCI-5367-T-01_muts.vcf > \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm.vcf


####################################################################
## matched norm-only call stats in altnorm
####################################################################
# add header line first
head -2 mutect-results-test2/DFCI-5367-T-01_call_stats.txt > \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm_call_stats.txt

awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm.vcf \
  mutect-results-test2/DFCI-5367-T-01_call_stats.txt >> \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm_call_stats.txt
grep -Ev '^(GL|NC)' \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm_call_stats.txt > \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm_call_stats2.txt
awk '{print $50}' \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm_call_stats2.txt > \
  mutect-results-test2/DFCI-5367-T-01_only_with_matchednorm_failurereasons.txt 


####################################################################
## altnorm-only call stats in matched norm
####################################################################
# add header line first
head -2 mutect-results/DFCI-5367-T-01_call_stats.txt > \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm_call_stats.txt

awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm.vcf \
  mutect-results/DFCI-5367-T-01_call_stats.txt >> \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm_call_stats.txt
grep -Ev '^(GL|NC)' \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm_call_stats.txt > \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm_call_stats2.txt
awk '{print $50}' \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm_call_stats2.txt > \
  mutect-results-test2/DFCI-5367-T-01_only_with_altnorm_failurereasons.txt 


# compare mutect-results-test3/DFCI-5367-T-01_muts.vcf (matched norm) and 
# mutect-results-test2/DFCI-5367-T-01_muts.vcf (unmatched norm)

####################################################################
## overlap between the two
####################################################################
awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' \
  mutect-results-test3/DFCI-5367-T-01_muts.vcf \
  mutect-results-test2/DFCI-5367-T-01_muts.vcf > \
  mutect-results-test2/DFCI-5367-T-01_overlap_with_2diff_altnorm.vcf

