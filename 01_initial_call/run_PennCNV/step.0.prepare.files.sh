
## on minerva
module load penncnv

# in penncnv working folder

### compile_pfb file
# run following script in penncnv perpared data folder ----------------------/
# format of list_pfb.txt: (each data file in one line)
# 1004044.txt
# 1004045.txt
# 1004046.txt
# 1004047.txt
# 1004048.txt
# 1004049.txt

# format of SNP_pos.txt: 
# Name    Chr     Position
# 1:10001102-G-T  1       10001102
# 1:100011159-T-G 1       100011159
# 1:10002775-GA   1       10002775
# 1:100122796-C-T 1       100122796

compile_pfb.pl -listfile ../list_pfb.txt -snpposfile ../SNP_pos.txt -output ../SNP.pfb

## gcmodel (Instead of SNP_pos.txt)
cal_gc_snp.pl /sc/orga/projects/haok01a/chengh04/shared_genomics_resources/gc5Base_hg19.txt.sorted \
SNP.pfb -output SNP.gcmodel




