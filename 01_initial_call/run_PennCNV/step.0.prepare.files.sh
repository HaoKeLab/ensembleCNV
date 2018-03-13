

module load penncnv

### compile_pfbfile
compile_pfb.pl -listfile ../list_pfb.txt -snpposfile ../SNP_pos.txt -output ../SNP.pfb

## gcmodel (Instead of SNP_pos.txt)
cal_gc_snp.pl /sc/orga/projects/haok01a/chengh04/shared_genomics_resources/gc5Base_hg19.txt.sorted SNP.pfb -output SNP.gcmodel




