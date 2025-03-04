
scripts=/share/home/zhangby/data/MK_resequencing/scripts
genegff=/share/home/zhangby/genomefile/rice/longi_genome/OLongi_RD23_PR25_BGI20221215/CX_mosaic/CX.1to12.gff


python ${scripts}/get_gene_bed.py --gene OL1g29489 --region promoter --gff ${genegff} --plink_input  ../CX/merged_snp_filtered --data_type bfile --plink_out Gene_region
python ${scripts}/get_gene_haplo.py --ped Gene_region.ped --output Gene_haplotype.txt