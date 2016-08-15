/lore/tyler.kent/Software/hapcut/extractHAIRS --VCF ./Hapcut/Data/VCFs/108.vcf --bam ./Hapcut/Data/Bams/108.sorted.bam --maxIS 600 > ./Hapcut/Results/108.hairs 2> hairs.err

/lore/tyler.kent/Software/hapcut/HAPCUT --fragments ./Hapcut/Results/108.hairs --maxmem 12000 --VCF ./Hapcut/Data/VCFs/108.vcf --output ./Hapcut/Results/108.hapcut --maxiter 100 > hapcut.log 2> hapcut.err
