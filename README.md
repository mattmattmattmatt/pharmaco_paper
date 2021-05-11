Supporting data for PeerJ publication based on pre-print

https://www.biorxiv.org/content/10.1101/780213v1

Commands run contained in pharmaco.R

PharmGKB annotation data downloaded:
wget https://api.pharmgkb.org/v1/download/file/data/annotations.zip
wget https://api.pharmgkb.org/v1/download/file/data/relationships.zip
wget https://api.pharmgkb.org/v1/download/file/data/drugLabels.zip
wget https://api.pharmgkb.org/v1/download/file/data/clinicalVariants.zip

Functional inference tools run using command line or web portal
VEP: variant_effect_predictor.pl --force_overwrite --cache --offline --species human --fork 2 --canonical -coding_only --poly b --sift b --plugin CADD,whole_genome_SNVs.tsv.gz --plugin dbNSFP,dbNSFP.gz,ALL --dir ./human/GRCh38/vep_index -o vep.out -i vep.in
Phylo: wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP20way/hg38.phyloP20way.wigFix.gz
DANN: wget https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz
MutationPredictor: http://mutpred.mutdb.org/
MutationTaster: http://www.mutationtaster.org/

