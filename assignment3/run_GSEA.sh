#!/bin/bash
topdir=/Users/andreeduquette/BCB420/student_apostonaut/assignment3
gmt=$topdir/Human_GOBP_AllPathways_no_GO_iea_February_01_2020_symbol.gmt
#rnk=$topdir/MesenvsImmuno_RNASeq_ranks.rnk
rnk=$topdir/sign_hits.rnk
outdir=$topdir/GSEA_outdir
/Users/andreeduquette/GSEA_4.0.3/gsea-cli.sh GSEAPreranked -gmx $gmt -nperm 1000 -set_max 200 -set_min 15 -rnk $rnk -collapse false -scoring_scheme weighted -rpt_label A3_GSEA -out $outdir 
