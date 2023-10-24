bedtools intersect -a ../../DRF_PaperApp/data/Illumina100_CHM13/Updated_output_07_17_2023/originalADSP.Illumina_OriginalADSP.T2T_CHM13_v2.0.low_mapq-merged.bed -b ../../DRF_PaperApp/data/Illumina250_CHM13/Updated_output_07_17_2023/illuminaRL250.IlluminaRL250.T2T_CHM13_v2.0.low_mapq-merged.bed > Illumina100_vs_Illumina250.CHM13.low-mapq.bed
bedtools intersect -a ../../DRF_PaperApp/data/PacBio_CHM13/Updated_output_07_17_2023/PacBio.PacBio.All1KG.T2T_CHM13_v2.0.low_depth-merged.bed -b ../../DRF_PaperApp/data/ONT_CHM13/TenSamples_10_17_2023/ONT.ONT_1KG.T2T_CHM13_v2.0.low_depth-merged.bed > PacBio_vs_ONT.CHM13.low-depth.bed
bedtools intersect -wao -a ../../DRF_PaperApp/data/PacBio_CHM13/Updated_output_07_17_2023/PacBio.PacBio.All1KG.T2T_CHM13_v2.0.low_depth-merged.bed -b ../../DRF_PaperApp/data/ONT_CHM13/TenSamples_10_17_2023/ONT.ONT_1KG.T2T_CHM13_v2.0.low_depth-merged.bed > PacBio_vs_ONT.CHM13.low-depth_allPacBio.bed
bedtools intersect -wao -b ../../DRF_PaperApp/data/PacBio_CHM13/Updated_output_07_17_2023/PacBio.PacBio.All1KG.T2T_CHM13_v2.0.low_depth-merged.bed -a ../../DRF_PaperApp/data/ONT_CHM13/TenSamples_10_17_2023/ONT.ONT_1KG.T2T_CHM13_v2.0.low_depth-merged.bed > PacBio_vs_ONT.CHM13.low-depth_allONT.bed
bedtools intersect -wao -a ../../DRF_PaperApp/data/Illumina100_CHM13/Updated_output_07_17_2023/originalADSP.Illumina_OriginalADSP.T2T_CHM13_v2.0.low_mapq-merged.bed -b ../../DRF_PaperApp/data/Illumina250_CHM13/Updated_output_07_17_2023/illuminaRL250.IlluminaRL250.T2T_CHM13_v2.0.low_mapq-merged.bed > Illumina100_vs_Illumina250.CHM13.low-mapq_allIllumina100.bed
bedtools intersect -wao -b ../../DRF_PaperApp/data/Illumina100_CHM13/Updated_output_07_17_2023/originalADSP.Illumina_OriginalADSP.T2T_CHM13_v2.0.low_mapq-merged.bed -a ../../DRF_PaperApp/data/Illumina250_CHM13/Updated_output_07_17_2023/illuminaRL250.IlluminaRL250.T2T_CHM13_v2.0.low_mapq-merged.bed > Illumina100_vs_Illumina250.CHM13.low-mapq_allIllumina250.bed
bedtools intersect -wao -b Illumina100_vs_Illumina250.CHM13.low-mapq.bed -a PacBio_vs_ONT.CHM13.low-depth.bed > LongRead.DbD_vs_ShortRead.DbM.bed


bedtools intersect -a PacBio_vs_ONT.CHM13.low-depth.bed -b Illumina100_vs_Illumina250.CHM13.low-mapq.bed | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$3-$2}' > PacBio_vs_ONT.CHM13.low-depth.vs.Illumina100.vs.Illumina250.CHM13.low-mapq.bed


