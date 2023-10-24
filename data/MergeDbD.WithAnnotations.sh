awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$3-$2}' PacBio_vs_ONT.CHM13.low-depth.bed> PacBio_vs_ONT.CHM13.low-depth.size.bed
bedtools intersect -wao -a PacBio_vs_ONT.CHM13.low-depth.size.bed -b chm13UniqueCompared2HG38.bed > PacBio_vs_ONT.CHM13.low-depth.size.CHM13Unique.bed
bedtools intersect -b PacBio_vs_ONT.CHM13.low-depth.bed -a chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed | sort -V | bedtools merge -i - -c 4,7 -o distinct,distinct > RMSK_IntersectWith_IntersectedLongReadDbD_Regions.txt
bedtools intersect -wao -a PacBio_vs_ONT.CHM13.low-depth.size.CHM13Unique.bed -b censat.bed > PacBio_vs_ONT.CHM13.low-depth.size.CHM13Unique.censat.bed 
bedtools intersect -wao -a PacBio_vs_ONT.CHM13.low-depth.size.CHM13Unique.censat.bed  -b sedefSegDups_CHM13.bed > PacBio_vs_ONT.CHM13.low-depth.size.CHM13Unique.censat.segDups.bed
