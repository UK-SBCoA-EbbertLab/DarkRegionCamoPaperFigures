



bedtools intersect -wao -a PacBio_vs_ONT.CHM13.low-depth.bed -b ../../DarkRegionApp/data/RMSK_CHM13/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed | awk -F"\t" '$6 != "."{print $6"\t"$7"\t"$8"\t"$8-$7"\t"$16}' | sort -V | bedtools merge -i - -c 5 -o sum | awk -F"\t" '{a+=$4}END{print a}'


