#!/bin/bash

REFERENCE="../data/reference_genome_hs37d5x/hs37d5x.fa"

temp_hdr="temp_hdr.txt"
echo '##fileformat=VCFv4.2' > $temp_hdr
echo '##INFO=<ID=GERP_SCORE,Number=1,Type=Float,Description="Genomic Evolutionary Rate Profiling (GERP) RS score S (from step 4)">' >> $temp_hdr
echo '##INFO=<ID=GERP_RATE,Number=1,Type=Float,Description="Genomic Evolutionary Rate Profiling (GERP) neutral rate N (from step 2)">' >> $temp_hdr
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> $temp_hdr

list_shards="list_shards.txt"
while IFS= read -r inline; do

	# Process each shard.
	# For each shard, only the nucleotide positions falling in the shard are extracted.

	IFS=' ' read -r -a array <<< "$inline"
	bit1="${array[0]}"
	shard=${bit1:5}
	bit2="${array[1]}"
	chrom=${bit2:5}
	bit3="${array[2]}"
	pos_start_end=${bit3:3}
	IFS='-' read -r -a array2 <<< "$pos_start_end"
	pos_start="${array2[0]}"
	pos_end="${array2[1]}"

	head_amt=$pos_end
	tail_amt=$(( $pos_end - $pos_start + 1 ))
	infile=gerp.chr"${chrom}".maf.rates.bed
	temp_outfile=temp_gerp.shard"${shard}".chrom"${chrom}".pos"${pos_start}"-"${pos_end}".vcf
	outfile=gerp.shard"${shard}".chrom"${chrom}".pos"${pos_start}"-"${pos_end}".vcf

	echo 'processing' $infile 'to produce' $temp_outfile
	head -n "${head_amt}" $infile | tail -n "${tail_amt}" | awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", ".", "\t", "N", "\t", "A", "\t", ".", "\t", ".", "\t", "GERP_SCORE=",$4,";GERP_RATE=",$3}' | grep -v 'GERP_SCORE=0;GERP_RATE=0$' > $temp_outfile

	cat $temp_hdr > $outfile

	while IFS= read -r inline2; do
		IFS=$'\t' read -r -a array2 <<< "$inline2"
		chrom2="${array2[0]}"
		pos2="${array2[1]}"
		id2="${array2[2]}"
		ref2="${array2[3]}"
		alt2="${array2[4]}"
		qual2="${array2[5]}"
		filter2="${array2[6]}"
		info2="${array2[7]}"
		result2=`samtools faidx $REFERENCE $chrom2:$pos2-$pos2`
		for line2 in $result2; do 
			real_ref2=$line2
		done
		if [ "$real_ref2" != "A" ]; then
			new_line="${chrom2}\t${pos2}\t${id2}\t${real_ref2}\tA\t${qual2}\t${filter2}\t${info2}"
			echo -e $new_line >> $outfile
		fi
                if [ "$real_ref2" != "C" ]; then
                        new_line="${chrom2}\t${pos2}\t${id2}\t${real_ref2}\tC\t${qual2}\t${filter2}\t${info2}"
                        echo -e $new_line >> $outfile
                fi
                if [ "$real_ref2" != "G" ]; then
                        new_line="${chrom2}\t${pos2}\t${id2}\t${real_ref2}\tG\t${qual2}\t${filter2}\t${info2}"
                        echo -e $new_line >> $outfile
                fi
                if [ "$real_ref2" != "T" ]; then
                        new_line="${chrom2}\t${pos2}\t${id2}\t${real_ref2}\tT\t${qual2}\t${filter2}\t${info2}"
                        echo -e $new_line >> $outfile
                fi
	done < "$temp_outfile"

	rm $temp_outfile

	echo 'bgzip and tabix of' $outfile
	bgzip -f "${outfile}"
	tabix -p vcf "${outfile}".gz
	mv "${outfile}".gz "${outfile}".bgz
	mv "${outfile}".gz.tbi "${outfile}".bgz.tbi

done < "$list_shards"

# cat list_shards.txt
# shard0001 chrom1 pos1-13077999
# shard0002 chrom1 pos13078000-17150659
# shard0003 chrom1 pos17150660-29953083
# shard0004 chrom1 pos29953084-103888907
# shard0005 chrom1 pos103888908-120722157
# shard0006 chrom1 pos120722158-132010435
# shard0007 chrom1 pos132010436-142756023
# shard0008 chrom1 pos142756024-149484646
# shard0009 chrom1 pos149484647-205997708
# shard0010 chrom1 pos205997709-223772847
# shard0011 chrom1 pos223772848-235217212
# shard0012 chrom1 pos235217213-248983211
# shard0013 chrom1 pos248983212-249250621
# shard0014 chrom2 pos1-5068789
# shard0015 chrom2 pos5068790-16304725
# shard0016 chrom2 pos16304726-21165614
# shard0017 chrom2 pos21165615-33092698
# shard0018 chrom2 pos33092699-33142193
# shard0019 chrom2 pos33142194-87693207
# shard0020 chrom2 pos87693208-93826172
# shard0021 chrom2 pos93826173-110180338
# shard0022 chrom2 pos110180339-149740583
# shard0023 chrom2 pos149740584-234028742
# shard0024 chrom2 pos234028743-243199373
# shard0025 chrom3 pos1-66220271
# shard0026 chrom3 pos66220272-92004855
# shard0027 chrom3 pos92004856-194044607
# shard0028 chrom3 pos194044608-198022430
# shard0029 chrom4 pos1-9299643
# shard0030 chrom4 pos9299644-31829168
# shard0031 chrom4 pos31829169-49413942
# shard0032 chrom4 pos49413943-59764334
# shard0033 chrom4 pos59764335-75439830
# shard0034 chrom4 pos75439831-191154276
# shard0035 chrom5 pos1-17555658
# shard0036 chrom5 pos17555659-47905642
# shard0037 chrom5 pos47905643-91661129
# shard0038 chrom5 pos91661130-138812074
# shard0039 chrom5 pos138812075-155163728
# shard0040 chrom5 pos155163729-180915260
# shard0041 chrom6 pos1-58112660
# shard0042 chrom6 pos58112661-62153590
# shard0043 chrom6 pos62153591-95755544
# shard0044 chrom6 pos95755545-157584468
# shard0045 chrom6 pos157584469-167992074
# shard0046 chrom6 pos167992075-171115067
# shard0047 chrom7 pos1-257485
# shard0048 chrom7 pos257486-50390632
# shard0049 chrom7 pos50390633-59554332
# shard0050 chrom7 pos59554333-74740725
# shard0051 chrom7 pos74740726-100581044
# shard0052 chrom7 pos100581045-130204524
# shard0053 chrom7 pos130204525-139391878
# shard0054 chrom7 pos139391879-154320635
# shard0055 chrom7 pos154320636-159138663
# shard0056 chrom8 pos1-12116855
# shard0057 chrom8 pos12116856-45338888
# shard0058 chrom8 pos45338889-48133050
# shard0059 chrom8 pos48133051-86651452
# shard0060 chrom8 pos86651453-142791516
# shard0061 chrom8 pos142791517-146364022
# shard0062 chrom9 pos1-39688687
# shard0063 chrom9 pos39688688-47110134
# shard0064 chrom9 pos47110135-56392680
# shard0065 chrom9 pos56392681-66429657
# shard0066 chrom9 pos66429658-70785469
# shard0067 chrom9 pos70785470-92393417
# shard0068 chrom9 pos92393418-92603797
# shard0069 chrom9 pos92603798-133148061
# shard0070 chrom9 pos133148062-141213431
# shard0071 chrom10 pos1-17999676
# shard0072 chrom10 pos17999677-38843836
# shard0073 chrom10 pos38843837-49145537
# shard0074 chrom10 pos49145538-51423846
# shard0075 chrom10 pos51423847-125894473
# shard0076 chrom10 pos125894474-135534747
# shard0077 chrom11 pos1-1187760
# shard0078 chrom11 pos1187761-50937354
# shard0079 chrom11 pos50937355-53144206
# shard0080 chrom11 pos53144207-69114802
# shard0081 chrom11 pos69114803-87713379
# shard0082 chrom11 pos87713380-96362585
# shard0083 chrom11 pos96362586-135006516
# shard0084 chrom12 pos1-7214877
# shard0085 chrom12 pos7214878-36356695
# shard0086 chrom12 pos36356696-109398471
# shard0087 chrom12 pos109398472-122555624
# shard0088 chrom12 pos122555625-132756993
# shard0089 chrom12 pos132756994-133851895
# shard0090 chrom13 pos1-86835325
# shard0091 chrom13 pos86835326-112428995
# shard0092 chrom13 pos112428996-115169878
# shard0093 chrom14 pos1-107349540
# shard0094 chrom15 pos1-20914855
# shard0095 chrom15 pos20914856-29184444
# shard0096 chrom15 pos29184445-82854646
# shard0097 chrom15 pos82854647-85009474
# shard0098 chrom15 pos85009475-102531392
# shard0099 chrom16 pos1-8661922
# shard0100 chrom16 pos8661923-34098151
# shard0101 chrom16 pos34098152-40835802
# shard0102 chrom16 pos40835803-88414384
# shard0103 chrom16 pos88414385-90354753
# shard0104 chrom17 pos1-21616609
# shard0105 chrom17 pos21616610-34700849
# shard0106 chrom17 pos34700850-62435761
# shard0107 chrom17 pos62435762-77571462
# shard0108 chrom17 pos77571463-81195210
# shard0109 chrom18 pos1-16960899
# shard0110 chrom18 pos16960900-52134137
# shard0111 chrom18 pos52134138-72308354
# shard0112 chrom18 pos72308355-78077248
# shard0113 chrom19 pos1-8712199
# shard0114 chrom19 pos8712200-20548416
# shard0115 chrom19 pos20548417-26181783
# shard0116 chrom19 pos26181784-59128983
# shard0117 chrom20 pos1-27869570
# shard0118 chrom20 pos27869571-34922086
# shard0119 chrom20 pos34922087-61116438
# shard0120 chrom20 pos61116439-63025520
# shard0121 chrom21 pos1-10059921
# shard0122 chrom21 pos10059922-12763130
# shard0123 chrom21 pos12763131-42980560
# shard0124 chrom21 pos42980561-48129895
# shard0125 chrom22 pos1-16772851
# shard0126 chrom22 pos16772852-20559432
# shard0127 chrom22 pos20559433-50389778
# shard0128 chrom22 pos50389779-51304566
# shard0129 chromX pos1-10763675
# shard0130 chromX pos10763676-37123257
# shard0131 chromX pos37123258-49267998
# shard0132 chromX pos49267999-60132013
# shard0133 chromX pos60132014-76678693
# shard0134 chromX pos76678694-113542669
# shard0135 chromX pos113542670-120038236
# shard0136 chromX pos120038237-143532325
# shard0137 chromX pos143532326-152302100
# shard0138 chromX pos152302101-155270560
# shard0139 chromY pos1-9266323
# shard0140 chromY pos9266324-20168886
# shard0141 chromY pos20168887-23926429
# shard0142 chromY pos23926430-43819362
# shard0143 chromY pos43819363-58942657
# shard0144 chromY pos58942658-59373566
# shard0145 chromMT pos1-16569
