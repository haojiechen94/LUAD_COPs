# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn


# preprocess for WES data
# step1: QC(trim galore)
# step2: mapping(bwa)
# step3: sort(samtools) and deduplicate(picard)
# step4: QC for bases(BQSR)
# step5: call CNV(cnvkit)


# take P12 as an example
meta_info=/data/WES/raw_data/file_info.xlsx
file_info=(`echo $meta_info |cut -d ',' -f 1,2,3,4,5 |tr ',' ' '`)
name="${file_info[0]}"
tumor_read1_fq_gz="$input_dir/${file_info[1]}"
tumor_read2_fq_gz="$input_dir/${file_info[2]}"
normal_read1_fq_gz="$input_dir/${file_info[3]}"
normal_read2_fq_gz="$input_dir/${file_info[4]}"

input_dir=/data/WES/raw_data/
output_dir=/data/WES/output/

genome_version='hg38'


# +++++++++++++++++++++++++++++++++
# step1: QC (trim galore)
trim_galore --paired -o "${output_dir}/step1_trim_galore_cutting_adapters/" ${normal_read1_fq_gz} ${normal_read2_fq_gz};
trim_galore --paired -o "${output_dir}/step1_trim_galore_cutting_adapters/" ${tumor_read1_fq_gz} ${tumor_read2_fq_gz};



# +++++++++++++++++++++++++++++++++
# step2: mapping (bwa)
reference_genome=hg38

tumor_read1_name=(`echo $input_dir/${file_info[1]}| tr '.' ' '`)
tumor_read2_name=(`echo $$input_dir/${file_info[2]}| tr '.' ' '`)
normal_read1_name=(`echo $input_dir/${file_info[3]}| tr '.' ' '`)
normal_read2_name=(`echo $$input_dir/${file_info[4]}| tr '.' ' '`)

tumor_read1_fq_gz="$output_dir/step1_trim_galore_cutting_adapters/${tumor_read1_name[0]}""_val_1.fq.gz"
tumor_read2_fq_gz="$output_dir/step1_trim_galore_cutting_adapters/${tumor_read2_name[0]}""_val_2.fq.gz"
normal_read1_fq_gz="$output_dir/step1_trim_galore_cutting_adapters/${normal_read1_name[0]}""_val_1.fq.gz"
normal_read2_fq_gz="$output_dir/step1_trim_galore_cutting_adapters/${normal_read2_name[0]}""_val_2.fq.gz"

bwa mem -t 10 -M -Y -R "@RG\tID:${name}tumorWES\tSM:${name}tumor\tLB:WES\tPL:Illumina" ${reference_genome} $tumor_read1_fq_gz $tumor_read2_fq_gz | samtools view -Sb - > "${output_dir}/step2_bwa_mapping/${name}tumor"".bam"
bwa mem -t 10 -M -Y -R "@RG\tID:${name}normalWES\tSM:${name}normal\tLB:WES\tPL:Illumina" ${reference_genome} $normal_read1_fq_gz $normal_read2_fq_gz | samtools view -Sb - > "${output_dir}/step2_bwa_mapping/${name}normal"".bam"





#++++++++++++++++++++++++++++++++++++
# step3: sort(samtools) and deduplicate(picard)
tumor_bam="${output_dir}/step2_bwa_mapping/${name}tumor"".bam"
normal_bam="${output_dir}/step2_bwa_mapping/${name}normal"".bam"

tumor_sorted_bam="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}tumor"".sorted.bam"
normal_sorted_bam="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}normal"".sorted.bam"

tumor_sorted_markdup_bam="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}tumor"".sorted.markdup.bam"
normal_sorted_markdup_bam="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}normal"".sorted.markdup.bam"


samtools sort -@ 5 "${tumor_bam}" -o "${tumor_sorted_bam}";
samtools sort -@ 5 "${normal_bam}" -o "${normal_sorted_bam}";

picard MarkDuplicates -Xmx64g I="${tumor_sorted_bam}" O="${tumor_sorted_markdup_bam}" M="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}tumor"".sorted.markdup.txt" REMOVE_DUPLICATES=true;
picard MarkDuplicates -Xmx64g I="${normal_sorted_bam}" O="${normal_sorted_markdup_bam}" M="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}norml"".sorted.markdup.txt" REMOVE_DUPLICATES=true;

picard BuildBamIndex -Xmx64g I="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${tumor_sorted_markdup_bam}";
picard BuildBamIndex -Xmx64g I="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${normal_sorted_markdup_bam}";






# ++++++++++++++++++++++++++++++++++++++++++++
# step4: QC for bases(BQSR)
genome_fa="/data/WES/${genome_version}_genome_fa/${genome_version}.fa"
if [ "$genome_version" = "hg38" ];
    then
        dbsnp=/data/WES/database/Homo_sapiens_assembly38.dbsnp138.vcf
        mills=/data/WES/database/Mills_and_1000G_gold_standard.indels.hg38.vcf
        G1000=/data/WES/database/1000G_phase1.snps.high_confidence.hg38.vcf
fi

tumor_sorted_markdup_bam="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}tumor"".sorted.markdup.bam"
normal_sorted_markdup_bam="${output_dir}/step3_sorting_alignments_and_removing_duplicates/${name}normal"".sorted.markdup.bam"

gatk BaseRecalibrator -R ${genome_fa} \
-I "${tumor_sorted_markdup_bam}" \
--known-sites ${dbsnp} --known-sites ${mills} --known-sites ${G1000} \
-O "${output_dir}/step4_applying_base_quality_score_recalibration/${name}""tumor.BQSR.table";

gatk ApplyBQSR --bqsr-recal-file "${output_dir}/step4_applying_base_quality_score_recalibration/${name}""tumor.BQSR.table" \
-R ${genome_fa} \
-I "${tumor_sorted_markdup_bam}" \
-O "${output_dir}/step4_applying_base_quality_score_recalibration/${name}""tumor.sorted.markdup.BQSR.bam";

gatk BaseRecalibrator -R ${genome_fa} \
-I "${normal_sorted_markdup_bam}" \
--known-sites ${dbsnp} --known-sites ${mills} --known-sites ${G1000} \
-O "${output_dir}/step4_applying_base_quality_score_recalibration/${name}""normal.BQSR.table";

gatk ApplyBQSR --bqsr-recal-file "${output_dir}/step4_applying_base_quality_score_recalibration/${name}""normal.BQSR.table" \
-R ${genome_fa} \
-I "${normal_sorted_markdup_bam}" \
-O "${output_dir}/step4_applying_base_quality_score_recalibration/${name}""normal.sorted.markdup.BQSR.bam";







# ++++++++++++++++++++++++++++++++++++
# call CNV(cnvkit)

genome_fa="/data/WES/${genome_version}_genome_fa/${genome_version}.fa"
access_bed="/data/WES/""access-5kb.${genome_version}.bed"
refflat="/data/WES""${genome_version}_refFlat.txt"

guess_baits.py -a ${access_bed} "${output_dir}/step4_applying_base_quality_score_recalibration/*"".sorted.markdup.BQSR.bam" -o "${outpur_dir}/step5_calling_somatic_mutation/""cnvkit_baits.bed"

cnvkit.py batch "${output_dir}/step4_applying_base_quality_score_recalibration/*""tumor.sorted.markdup.BQSR.bam" \
--normal "${output_dir}/step4_applying_base_quality_score_recalibration/*""normal.sorted.markdup.BQSR.bam" \
--targets "${outpur_dir}/step5_calling_somatic_mutation/""cnvkit_baits.bed" \
--annotate ${refflat} \
--fasta ${genome_fa} \
--access ${access_bed} \
--output-reference "${outpur_dir}/step5_calling_somatic_mutation/""cnvkit_reference.cnn" \
--output-dir "${outpur_dir}/step5_calling_somatic_mutation/""cnvkit_results/" \
--diagram --scatter
