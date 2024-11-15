

rule centrifuge_all:
    input:
        i1 = dir_in+"/{sample}_R1.fastq.gz",
        i2 = dir_in+"/{sample}_R2.fastq.gz"
    output:
        classification = dir_out+"/1_centrifuge_all/{sample}.classification",
        report = dir_out+"/1_centrifuge_all/{sample}.report"
    log:        dir_out+"/logs/1_centrifuge_all/{sample}.log"
    benchmark:  dir_out+"/benchmarks/1_centrifuge_all/{sample}.bmk"
    params:
        centrifuge=centrifuge,
        ref_=ref_centrifuge_all
    threads: 8
    shell:
        '''
        {params.centrifuge} -x {params.ref_} -q -t -p {threads} -1 {input.i1} -2 {input.i2} \
        -S {output.classification} --report-file {output.report} &> {log}
        '''

rule write_code:
    input:
        key_smk=key_smk
    output:
        code=   dir_code +"/{sample}.sh"
    params:
        cores=      cores,
        mem=        mem,
        walltime=   walltime,
        job_name=   job_prefix +"_{sample}",
        dir_log=    dir_log,
        log=        dir_log +"/J%j-{sample}.log",
        err=        dir_log +"/J%j-{sample}.err",
        key_file=   key_file
    shell:
        '''
        if [ ! -d "{params.dir_log}" ]; then mkdir {params.dir_log}; fi

        echo "#!/bin/bash"                                                          > {output.code}
        echo "#SBATCH --job-name={params.job_name}"                                 >> {output.code}
        echo "#SBATCH -n {params.cores}"                                            >> {output.code}
        echo "#SBATCH -N 1-1"                                                       >> {output.code}
        echo "#SBATCH -p all"                                                      >> {output.code}
        echo "#SBATCH --mem={params.mem}"                                           >> {output.code}
        echo "#SBATCH --time={params.walltime}"                                     >> {output.code}
        echo "#SBATCH --output={params.log}"                                        >> {output.code}
        echo "#SBATCH --error={params.err}"                                         >> {output.code}

        echo "snakemake -s {input.key_smk} -p {params.key_file}  -j{params.cores}"  >> {output.code}
        '''


# rule centrifuge_bacteria:
#     input:
#         i1 = dir_in+"/{sample}_R1.fastq",
#         i2 = dir_in+"/{sample}_R2.fastq"
#     output:
#         classification = dir_out+"/1_centrifuge/{sample}.classification",
#         report = dir_out+"/1_centrifuge/{sample}.report"
#     log:        dir_out+"/logs/1_centrifuge/{sample}.log"
#     benchmark:  dir_out+"/benchmarks/1_centrifuge/{sample}.tsv"
#     params:
#         centrifuge=centrifuge,
#         ref_bacteria=ref_bacteria
#     threads: 4
#     shell:
#         '''
#         {params.centrifuge} -x {params.ref_bacteria} -q -t -p {threads} -1 {input.i1} -2 {input.i2} \
#         -S {output.classification} --report-file {output.report} &> {log}
#         '''






# rule unmapped_reads_from_ctf:
#     input:
#         i1 = dir_in+"/{sample}_R1.fastq",
#         i2 = dir_in+"/{sample}_R2.fastq",
#         ctf = dir_out+"/1_centrifuge/{sample}.classification",
#     output:
#         o1 = temp(dir_out+"/2_centrifuge_unmapped/{sample}.unmapped.R1.fastq"),
#         o2 = temp(dir_out+"/2_centrifuge_unmapped/{sample}.unmapped.R2.fastq"),
#     log:        dir_out+"/logs/2_centrifuge_unmapped/{sample}.log"
#     benchmark:  dir_out+"/benchmarks/2_centrifuge_unmapped/{sample}.tsv"
#     params:
#         URC=URC
#     shell:
#         "{params.URC} {input.ctf} {input.i1} {input.i2} {output.o1} {output.o2} &> {log}"
#
# rule cat_unmapped_reads_1:
#     input:
#         i1 = expand(dir_out+"/2_centrifuge_unmapped/{sample}.unmapped.R1.fastq",sample=samplelist)
#     output:
#         o1 = dir_out+"/3_cross_assembly/1_unmapped.fastq",
#         unmapped_reads_number=dir_out+"/docs/wc_unmapped"
#     log:        dir_out+"/logs/3_cross_assembly/cat_1.log"
#     benchmark:  dir_out+"/benchmarks/3_cross_assembly/cat_1.tsv"
#     shell:
#         '''
#         cat {input.i1} > {output.o1} 2> {log}
#         wc -l {output.o1} > {output.unmapped_reads_number} 2> {log}
#         '''
#
# rule cat_unmapped_reads_2:
#     input:
#         i2 = expand(dir_out+"/2_centrifuge_unmapped/{sample}.unmapped.R2.fastq",sample=samplelist)
#     output:
#         o2 = dir_out+"/3_cross_assembly/2_unmapped.fastq"
#     log:        dir_out+"/logs/3_cross_assembly/cat_2.log"
#     benchmark:  dir_out+"/benchmarks/3_cross_assembly/cat_2.tsv"
#     shell:
#         "cat {input.i2} > {output.o2} 2> {log}"
#
# rule megahit_assembly:
#     input:
#         i1 = dir_out+"/3_cross_assembly/1_unmapped.fastq",
#         i2 = dir_out+"/3_cross_assembly/2_unmapped.fastq"
#     output:
#         dir=directory(dir_out+"/3_cross_assembly/megahit_out"),
#         #final_contigs_fa= dir_out + "/3_cross_assembly/megahit_out/final.contigs.fa",
#     log:        dir_out+"/logs/3_cross_assembly/assembly.log"
#     benchmark:  dir_out+"/benchmarks/3_cross_assembly/assembly.tsv"
#     params:
#         MIN_CONTIG_LENGTH=MIN_CONTIG_LENGTH
#     threads: 32
#     shell:
#         '''
#         megahit -t {threads} -1 {input.i1} -2 {input.i2} -o {output.dir} --min-contig-len  {params.MIN_CONTIG_LENGTH} &> {log}
#         '''
#
# rule kraken:
#     input:
#         dir=dir_out+"/3_cross_assembly/megahit_out",
#     output:
#         labels= dir_out+"/4_kraken/final.contigs.labels",
#         report= dir_out+"/4_kraken/final.contigs.report"
#     log:        dir_out+"/logs/3_kraken/kraken.log"
#     benchmark:  dir_out+"/benchmarks/3_kraken/kraken.tsv"
#     params:
#         ref=ref_kraken,
#         dir_kraken=dir_kraken,
#         final_contigs_fa = dir_out+"/3_cross_assembly/megahit_out/final.contigs.fa",
#         final_contigs_kraken=dir_out+"/4_kraken/final.contigs.kraken",
#     threads: 8
#     shell:
#         '''
#         {params.dir_kraken}/kraken --threads {threads} --db {params.ref} --fasta-input {params.final_contigs_fa} --output  {params.final_contigs_kraken}
#         {params.dir_kraken}/kraken-translate --mpa-format --db {params.ref} {params.final_contigs_kraken} > {output.labels}
#         {params.dir_kraken}/kraken-mpa-report --db {params.ref} {params.final_contigs_kraken} > {output.report}
#         '''

# rule bowtie2_build:
#     input:
#         final_contigs_fa=dir_out+"/3_cross_assembly/megahit_out/final.contigs.fa",
#     output:
#         bmk=dir_out+"/logs/bowtie2_build.bmk"
#     log:        dir_out+"/logs/bowtie2_build.log"
#     benchmark:  dir_out+"/logs/bowtie2_build.bmk"
#     params:
#         basename=dir_out+"/bowtie2_build/final.contigs.fa",
#         dir=dir_out+"/bowtie2_build"
#     threads: 1
#     shell:
#         '''
#         if [ ! -d {params.dir} ];then mkdir {params.dir}; fi
#         bowtie2-build {input.final_contigs_fa} {params.basename} &> {log}
#         '''

# rule bowtie2_align:
#     input:
#         i1 = dir_in+"/{sample}_R1.fastq",
#         i2 = dir_in+"/{sample}_R2.fastq",
#         #bmk = dir_out + "/logs/bowtie2_build.bmk"
#     output:
#         sam1=temp(dir_out+"/bowtie2_align/{sample}/{sample}.sam"),
#         bam1=temp(dir_out+"/bowtie2_align/{sample}/{sample}.bam"),
#         bam2=dir_out+"/bowtie2_align/{sample}/{sample}.sort.bam",
#         idx=dir_out+"/bowtie2_align/{sample}/{sample}.idx"
#     log:        dir_out+"/logs/bowtie2_align/{sample}.log"
#     benchmark:  dir_out+"/logs/bowtie2_align/{sample}.bmk"
#     params:
#         ref_basename=dir_out+"/bowtie2_build/final.contigs.fa",
#     threads: 4
#     shell:
#         '''
#         bowtie2 -p {threads} -q -x {params.ref_basename} -1 {input.i1} -2 {input.i2} -S {output.sam1} &> {log}
#         samtools view -b -S {output.sam1} > {output.bam1}
#         samtools sort {output.bam1} -o {output.bam2}
#         samtools index {output.bam2}
#         samtools idxstats {output.bam2} > {output.idx}
#         '''

# rule calculate_abundance:
#     input:
#         dir_out + "/logs/bowtie2_align/{sample}.bmk"
#     output:
#         "2.data_temp/abundance_out_EachFile/s_abundance_{sample}.csv"
#     log:        dir_out+"/logs/calculate_abundance/{sample}.log"
#     benchmark:  dir_out+"/logs/calculate_abundance/{sample}.bmk"
#     params:
#         sample="{sample}"
#     threads: 1
#     shell:
#         '''
#         Rscript id/R_code/3.calculate_abundance.R {params.sample}
#         '''





# Exclude bacteria
# bowtie2-build ../ViralPip/3_cross_assembly/megahit_out/final.contigs.fa final.contigs.fa_bowtie2Index
# bowtie2 -q -x final.contigs.fa_bowtie2Index -1 ../../data_fastq/Shen-MGS12002_1.fastq -2 ../../data_fastq/Shen-MGS12002_2.fastq -S sam/Shen-MGS12002.sam
#
# /lustre/project/hdeng2/jgreenb8/bowtie2/bowtie2 -q -x ./final.contigs.fa -1 ADD05009_Stool.rmhost.1.fq -2 ADD05009_Stool.rmhost.2.fq -S SAM/ADD05009.sam
#
#
#
# samtools view -b -S SAM/ADD05009.sam > SAM/ADD05009.sam.bam
# samtools sort SAM/ADD05009.sam.bam -o SAM/ADD05009.sam.bam.sorted.bam
# samtools view -h SAM/ADD05009.sam.bam.sorted.bam > SAM/ADD05009.sam.sorted.sam
# samtools index SAM/ADD05009.sam.bam.sorted.bam
# samtools idxstats SAM/ADD05009.sam.bam.sorted.bam > SAM/ADD05009.sam.bam.sorted.bam.idxstats
#
# rm SAM/ADD05009.sam.bam
# rm SAM/ADD05009.sam.bam.sorted.bam
# rm SAM/ADD05009.sam.sorted.sam
