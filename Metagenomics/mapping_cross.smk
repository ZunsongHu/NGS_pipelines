rule bowtieBuild_prokka_cross:
    input: dir_out + "/co_assembly/log/prokka_cross.bmk"
    output: dir_out + "/co_assembly/log/bowtieBuild_prokka_cross.bmk"
    log: dir_out + "/co_assembly/bowtieBuild_prokka_cross/bowtieBuild_prokka_cross.log"
    benchmark: dir_out + "/co_assembly/log/bowtieBuild_prokka_cross.bmk"
    params:
        fa=dir_out + "/co_assembly/prokka_cross/co_assembly.ffn",
        base=dir_out + "/co_assembly/bowtieBuild_prokka_cross/prokka"
    threads: 24
    shell:
        '''
        bowtie2-build --threads {threads} {params.fa} {params.base} &> {log}
        '''

rule bowtieMapping_prokka_cross:
    input:
        fq1=    dir_fq+"/{sample}.1.fastq",fq2=dir_fq+"/{sample}.2.fastq",
        bmk=    dir_out +"/co_assembly/log/bowtieBuild_prokka_cross.bmk"
    output:     dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross.bmk"
    benchmark:  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross.bmk"
    log:        dir_out +"/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.bowtieMapping_prokka_cross.log"
    params:
        ref=    dir_out +"/co_assembly/bowtieBuild_prokka_cross/prokka",
        dir=    dir_out +"/{sample}/get_abundance/bowtieMapping_prokka_cross/",
        sam=    dir_out +"/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.sam",
        bam=    dir_out +"/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.bam",
        bai=    dir_out +"/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.bam.bai",
    threads: 4
    shell:
        '''
        bowtie2 -x {params.ref} -1 {input.fq1} -2 {input.fq2} --very-sensitive -S {params.sam} &> {log}
        samtools sort -m 10G -@ {threads} -O bam -T {params.dir} -o {params.bam} {params.sam} &>> {log}
        samtools index {params.bam} {params.bai} &>> {log}
        rm -rf {params.sam}
        '''

rule bowtieMapping_prokka_cross_rmDup:
    input:      dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross.bmk"
    output:     dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross_rmDup.bmk"
    benchmark:  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross_rmDup.bmk"
    log:        dir_out +"/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.bowtieMapping_prokka_cross_rmDup.log"
    params:
        bam_in=     dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.bam",
        bam_out=    dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam",
        bai_out=    dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam.bai",
        metrics=    dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.MarkDuplicates",
        idx = dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam.idx"
    threads: 1
    shell:
        '''
        gatk.426 MarkDuplicates -I {params.bam_in} -O {params.bam_out} --REMOVE_DUPLICATES true -M {params.metrics} &> {log}
        samtools index {params.bam_out} {params.bai_out} &> {log}
        samtools idxstats {params.bam_out} > {params.idx} 2> {log}
        rm {params.bam_in} {params.bam_in}.bai
        '''

# rule bowtieBuild_megahit_cross:
#     input:      dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
#     output:     dir_out +"/co_assembly/log/bowtieBuild_megahit_cross.bmk"
#     log:        dir_out +"/co_assembly/bowtieBuild_megahit_cross/bowtieBuild_megahit_cross.log"
#     benchmark:  dir_out +"/co_assembly/log/bowtieBuild_megahit_cross.bmk"
#     params:
#         fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
#         base=   dir_out +"/co_assembly/bowtieBuild_megahit_cross/megahit"
#     threads: 24
#     shell:
#         '''
#         bowtie2-build --threads {threads} {params.fa} {params.base} &> {log}
#         '''

# rule bowtieMapping_megahit_cross:
#     input:
#         fq1=    dir_fq+"/{sample}.1.fastq",fq2=dir_fq+"/{sample}.2.fastq",
#         bmk=    dir_out +"/co_assembly/log/bowtieBuild_megahit_cross.bmk"
#     output:     dir_out +"/{sample}/log/{sample}.bowtieMapping_megahit_cross.bmk"
#     log:        dir_out +"/{sample}/get_abundance/bowtieMapping_megahit_cross/{sample}.bowtieMapping_megahit_cross.log"
#     benchmark:  dir_out +"/{sample}/log/{sample}.bowtieMapping_megahit_cross.bmk"
#     params:
#         ref=    dir_out +"/co_assembly/bowtieBuild_megahit_cross/megahit",
#         dir=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit_cross/",
#         sam=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit_cross/{sample}.sam",
#         bam=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit_cross/{sample}.bam",
#         bai=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit_cross/{sample}.bam.bai",
#         idx=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit_cross/{sample}.bam.idx",
#     threads: 4
#     shell:
#         '''
#         bowtie2 -x {params.ref} -1 {input.fq1} -2 {input.fq2} --very-sensitive -S {params.sam} &> {log}
#         samtools sort -m 10G -@ {threads} -O bam -T {params.dir} -o {params.bam} {params.sam} &>> {log}
#         samtools index {params.bam} {params.bai} &>> {log}
#         samtools idxstats {params.bam} > {params.idx} 2>> {log}
#         rm {params.sam}
#         '''

# rule bowtieBuild_megahit:
#     input:
#         bmk=    dir_out +"/{sample}/log/{sample}.cd_hit_est.bmk",
#     output:
#         bmk=    dir_out +"/{sample}/log/{sample}.bowtieBuild_megahit.bmk"
#     log:        dir_out +"/{sample}/get_abundance/bowtieBuild_megahit/bowtieBuild_megahit.log",
#     benchmark:  dir_out +"/{sample}/log/{sample}.bowtieBuild_megahit.bmk"
#     params:
#         fa=     dir_out +"/{sample}/cd_hit_est/{sample}.cd_hit_est.fa",
#         dir1=   dir_out +"/{sample}/get_abundance",
#         dir2=   dir_out +"/{sample}/get_abundance/bowtieBuild_megahit",
#         base=   dir_out +"/{sample}/get_abundance/bowtieBuild_megahit/megahit"
#     threads: 4
#     shell:
#         '''
#         bowtie2-build --threads {threads} {params.fa} {params.base} &> {log}
#         '''

# rule bowtieMapping_megahit:
#     input:
#         fq1=    dir_fq_in+"/{sample}.1.fastq",fq2=dir_fq_in+"/{sample}.2.fastq",
#         bmk=    dir_out +"/{sample}/log/{sample}.bowtieBuild_megahit.bmk"
#     output:     dir_out +"/{sample}/log/{sample}.bowtieMapping_megahit.bmk"
#     log:        dir_out +"/{sample}/get_abundance/bowtieMapping_megahit/{sample}.bowtieMapping_megahit.log"
#     benchmark:  dir_out +"/{sample}/log/{sample}.bowtieMapping_megahit.bmk"
#     params:
#         ref=    dir_out +"/{sample}/get_abundance/bowtieBuild_megahit/megahit",
#         dir=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit/",
#         sam=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit/{sample}.sam",
#         bam=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit/{sample}.bam",
#         bai=    dir_out +"/{sample}/get_abundance/bowtieMapping_megahit/{sample}.bam.bai",
#     threads: 4
#     shell:
#         '''
#         bowtie2 -x {params.ref} -1 {input.fq1} -2 {input.fq2} --very-sensitive -S {params.sam} &> {log}
#         samtools sort -m 10G -@ {threads} -O bam -T {params.dir} -o {params.bam} {params.sam} &>> {log}
#         samtools index {params.bam} {params.bai} &>> {log}
#         '''



# samtools faidx {params.fa} &> {log}
# samtools import {params.faidx} {params.sam} {output.bam_unsort} &> {log}
# samtools sort -m 10G -@ {threads} -O bam -T {params.dir} -o {params.bam} {output.bam_unsort} &> {log}
# samtools index {params.bam} {params.bai} &>> {log}

# rule bowtieBuild_prokka:
#     input:
#         bmk=dir_out + "/{sample}/log/{sample}.prokka.bmk",
#     output:
#         bmk=dir_out + "/{sample}/log/{sample}.bowtieBuild_prokka.bmk"
#     log: dir_out + "/{sample}/log/{sample}.bowtieBuild_prokka.log"
#     benchmark: dir_out + "/{sample}/log/{sample}.bowtieBuild_prokka.bmk"
#     params:
#         fa=dir_out + "/{sample}/prokka/{sample}.ffn",
#         dir1=dir_out + "/{sample}/get_abundance",
#         dir2=dir_out + "/{sample}/get_abundance/bowtieBuild",
#         base=dir_out + "/{sample}/get_abundance/bowtieBuild/prokka"
#     threads: 4
#     shell:
#         '''
#         if [ ! -d "{params.dir1}" ]; then mkdir {params.dir1}; fi
#         if [ ! -d "{params.dir2}" ]; then mkdir {params.dir2}; fi
#         bowtie2-build --threads {threads} {params.fa} {params.base} &> {log}
#         '''

# rule bowtieMapping_prokka:
#     input:
#         fq1=    dir_fq_in+"/{sample}.1.fastq",fq2=dir_fq_in+"/{sample}.2.fastq",
#         bmk=    dir_out +"/{sample}/log/{sample}.bowtieBuild_prokka.bmk"
#     output:
#         sam=    dir_out +"/{sample}/get_abundance/bowtieMapping/{sample}.sam"
#     log:        dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka.log"
#     benchmark:  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka.bmk"
#     params:
#         ref=    dir_out +"/{sample}/get_abundance/bowtieBuild/prokka",
#     threads: 4
#     shell:
#         '''
#         bowtie2 -x {params.ref} -1 {input.fq1} -2 {input.fq2} --very-sensitive -S {output.sam} &> {log}
#         samtools sort -m 10G -@ {threads} -O bam -T {params.dir} -o {params.bam} {params.sam} &>> {log}
#         samtools index {params.bam} {params.bai} &>> {log}
#         samtools idxstats {params.bam} > {params.idx} 2>> {log}
#         rm -rf {params.sam}

#         '''
