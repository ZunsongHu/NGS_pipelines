rule getAbundance_BowtieProkka_cross:
    input:
        bmk =   dir_out +"/{id}/log/{id}.convertBowtieSam2Bam_cross.bmk",
    output:
        bmk =   dir_out +"/{id}/log/{id}.getAbundance_BowtieProkka_cross.bmk",
    log:        dir_out +"/{id}/get_abundance/getAbundance_BowtieProkka_cross/{id}.getAbundance_BowtieProkka_cross.log"
    benchmark:  dir_out +"/{id}/log/{id}.getAbundance_BowtieProkka_cross.bmk",
    params:
        bam=    dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/{id}.bam",
        idx=    dir_out +"/{id}/get_abundance/getAbundance_BowtieProkka_cross/{id}.bam.idx",
    shell:
        '''
        samtools idxstats {params.bam} > {params.idx} &> {log}
        '''



#module load SAMtools/1.9-foss-2018b
rule convertBowtieSam2Bam_cross:
    input:
        bmk=    dir_out +"/{id}/log/{id}.bowtieMapping_prokka_cross.bmk"
    output:
        bmk=    dir_out +"/{id}/log/{id}.convertBowtieSam2Bam_cross.bmk",
        bam_unsort=temp(dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/{id}.unsort.bam")
    log:        dir_out +"/{id}/log/{id}.convertBowtieSam2Bam_cross.log"
    benchmark:  dir_out +"/{id}/log/{id}.convertBowtieSam2Bam_cross.bmk"
    params:
        fa=     dir_out +"/co_assembly/prokka_cross/co_assembly.ffn",
        faidx=  dir_out +"/co_assembly/prokka_cross/co_assembly.ffn.fai",
        sam=    dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/{id}.sam",
        dir=    dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/",
        bam=    dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/{id}.bam",
        bai=    dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/{id}.bam.bai",
    threads: 4
    shell:
        '''
        samtools faidx {params.fa} &> {log}
        samtools import {params.faidx} {params.sam} {output.bam_unsort} &> {log}        
        samtools sort -m 10G -@ {threads} -O bam -T {params.dir} -o {params.bam} {output.bam_unsort} &> {log}
        samtools index {params.bam} {params.bai} &>> {log}
        '''



rule bowtieMapping_prokka_cross:
    input:
        fq1=    dir_fq+"/{id}.1.fastq",fq2=dir_fq+"/{id}.2.fastq",
        bmk=    dir_out +"/co_assembly/log/bowtieBuild_prokka_cross.bmk"
    output:
        sam=    dir_out +"/{id}/get_abundance/bowtieMapping_prokka_cross/{id}.sam"
    log:        dir_out +"/{id}/log/{id}.bowtieMapping_prokka_cross.log"
    benchmark:  dir_out +"/{id}/log/{id}.bowtieMapping_prokka_cross.bmk"
    params:
        ref=    dir_out +"/co_assembly/bowtieBuild/prokka"
    threads: 4
    shell:
        '''
        bowtie2 -x {params.ref} -1 {input.fq1} -2 {input.fq2} --very-sensitive -S {output.sam} &> {log}
        '''

rule kraken2_megahit_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/megahit_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/kraken2_megahit_cross.bmk"
    log:        dir_out +"/co_assembly/kraken2_megahit_cross/kraken2_megahit_cross.log"
    benchmark:  dir_out +"/co_assembly/log/kraken2_megahit_cross.bmk"
    params:
        fa=         dir_out +"/co_assembly/megahit_cross/final.contigs.fa",
        kraken2_db= kraken2_db,
        output_tsv= dir_out +"/co_assembly/kraken2_megahit_cross/kraken2_megahit_cross.output.tsv",
        report_tsv= dir_out +"/co_assembly/kraken2_megahit_cross/kraken2_megahit_cross.report.tsv",
    threads: 64
    shell:
        '''
        kraken2 --db {params.kraken2_db} \
        --threads {threads} \
        --output {params.output_tsv} \
        --report {params.report_tsv} \
        --use-names --use-mpa-style\
        {params.fa} &> {log}
        '''

rule diamond_nr_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/diamond_nr_cross.bmk"
    log:        dir_out +"/co_assembly/diamond_nr_cross/diamond_nr_cross.log",
    benchmark:  dir_out +"/co_assembly/log/diamond_nr_cross.bmk"
    params:
        out=    dir_out +"/co_assembly/diamond_nr_cross/diamond_nr_cross.diamond",
        fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
        ref=    diamond_nr,
        diamond=diamond
    threads: 24
    shell:
        '''
        {params.diamond} blastp -q {params.fa} -p {threads} -d {params.ref} \
        -e 0.001 --id 40 -f tab -b 50 \
        -o {params.out} &> {log}
        '''

rule diamond_eggnog_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/diamond_eggnog_cross.bmk"
    log:        dir_out +"/co_assembly/diamond_eggnog_cross/diamond_eggnog_cross.log",
    benchmark:  dir_out +"/co_assembly/log/diamond_eggnog_cross.bmk"
    params:
        out=    dir_out +"/co_assembly/diamond_eggnog_cross/diamond_eggnog_cross.diamond",
        fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
        ref=    diamond_eggnog,
        diamond=diamond
    threads: 24
    shell:
        '''
        {params.diamond} blastp -q {params.fa} -p {threads} -d {params.ref} \
        -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
        -o {params.out} &> {log}
        '''

rule diamond_kegg_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/diamond_kegg_cross.bmk"
    log:        dir_out +"/co_assembly/diamond_kegg_cross/diamond_kegg_cross.log",
    benchmark:  dir_out +"/co_assembly/log/diamond_kegg_cross.bmk"
    params:
        out=    dir_out +"/co_assembly/diamond_kegg_cross/diamond_kegg_cross.diamond",
        fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
        ref=    diamond_kegg,
        diamond=diamond
    threads: 24
    shell:
        '''
        {params.diamond} blastp -q {params.fa} -p {threads} -d {params.ref} \
        -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
        -o {params.out} &> {log}
        '''

rule hmmsearch_cross:
    input:
        hmm=        dir_out +"/co_assembly/log/prokka_cross.bmk"
    output:
        out_hmm=    dir_out +"/co_assembly/log/prokka_cross.{hmmId}.bmk"
    log:            dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.log"
    benchmark:      dir_out +"/co_assembly/log/prokka_cross.{hmmId}.bmk"
    params:
        fa=         dir_out +"/co_assembly/prokka_cross/co_assembly.faa",
        out_hmm=    dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.hmmsearch",
        all=        dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.all",
        tblout=     dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.tblout",
        domtblout=  dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.domtblout",
        pfamtblout= dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.pfamtblout"
    threads: 4
    shell:
        '''
        hmmsearch -A {params.all} --tblout {params.tblout} --domtblout {params.domtblout} \
        --pfamtblout {params.pfamtblout} -E 10 -o {params.out_hmm} --cpu {threads} {input.hmm} {params.fa}        
        '''

rule bowtieBuild_prokka_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/bowtieBuild_prokka_cross.bmk"
    log:        dir_out +"/co_assembly/bowtieBuild/bowtieBuild_prokka_cross.log"
    benchmark:  dir_out +"/co_assembly/log/bowtieBuild_prokka_cross.bmk"
    params:
        fa=     dir_out +"/co_assembly/prokka_cross/co_assembly.ffn",
        base=   dir_out +"/co_assembly/bowtieBuild/prokka"
    threads: 24
    shell:
        '''
        bowtie2-build --threads {threads} {params.fa} {params.base} &> {log}
        '''

rule prokka_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
    log:        dir_out +"/co_assembly/log/prokka_cross.log"
    benchmark:  dir_out +"/co_assembly/log/prokka_cross.bmk"
    params:
        outprefix=  dir_out +"/co_assembly/prokka_cross",
        fa= dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
        id="co_assembly",
        prokka=prokka
    threads: 24
    shell:
        '''
        {params.prokka} --force --cpus {threads} --outdir {params.outprefix} --prefix {params.id} --metagenome --kingdom Bacteria {params.fa}  &> {log}
        '''


rule cd_hit_est_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/megahit_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    log:        dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.log"
    benchmark:  dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    params:
        fa_in=  dir_out +"/co_assembly/megahit_cross/final.contigs.fa",
        fa_out= dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa"
    threads: 24
    shell:
        '''
        cd-hit-est -i {params.fa_in} -o {params.fa_out} -T {threads} -M 0 -c 0.99 -d 100 -aS 0.9 &> {log}
        '''

rule quast_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/megahit_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/quast_cross.bmk"
    log:        dir_out +"/co_assembly/log/quast_cross.log"
    benchmark:  dir_out +"/co_assembly/log/quast_cross.bmk"
    params:
        fa=     dir_out +"/co_assembly/megahit_cross/final.contigs.fa",
        outprefix =dir_out +"/co_assembly/quast_cross"
    threads: 24
    shell:
        '''
        quast --threads {threads} {params.fa} -o {params.outprefix} &> {log}
        '''

rule megahit_cross:
    input:
        expand(dir_out +"/{id}/log/{id}.fastqp.bmk",id=idlist)
    output:
        bmk=    dir_out +"/co_assembly/log/megahit_cross.bmk",
    log:        dir_out +"/co_assembly/log/megahit_cross.log"
    benchmark:  dir_out +"/co_assembly/log/megahit_cross.bmk"
    params:
        fq1=    fq1_all,
        fq2=    fq2_all,
        outprefix=  dir_out +"/co_assembly/megahit_cross",
        dir_tmp=    dir_out +"/co_assembly/megahit_cross_tmp"
    threads: 24
    shell:
        '''
        if [ ! -d "{params.dir_tmp}" ]; then mkdir {params.dir_tmp}; fi
        if [ -d "{params.outprefix}" ]; then rm -rf {params.outprefix}; fi
        megahit -t {threads} -1 {params.fq1} -2 {params.fq2} -o {params.outprefix} --tmp-dir {params.dir_tmp} &> {log}
        if [ -d "{params.dir_tmp}" ]; then rm -rf {params.dir_tmp}; fi
        '''
