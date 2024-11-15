rule diamond_nr_prokka:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.prokka.bmk"
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.diamond_nr_prokka.bmk",
    log:        dir_out +"/{sample}/diamond/{sample}.diamond_nr_prokka.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.diamond_nr_prokka.bmk"
    params:
        out=    dir_out +"/{sample}/diamond/{sample}.nr_prokka.diamond",
        fa=     dir_out +"/{sample}/prokka/{sample}.faa",
        ref=    diamond_nr,
    threads: 24
    shell:
        '''
        diamond blastp -q {params.fa} -p {threads} -d {params.ref} \
        -e 0.001 --id 40 -f tab -b 50 \
        -o {params.out} &> {log}
        '''

rule diamond_eggnog_prokka:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.prokka.bmk"
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.diamond_eggnog_prokka.bmk",
    log:        dir_out +"/{sample}/diamond/{sample}.diamond_eggnog_prokka.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.diamond_eggnog_prokka.bmk"
    params:
        out=    dir_out +"/{sample}/diamond/{sample}.eggnog_prokka.diamond",
        fa=     dir_out +"/{sample}/prokka/{sample}.faa",
        ref=    diamond_eggnog,
    threads: 24
    shell:
        '''
        diamond blastp -q {params.fa} -p {threads} -d {params.ref} \
        -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
        -o {params.out} &> {log}
        '''

rule diamond_kegg_prokka:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.prokka.bmk"
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.diamond_kegg_prokka.bmk",
    log:        dir_out +"/{sample}/diamond/{sample}.diamond_kegg_prokka.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.diamond_kegg_prokka.bmk"
    params:
        out=    dir_out +"/{sample}/diamond/{sample}.kegg_prokka.diamond",
        fa=     dir_out +"/{sample}/prokka/{sample}.faa",
        ref=    diamond_kegg,
    threads: 24
    shell:
        '''
        diamond blastp -q {params.fa} -p {threads} -d {params.ref} \
        -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
        -o {params.out} &> {log}
        '''

rule hmmsearch:
    input:
        hmm=        "gene_hmm/{hmmId}.hmm",
        bmk=        dir_out +"/{sample}/log/{sample}.prokka.bmk",
    output:
        out_hmm=    dir_out +"/{sample}/hmmsearch/{sample}.{hmmId}.hmmsearch",
    log:            dir_out + "/{sample}/log/{sample}.{hmmId}.hmmsearch.log",
    benchmark:      dir_out + "/{sample}/log/{sample}.{hmmId}.hmmsearch.bmk",
    params:
        fa=         dir_out +"/{sample}/prokka/{sample}.faa",
        all=        dir_out +"/{sample}/hmmsearch/{sample}.{hmmId}.all",
        tblout=     dir_out +"/{sample}/hmmsearch/{sample}.{hmmId}.tblout",
        domtblout=  dir_out +"/{sample}/hmmsearch/{sample}.{hmmId}.domtblout",
        pfamtblout= dir_out +"/{sample}/hmmsearch/{sample}.{hmmId}.pfamtblout"
    threads: 4
    shell:
        '''
        hmmsearch -A {params.all} --tblout {params.tblout} --domtblout {params.domtblout} \
        --pfamtblout {params.pfamtblout} -E 10 -o {output.out_hmm} --cpu {threads} {input.hmm} {params.fa}        
        '''

rule hmmsearch_done:
    input: lambda wildcards: expand(dir_out + "/{sample}/log/{sample}.{hmmId}.hmmsearch.bmk",hmmId=hmmid_list,sample=wildcards.sample)
    output: dir_out + "/{sample}/log/{sample}.hmmsearch_ALLDone.txt"
    shell:
        '''
        touch {output}
        '''

rule kraken2_prokka:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.prokka.bmk"
    output:     dir_out +"/{sample}/log/{sample}.kraken2_prokka.bmk"
    log:        dir_out +"/{sample}/kraken2_prokka/{sample}.kraken2_prokka.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.kraken2_prokka.bmk"
    params:
        kraken2_db= kraken2_db,
        fa=         dir_out +"/{sample}/prokka/{sample}.ffn",
        output_tsv= dir_out +"/{sample}/kraken2_prokka/{sample}.kraken2_prokka.output.tsv",
        report_tsv= dir_out +"/{sample}/kraken2_prokka/{sample}.kraken2_prokka.report.tsv",
    threads: 4
    shell:
        '''
        kraken2 --db {params.kraken2_db} \
        --threads {threads} \
        --output {params.output_tsv} \
        --report {params.report_tsv} \
        --use-names --report-zero-counts  --use-mpa-style  {params.fa} &> {log}
        '''

rule prokka:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.cd_hit_est.bmk",
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.prokka.bmk",
    log:        dir_out +"/{sample}/prokka/{sample}.prokka.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.prokka.bmk"
    params:
        outprefix=  dir_out +"/{sample}/prokka",
        fa= dir_out +"/{sample}/cd_hit_est/{sample}.cd_hit_est.fa",
        id="{sample}"
    threads: 4
    shell:
        '''
        prokka --force --cpus {threads} --outdir {params.outprefix} --prefix {params.id} --metagenome --kingdom Bacteria {params.fa}  &> {log}
        '''

rule kraken2_megahit:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.cd_hit_est.bmk"
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.kraken2_megahit.bmk"
    log:        dir_out +"/{sample}/kraken2_megahit/{sample}.kraken2_megahit.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.kraken2_megahit.bmk"
    params:
        kraken2_db= kraken2_db,
        fa=         dir_out +"/{sample}/cd_hit_est/{sample}.cd_hit_est.fa",
        output_tsv= dir_out +"/{sample}/kraken2_megahit/{sample}.kraken2_megahit.output.tsv",
        report_tsv= dir_out +"/{sample}/kraken2_megahit/{sample}.kraken2_megahit.report.tsv",
    threads: 4
    shell:
        '''
        kraken2 --db {params.kraken2_db} \
        --threads {threads} \
        --output {params.output_tsv} \
        --report {params.report_tsv} \
        --use-names --report-zero-counts  --use-mpa-style  {params.fa} &> {log}
        '''

# rule cd_hit_est:
#     input:      dir_out +"/{sample}/log/{sample}.megahit.bmk",
#     output:     dir_out +"/{sample}/log/{sample}.cd_hit_est.bmk",
#     log:        dir_out +"/{sample}/cd_hit_est/{sample}.cd_hit_est.log"
#     benchmark:  dir_out +"/{sample}/log/{sample}.cd_hit_est.bmk"
#     params:
#         fa_in=  dir_out +"/{sample}/megahit/final.contigs.fa",
#         fa_out= dir_out +"/{sample}/cd_hit_est/{sample}.cd_hit_est.fa",
#     threads: 4
#     shell:
#         '''
#         cd-hit-est -i {params.fa_in} -o {params.fa_out} -T 16 -M 0 -c 0.99 -d 100 -aS 0.9 &> {log}
#         '''

rule quast:
    input:      dir_out +"/{sample}/log/{sample}.megahit.bmk",
    output:     dir_out +"/{sample}/quast/{sample}.quast.bmk",
    log:        dir_out +"/{sample}/quast/{sample}.quast.log"
    benchmark:  dir_out +"/{sample}/quast/{sample}.quast.bmk"
    params:
        fa=     dir_out +"/{sample}/megahit/final.contigs.fa",
        outprefix = dir_out + "/{sample}/quast"
    threads: 4
    shell:
        '''
        quast --threads {threads} {params.fa} -o {params.outprefix} &> {log}
        '''

# rule megahit:
#     input:
#         bmk=    key_smk
#     output:
#         bmk=    dir_out +"/{sample}/log/{sample}.megahit.bmk",
#     log:        dir_out +"/{sample}/log/{sample}.megahit.log"
#     benchmark:  dir_out +"/{sample}/log/{sample}.megahit.bmk"
#     params:
#         outprefix=  dir_out +"/{sample}/megahit",
#         dir_tmp=    dir_out +"/{sample}/megahit_tmp",
#         # fq1=dir_out +"/{sample}/fastqp/{sample}.1.fastq",
#         # fq2=dir_out +"/{sample}/fastqp/{sample}.2.fastq",
#         fq1=dir_fq_in+"/{sample}.1.fastq",
#         fq2=dir_fq_in+"/{sample}.2.fastq"
#     threads: 4
#     shell:
#         '''
#         if [ ! -d "{params.dir_tmp}" ]; then mkdir {params.dir_tmp}; fi
#         if [ -d "{params.outprefix}" ]; then rm -rf {params.outprefix}; fi
#         megahit -t {threads} -m 85899345920 -1 {params.fq1} -2 {params.fq2} -o {params.outprefix} --tmp-dir {params.dir_tmp} &> {log}
#         if [ -d "{params.dir_tmp}" ]; then rm -rf {params.dir_tmp}; fi
#         '''


