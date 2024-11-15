rule diamond_build:
    input:
        fa= fa_for_diamond_build
    output:         dir_ref+"diamond/log/{diamond_db}.bmk"
    benchmark:      dir_ref+"diamond/log/{diamond_db}.bmk"
    log:            dir_ref+"diamond/log/{diamond_db}.log"
    params: dmnd=   dir_ref+"diamond/{diamond_db}.dmnd"
    threads: 32
    shell:
        '''
        diamond makedb --threads {threads} --in {input.fa} --db {params.dmnd} &> {log}
        '''

rule diamond_prokka_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/diamond_prokka_cross.{diamond_db}.bmk"
    log:        dir_out +"/co_assembly/diamond_prokka_cross/{diamond_db}/diamond_prokka_cross.log",
    benchmark:  dir_out +"/co_assembly/log/diamond_prokka_cross.{diamond_db}.bmk"
    params:
        out=    dir_out +"/co_assembly/diamond_prokka_cross/{diamond_db}/diamond_prokka_cross.{diamond_db}.txt",
        fa=     dir_out +"/co_assembly/prokka_cross/co_assembly.faa",
        ref=    lambda wildcards: dict_diamond_ref_db[wildcards.diamond_db],
        para=   lambda wildcards: dict_diamond_ref_db_para[wildcards.diamond_db]
    threads: 24
    shell:
        '''
        diamond blastp --threads {threads} --db {params.ref} \
        --query {params.fa} \
        {params.para} \
        --evalue 0.001 --outfmt 6 -o {params.out} &> {log}
        '''

# rule diamond_megahit_cross:
#     input: bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"

# rule diamond_nr_cross_prokka:
#     input:
#         bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
#     output:
#         bmk=    dir_out +"/co_assembly/log/diamond_nr_cross_prokka.bmk"
#     log:        dir_out +"/co_assembly/diamond_nr_cross_prokka/diamond_nr_cross_prokka.log",
#     benchmark:  dir_out +"/co_assembly/log/diamond_nr_cross_prokka.bmk"
#     params:
#         out=    dir_out +"/co_assembly/diamond_nr_cross_prokka/diamond_nr_cross_prokka",
#         fa=     dir_out +"/co_assembly/prokka_cross/co_assembly.faa",
#         ref=    diamond_nr,
#     threads: 24
#     shell:
#         '''
#         diamond blastp -q {params.fa} -p {threads} -d {params.ref} \
#         -e 0.001 --id 40 -f tab -b 50 \
#         -o {params.out} &> {log}
#         '''
#
# rule diamond_eggnog_cross_prokka:
#     input:
#         bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
#     output:
#         bmk=    dir_out +"/co_assembly/log/diamond_eggnog_cross_prokka.bmk"
#     log:        dir_out +"/co_assembly/diamond_eggnog_cross_prokka/diamond_eggnog_cross_prokka.log",
#     benchmark:  dir_out +"/co_assembly/log/diamond_eggnog_cross_prokka.bmk"
#     params:
#         out=    dir_out +"/co_assembly/diamond_eggnog_cross_prokka/diamond_eggnog_cross_prokka",
#         fa=     dir_out +"/co_assembly/prokka_cross/co_assembly.faa",
#         ref=    diamond_eggnog,
#     threads: 24
#     shell:
#         '''
#         diamond blastp -q {params.fa} -p {threads} -d {params.ref} \
#         -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
#         -o {params.out} &> {log}
#         '''
#
# rule diamond_kegg_cross_prokka:
#     input:
#         bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
#     output:
#         bmk=    dir_out +"/co_assembly/log/diamond_kegg_cross_prokka.bmk"
#     log:        dir_out +"/co_assembly/diamond_kegg_cross_prokka/diamond_kegg_cross_prokka.log",
#     benchmark:  dir_out +"/co_assembly/log/diamond_kegg_cross_prokka.bmk"
#     params:
#         out=    dir_out +"/co_assembly/diamond_kegg_cross_prokka/diamond_kegg_cross_prokka",
#         fa=     dir_out +"/co_assembly/prokka_cross/co_assembly.faa",
#         ref=    diamond_kegg,
#     threads: 24
#     shell:
#         '''
#         diamond blastp -q {params.fa} -p {threads} -d {params.ref} \
#         -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
#         -o {params.out} &> {log}
#         '''
#
#
# rule diamond_nr_cross_megahit:
#     input:
#         bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
#     output:
#         bmk=    dir_out +"/co_assembly/log/diamond_nr_cross_megahit.bmk"
#     log:        dir_out +"/co_assembly/diamond_nr_cross_megahit/diamond_nr_cross_megahit.log",
#     benchmark:  dir_out +"/co_assembly/log/diamond_nr_cross_megahit.bmk"
#     params:
#         out=    dir_out +"/co_assembly/diamond_nr_cross_megahit/diamond_nr_cross_megahit",
#         fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
#         ref=    diamond_nr,
#     threads: 24
#     shell:
#         '''
#         diamond blastx -q {params.fa} -p {threads} -d {params.ref} \
#         -e 0.001 --id 40 -f tab -b 50 \
#         -o {params.out} &> {log}
#         '''
#
# rule diamond_eggnog_cross_megahit:
#     input:
#         bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
#     output:
#         bmk=    dir_out +"/co_assembly/log/diamond_eggnog_cross_megahit.bmk"
#     log:        dir_out +"/co_assembly/diamond_eggnog_cross_megahit/diamond_eggnog_cross_megahit.log",
#     benchmark:  dir_out +"/co_assembly/log/diamond_eggnog_cross_megahit.bmk"
#     params:
#         out=    dir_out +"/co_assembly/diamond_eggnog_cross_megahit/diamond_eggnog_cross_megahit",
#         fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
#         ref=    diamond_eggnog,
#     threads: 24
#     shell:
#         '''
#         diamond blastx -q {params.fa} -p {threads} -d {params.ref} \
#         -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
#         -o {params.out} &> {log}
#         '''
#
# rule diamond_kegg_cross_megahit:
#     input:
#         bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
#     output:
#         bmk=    dir_out +"/co_assembly/log/diamond_kegg_cross_megahit.bmk"
#     log:        dir_out +"/co_assembly/diamond_kegg_cross_megahit/diamond_kegg_cross_megahit.log",
#     benchmark:  dir_out +"/co_assembly/log/diamond_kegg_cross_megahit.bmk"
#     params:
#         out=    dir_out +"/co_assembly/diamond_kegg_cross_megahit/diamond_kegg_cross_megahit",
#         fa=     dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
#         ref=    diamond_kegg,
#     threads: 24
#     shell:
#         '''
#         diamond blastx -q {params.fa} -p {threads} -d {params.ref} \
#         -e 0.001 --id 30 -b 50 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send \
#         -o {params.out} &> {log}
#         '''

rule hmmsearch_cross:
    input:
        bmk=        dir_out +"/co_assembly/log/prokka_cross.bmk"
    output:
        out_hmm=    dir_out +"/co_assembly/log/hmmsearch_cross.{hmmId}.bmk"
    log:            dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.{hmmId}.log"
    benchmark:      dir_out +"/co_assembly/log/hmmsearch_cross.{hmmId}.bmk"
    params:
        hmm=        "gene_hmm/{hmmId}.hmm",
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
        --pfamtblout {params.pfamtblout} -E 10 -o {params.out_hmm} --cpu {threads} {params.hmm} {params.fa}        
        '''

rule rgi_prokka_cross:
    input:
        bmk=        dir_out +"/co_assembly/log/prokka_cross.bmk"
    output:         dir_out +"/co_assembly/log/rgi_prokka_cross.bmk"
    log:            dir_out +"/co_assembly/rgi_prokka_cross/rgi_prokka_cross.log"
    benchmark:      dir_out +"/co_assembly/log/rgi_prokka_cross.bmk"
    params:
        fa=         dir_out +"/co_assembly/prokka_cross/co_assembly.faa",
        out=        dir_out +"/co_assembly/rgi_prokka_cross/rgi_prokka_cross.out"
    threads: 32
    shell:
        '''
        # rgi load --card_json /scratch/zuhu/tmp/ref/metagenomics/CARD/Data_3.2.4/card.json --local
        # 
        # rgi card_annotation -i /scratch/zuhu/tmp/ref/metagenomics/CARD/Data_3.2.4/card.json 
        # 
        # rgi load \
        # -i /scratch/zuhu/tmp/ref/metagenomics/CARD/Data_3.2.4/card.json  \
        # --card_annotation /scratch/zuhu/tmp/ref/metagenomics/CARD/Data_3.2.4/card_database_v3.2.4.fasta --local
        # 
        # rgi wildcard_annotation -i ./ \
        # --card_json /scratch/zuhu/tmp/ref/metagenomics/CARD/Data_3.2.4/card.json -v 3.2.4
        
        rgi load \
        --wildcard_annotation /scratch/zuhu/tmp/ref/metagenomics/CARD/variants_4.0.0/wildcard_database_v3.2.4.fasta \
        --wildcard_index /scratch/zuhu/tmp/ref/metagenomics/CARD/variants_4.0.0/index-for-model-sequences.txt \
        --card_annotation /scratch/zuhu/tmp/ref/metagenomics/CARD/Data_3.2.4/card_database_v3.2.4.fasta --local
        
        rgi main --input_sequence {params.fa} \
          --output_file {params.out} --input_type protein --local \
          --alignment_tool DIAMOND --num_threads {threads} --clean --include_loose        
        '''

rule kraken2_prokka_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/prokka_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/kraken2_prokka_cross.bmk"
    log:        dir_out +"/co_assembly/kraken2_prokka_cross/kraken2_prokka_cross.log"
    benchmark:  dir_out +"/co_assembly/log/kraken2_prokka_cross.bmk"
    params:
        fa=         dir_out +"/co_assembly/prokka_cross/co_assembly.ffn",
        kraken2_db= kraken2_db,
        output_tsv= dir_out +"/co_assembly/kraken2_prokka_cross/kraken2_prokka_cross.output.tsv",
        report_tsv= dir_out +"/co_assembly/kraken2_prokka_cross/kraken2_prokka_cross.report.tsv",
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

rule prokka_cross:
    input:      dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    output:     dir_out +"/co_assembly/log/prokka_cross.bmk"
    log:        dir_out +"/co_assembly/log/prokka_cross.log"
    benchmark:  dir_out +"/co_assembly/log/prokka_cross.bmk"
    params:
        outprefix=  dir_out +"/co_assembly/prokka_cross",
        fa= dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
        id="co_assembly",
    threads: 24
    shell:
        '''
        prokka --force --cpus {threads} --outdir {params.outprefix} --prefix {params.id} --metagenome --kingdom Bacteria {params.fa}  &> {log}
        '''



rule kraken2_megahit_cross:
    input:
        bmk=    dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
    output:
        bmk=    dir_out +"/co_assembly/log/kraken2_megahit_cross.bmk"
    log:        dir_out +"/co_assembly/kraken2_megahit_cross/kraken2_megahit_cross.log"
    benchmark:  dir_out +"/co_assembly/log/kraken2_megahit_cross.bmk"
    params:
        fa=         dir_out +"/co_assembly/cd_hit_est_cross/cd_hit_est_cross.fa",
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

rule cd_hit_est_cross:
    input:      dir_out +"/co_assembly/log/megahit_cross.bmk"
    output:     dir_out +"/co_assembly/log/cd_hit_est_cross.bmk"
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
    input:      dir_out +"/co_assembly/log/megahit_cross.bmk"
    output:     dir_out +"/co_assembly/log/quast_cross.bmk"
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

# fq1_all=",".join(dir_fq_in+"/"+x+".1.fastq" for x in samplelist)
# fq2_all=",".join(dir_fq_in+"/"+x+".2.fastq" for x in samplelist)

rule megahit_cross:
    input:
        fq1=expand(dir_fq_in+"/{sample}.1.fastq",sample=samplelist),
        fq2=expand(dir_fq_in+"/{sample}.2.fastq",sample=samplelist)
    output:
        bmk=    dir_out +"/co_assembly/log/megahit_cross.bmk",
    log:        dir_out +"/co_assembly/log/megahit_cross.log"
    benchmark:  dir_out +"/co_assembly/log/megahit_cross.bmk"
    params:
        fq1=    ",".join(dir_fq_in+"/"+x+".1.fastq" for x in samplelist),
        fq2=    ",".join(dir_fq_in+"/"+x+".2.fastq" for x in samplelist),
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
