rule wget:
    input:      file_url
    output:     dir_wget+"{sample}/{sample}.wget.bmk"
    benchmark:  dir_wget+"{sample}/{sample}.wget.bmk"
    log:        dir_wget+"{sample}/{sample}.wget.log"
    params:
        url=lambda wildcards: dict_url[wildcards.sample],
        filename=lambda wildcards: os.path.basename(dict_url[wildcards.sample]),
        dir_out=dir_wget+"{sample}/"
    threads: 1
    shell:
        '''
        if [ -f {params.filename} ]; then rm {params.filename}; fi
        wget {params.url}
        touch {params.filename} 
        mv {params.filename} {params.dir_out}
        '''

# rule tar:
#     input: dir_fqRaw+"/{sample}.tar.bz2"
#     output: dir_out+"/{sample}/log/{sample}.tar.bmk"
#     benchmark: dir_out+"/{sample}/log/{sample}.tar.bmk"
#     shell:
#         '''
#         tar -xf {input}
#         '''

rule get_fqLink:
    input: file_intake
    output: dir_fq_in+"/get_fqLink.bmk"
    benchmark: dir_fq_in+"/get_fqLink.bmk"
    params:
        dir_fqRaw=dir_fqRaw,
        dir_fqLink=dir_fq_in,
        rscript="/home/zgu_labs/bin/R/GeneralUsage/get_fq_softlink.R"
    shell:
        '''
        Rscript {params.rscript} {input} {params.dir_fqRaw} {params.dir_fqLink}
        cat ln.code|bash
        rm ln.code
        '''
    
rule bowtie_host:
    input:
        fq1 = dir_fq_in + "/{sample}.R1.fq.gz", fq2 = dir_fq_in + "/{sample}.R2.fq.gz",
    output:     dir_out + "/{sample}/log/{sample}.bowtie_host.bmk"
    benchmark:  dir_out + "/{sample}/log/{sample}.bowtie_host.bmk"
    log:        dir_out + "/{sample}/bowtie_host/{sample}.bowtie_host.log"
    params:
        ref= bowtie_ref_host,
        sam= dir_out +"/{sample}/bowtie_host/{sample}.sam"
    threads: 4
    shell:
        '''
        bowtie2 -x {params.ref} -1 {input.fq1} -2 {input.fq2} --very-sensitive -S {params.sam} &> {log}
        rm {params.sam}
        '''

rule kneaddata:
    input:
        fq1 = dir_fq_in + "/{sample}.R1.fq.gz", fq2 = dir_fq_in + "/{sample}.R2.fq.gz",
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.kneaddata.bmk",
    log:        dir_out +"/{sample}/kneaddata/{sample}.kneaddata.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.kneaddata.bmk",
    params:
        ref=kneaddata_ref,
        dir_out=dir_out +"/{sample}/kneaddata",
        output1=dir_out +"/{sample}/kneaddata/{sample}.R1.fq",
        output2=dir_out +"/{sample}/kneaddata/{sample}.R2.fq",
        dir_trimmomatic=kneaddata_trimmomatic,
        dir_fastqc=kneaddata_fastqc,
        output_prefix="{sample}"
    threads: 4
    shell:
        '''
        kneaddata --threads {threads} --max-memory 10000m \
        --input1 {input.fq1} --input2 {input.fq2} --output {params.dir_out} --reference-db {params.ref} \
        --output-prefix {params.output_prefix} \
        --bowtie2-options "--very-sensitive --dovetail" \
        --trimmomatic {params.dir_trimmomatic} --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
        --fastqc {params.dir_fastqc} --run-trim-repetitive --run-fastqc-start --run-fastqc-end \
        --remove-intermediate-output &> {log}
        # gzip {params.output1};gzip {params.output2}
        '''

rule fastqp:
    input:
        bmk=    dir_out +"/{sample}/log/{sample}.kneaddata.bmk"
    output:
        bmk=    dir_out +"/{sample}/log/{sample}.fastqp.bmk"
    log:        dir_out +"/{sample}/fastqp/{sample}.fastqp.log"
    benchmark:  dir_out +"/{sample}/log/{sample}.fastqp.bmk"
    params:
        fq_in1= dir_out +"/{sample}/kneaddata/{sample}.1_kneaddata_paired_1.fastq",
        fq_in2= dir_out +"/{sample}/kneaddata/{sample}.1_kneaddata_paired_2.fastq",
        fq_out1=dir_out +"/{sample}/fastqp/{sample}.1.fastq",
        fq_out2=dir_out +"/{sample}/fastqp/{sample}.2.fastq",
        html=dir_out +"/{sample}/fastqp/{sample}.fastqp.html",
        json=dir_out +"/{sample}/fastqp/{sample}.fastqp.json"
    threads: 4
    shell:
        '''
        fastp --thread {threads} -i {params.fq_in1} -I {params.fq_in2} -o {params.fq_out1} -O {params.fq_out2} \
        --html {params.html} --json {params.json} &> {log}
        '''



