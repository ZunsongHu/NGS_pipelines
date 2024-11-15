

rule centrifuge_virus:
    input:
        # i1 = dir_fq_in+"/{sample}_R1.fastq.gz",
        # i2 = dir_fq_in+"/{sample}_R2.fastq.gz",
        i1= dir_fq_in + "/{sample}.R1.fq.gz",
        i2=dir_fq_in + "/{sample}.R2.fq.gz",
    output:
        classification = dir_out+"/{sample}/centrifuge_virus/{sample}.classification",
        report = dir_out+"/{sample}/centrifuge_virus/{sample}.report"
    log:        dir_out+"/{sample}/centrifuge_virus/{sample}.centrifuge_virus.log"
    benchmark:  dir_out+"/{sample}/log/{sample}.centrifuge_virus.bmk"
    params:
        centrifuge=centrifuge,
        ref_=dict_centrifuge_ref["virus"]
    threads: 8
    shell:
        '''
        {params.centrifuge} -x {params.ref_} -q -t -p {threads} -1 {input.i1} -2 {input.i2} \
        -S {output.classification} --report-file {output.report} &> {log}
        '''