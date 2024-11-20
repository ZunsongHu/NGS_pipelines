singleRref = pd.read_table(config["file_singleRref"])
dict_singleRfile=dict(zip(singleRref.label,singleRref.ref_file))
dict_singleRcol=dict(zip(singleRref.label,singleRref.col_file))
dict_singleRannolable=dict(zip(singleRref.label,singleRref.anno_label))

rule SingleR_integration:
    input:      dir_analysis+"obj/obj_seurat.SCT.{label_integration}.rds"
    output:     dir_analysis+"singleR/sct/SingleR_integration.{label_integration}.{singleR_label}.bmk"
    benchmark:  dir_analysis+"singleR/sct/SingleR_integration.{label_integration}.{singleR_label}.bmk"
    log:        dir_analysis+"singleR/sct/SingleR_integration.{label_integration}.{singleR_label}.log"
    params:
        sample=     "{label_integration}",
        dir_out=    dir_analysis+"singleR/sct/",
        label=      "{singleR_label}",
        ref_rds=    lambda wildcards: dict_singleRfile[wildcards.singleR_label],
        anno_label= lambda wildcards: dict_singleRannolable[wildcards.singleR_label],
        col_file=   lambda wildcards: dict_singleRcol[wildcards.singleR_label],
        transform=  "SCT"
    threads: 1
    shell:
        '''
        Rscript_seurat5  /net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/singleCell/scripts/R/singleR/run_singleR.R \
        {input} {params.dir_out} {params.sample} \
        {params.label} {params.ref_rds} {params.anno_label} {params.col_file} {params.transform} 
        '''

rule SingleR_annotation_1Ref_done:
    input:      lambda wildcards: expand(dir_out+"{sample}/log/{sample}.{singleR_label}.SingleR_annotation_1Ref.bmk",sample=wildcards.sample,singleR_label=singleR_label_list)
    output:     dir_out+"{sample}/log/SingleR_annotation_1Ref.done"
    shell:
        '''
        touch {output}
        '''