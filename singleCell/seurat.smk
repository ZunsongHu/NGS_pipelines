rule seurat_integration_sct:
    input:
        expand(dir_analysis+"obj/individual/{sample}.SeuratObj.QC.rds",sample=samplelist)
    output:
        dir_analysis+"obj/obj_seurat.SCT.{label_integration}.rds"
    log:    dir_analysis+"obj/obj_seurat.SCT.{label_integration}.log"
    params:
        script=dir_R + "seurat/integration/seurat_integration.R",
        files_in=lambda wildcards: ",".join([dir_analysis+"obj/individual/"+i+".SeuratObj.QC.rds" for i in samplelist]),
        method="sct_cca"

    shell:
        """
        Rscript_seurat5 {params.script} \
        {params.files_in} \
        {output} \
        {params.method} &> {log}
        """

rule seurat_integration_harmony:
    input:
        expand(dir_analysis+"obj/individual/{sample}.SeuratObj.QC.rds",sample=samplelist)
    output:
        dir_analysis+"obj/obj_seurat.harmony.{label_integration}.rds"
    log:    dir_analysis+"obj/obj_seurat.harmony.{label_integration}.log"
    params:
        script=dir_R + "seurat/integration/seurat_integration.R",
        files_in=lambda wildcards: ",".join([dir_analysis+"obj/individual/"+i+".SeuratObj.QC.rds" for i in samplelist]),
        method="harmony"

    shell:
        """
        Rscript_seurat5 {params.script} \
        {params.files_in} \
        {output} \
        {params.method} &> {log}
        """

rule done_seurat_integration:
    input:
        expand(dir_analysis+"obj/obj_seurat.harmony.{label_integration}.rds",label_integration=label_integration),
        expand(dir_analysis+"obj/obj_seurat.SCT.{label_integration}.rds",label_integration=label_integration)
    output:
        dir_analysis+"obj/done_seurat_integration"
    shell:
        'touch {output}'
