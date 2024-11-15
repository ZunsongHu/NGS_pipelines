rule ContigMapping_kraken2_Prokka_cross:
    input:      dir_out +"/co_assembly/log/kraken2_prokka_cross.bmk"
    output:     dir_out +"/co_assembly/log/ContigMapping_kraken2_Prokka_cross.bmk"
    benchmark:  dir_out +"/co_assembly/log/ContigMapping_kraken2_Prokka_cross.bmk"
    log:        dir_out +"/co_assembly/kraken2_prokka_cross/anno/ContigMapping_kraken2_Prokka_cross.log"
    params:
        output_tsv= dir_out +"/co_assembly/kraken2_prokka_cross/kraken2_prokka_cross.output.tsv",
        report_tsv= dir_out +"/co_assembly/kraken2_prokka_cross/kraken2_prokka_cross.report.tsv",
        dir_out=    dir_out +"/co_assembly/kraken2_prokka_cross/anno/",
        rscript=dir_bin+"R/metagenomics/get_ContigMapping_kraken2.R"
    shell:
        '''
        /home/zuhu/bin/krakenReport_extraction.sh {params.report_tsv} {params.output_tsv} {params.dir_out}
        Rscript {params.rscript} {params.dir_out}output.mini.txt {params.dir_out}report.mini.txt {params.dir_out} &> {log}
        '''

rule ContigMapping_kraken2_Prokka:
    input:      dir_out +"/{sample}/log/{sample}.kraken2_prokka.bmk"
    output:     dir_out +"/{sample}/log/{sample}.ContigMapping_kraken2_Prokka.bmk"
    benchmark:  dir_out +"/{sample}/log/{sample}.ContigMapping_kraken2_Prokka.bmk"
    log:        dir_out +"/{sample}/kraken2_prokka/anno/{sample}.ContigMapping_kraken2_Prokka.log"
    params:
        output_tsv= dir_out +"/{sample}/kraken2_prokka/{sample}.kraken2_prokka.output.tsv",
        report_tsv= dir_out +"/{sample}/kraken2_prokka/{sample}.kraken2_prokka.report.tsv",
        dir_out=    dir_out +"/{sample}/kraken2_prokka/anno/",
        rscript=dir_bin+"R/metagenomics/get_ContigMapping_kraken2.R"
    shell:
        '''
        /home/zuhu/bin/krakenReport_extraction.sh {params.report_tsv} {params.output_tsv} {params.dir_out}
        Rscript {params.rscript} {params.dir_out}output.mini.txt {params.dir_out}report.mini.txt {params.dir_out} &> {log}
        '''

rule count_kraken2_Prokka_cross:
    input:
        bowtieMapping=  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross_rmDup.bmk",
        ContigMapping=  dir_out + "/co_assembly/log/ContigMapping_kraken2_Prokka_cross.bmk"
    output:             dir_out +"/{sample}/log/{sample}.count_kraken2_Prokka_cross.bmk",
    benchmark:          dir_out +"/{sample}/log/{sample}.count_kraken2_Prokka_cross.bmk",
    log:                dir_out +"/{sample}/count/kraken2_Prokka_cross/{sample}.kraken2_Prokka_cross.log"
    params:
        idx=            dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam.idx",
        ContigMapping=  dir_out +"/co_assembly/kraken2_prokka_cross/anno/df_contig_taxa.xlsx",
        dir_out=        dir_out +"/{sample}/count/kraken2_Prokka_cross/",
        rscript=        dir_bin+"R/metagenomics/count_kraken2_idx.R",
        sampleid=       "{sample}"
    threads: 1
    shell:
        '''
        Rscript {params.rscript} {params.idx} {params.ContigMapping} {params.dir_out} {params.sampleid} &> {log}
        '''

rule CountSummarise_kraken2_Prokka_cross:
    input:      expand(dir_out +"/{sample}/log/{sample}.count_kraken2_Prokka_cross.bmk",sample=samplelist)
    output:     dir_out +"/Summarise/log/CountSummarise_kraken2_Prokka_cross_{taxa}.bmk"
    benchmark:  dir_out +"/Summarise/log/CountSummarise_kraken2_Prokka_cross_{taxa}.bmk"
    log:        dir_out +"/Summarise/matrix/kraken2_Prokka_cross/CountSummarise_kraken2_Prokka_cross_{taxa}.log"
    params:
        prefix_in=   dir_out +"/*/count/kraken2_Prokka_cross/df_count_{taxa}.xlsx",
        file_list=  dir_out +"/Summarise/matrix/kraken2_Prokka_cross/{taxa}.list",
        prefix_out=    dir_out +"/Summarise/matrix/kraken2_Prokka_cross/matrix_{taxa}",
        rscript=    dir_bin+"R/metagenomics/get_count_matrix.R",
    shell:
        '''
        ls {params.prefix_in} > {params.file_list}
        Rscript {params.rscript} {params.file_list} {params.prefix_out} &> {log}
        rm {params.file_list}
        '''

rule ContigMapping_hmmsearch_cross:
    input:
        hmmsearch=      expand(dir_out +"/co_assembly/log/hmmsearch_cross.{hmmId}.bmk",hmmId=hmmid_list),
        ContigMapping_kraken=  dir_out +"/co_assembly/log/ContigMapping_kraken2_Prokka_cross.bmk"
    output:     dir_out +"/co_assembly/log/ContigMapping_hmmsearch_cross.bmk"
    benchmark:  dir_out +"/co_assembly/log/ContigMapping_hmmsearch_cross.bmk"
    log:        dir_out +"/co_assembly/hmmsearch_cross/ContigMapping_hmmsearch_cross.log"
    params:
        prefix_in=dir_out +"/co_assembly/hmmsearch_cross/hmmsearch_cross.*.tblout",
        dir_out=dir_out +"/co_assembly/hmmsearch_cross/",
        ContigMapping_kraken=dir_out +"/co_assembly/kraken2_prokka_cross/anno/df_contig_taxa.xlsx",
        rscript = dir_bin + "R/metagenomics/get_ContigMapping_hmmsearch.R",
    shell:
        '''
        ls {params.prefix_in} > {params.dir_out}tblout.list
        Rscript {params.rscript} {params.dir_out}tblout.list {params.ContigMapping_kraken} {params.dir_out} &> {log}
        rm {params.dir_out}tblout.list
        '''

rule count_hmmsearch_cross:
    input:
        bowtieMapping=  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross_rmDup.bmk",
        ContigMapping=  dir_out +"/co_assembly/log/ContigMapping_hmmsearch_cross.bmk"
    output:             dir_out +"/{sample}/log/{sample}.count_hmmsearch_cross.bmk",
    benchmark:          dir_out +"/{sample}/log/{sample}.count_hmmsearch_cross.bmk",
    log:                dir_out +"/{sample}/count/hmmsearch_cross/{sample}.count_hmmsearch_cross.log"
    params:
        idx=            dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam.idx",
        ContigMapping=  dir_out +"/co_assembly/hmmsearch_cross/df_contig_hmmsearchGene_taxa.xlsx",
        dir_out=        dir_out +"/{sample}/count/hmmsearch_cross/",
        rscript=        dir_bin+"R/metagenomics/count_hmmsearch_idx.R",
        sampleid=       "{sample}"
    threads: 1
    shell:
        '''
        Rscript {params.rscript} {params.idx} {params.ContigMapping} {params.dir_out} {params.sampleid} &> {log}
        '''

rule CountSummarise_hmmsearch_cross:
    input:      expand(dir_out +"/{sample}/log/{sample}.count_hmmsearch_cross.bmk",sample=samplelist)
    output:     dir_out +"/Summarise/log/CountSummarise_hmmsearch_cross.bmk"
    benchmark:  dir_out +"/Summarise/log/CountSummarise_hmmsearch_cross.bmk"
    log:        dir_out +"/Summarise/matrix/hmmsearch_cross/CountSummarise_hmmsearch_cross.log"
    params:
        prefix_in=   dir_out +"/*/count/hmmsearch_cross/df_count_hmmsearchGene.xlsx",
        file_list=  dir_out +"/Summarise/matrix/hmmsearch_cross/hmmsearchGene.list",
        prefix_out=    dir_out +"/Summarise/matrix/hmmsearch_cross/matrix_hmmsearchGene",
        rscript=    dir_bin+"R/metagenomics/get_count_matrix.R",
    shell:
        '''
        ls {params.prefix_in} > {params.file_list}
        Rscript {params.rscript} {params.file_list} {params.prefix_out} &> {log}
        rm {params.file_list}
        '''

rule ContigMapping_rgi_cross:
    input:
        rgi=      dir_out +"/co_assembly/log/rgi_prokka_cross.bmk",
        ContigMapping_kraken=  dir_out +"/co_assembly/log/ContigMapping_kraken2_Prokka_cross.bmk"
    output:     dir_out +"/co_assembly/log/ContigMapping_rgi_cross.bmk"
    benchmark:  dir_out +"/co_assembly/log/ContigMapping_rgi_cross.bmk"
    log:        dir_out +"/co_assembly/rgi_prokka_cross/ContigMapping_rgi_cross.log"
    params:
        rgi=dir_out +"/co_assembly/rgi_prokka_cross/rgi_prokka_cross.out.txt",
        rgi_mini=dir_out +"/co_assembly/rgi_prokka_cross/rgi_prokka_cross.mini.out",
        dir_out=dir_out +"/co_assembly/rgi_prokka_cross/",
        ContigMapping_kraken=dir_out +"/co_assembly/kraken2_prokka_cross/anno/df_contig_taxa.xlsx",
        rscript = dir_bin + "R/metagenomics/get_ContigMapping_rgi.R",
    shell:
        '''
        cat {params.rgi} | cut -f1,11,12,15,16,17,22,23 > {params.rgi_mini}
        Rscript {params.rscript} {params.rgi_mini} {params.ContigMapping_kraken} {params.dir_out}  &> {log}
        '''

rule count_rgi_cross:
    input:
        bowtieMapping=  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross_rmDup.bmk",
        ContigMapping=  dir_out +"/co_assembly/log/ContigMapping_rgi_cross.bmk"
    output:             dir_out +"/{sample}/log/{sample}.count_rgi_cross.bmk"
    benchmark:          dir_out +"/{sample}/log/{sample}.count_rgi_cross.bmk"
    log:                dir_out +"/{sample}/count/rgi_prokka_cross/{sample}.count_rgi_cross.log"
    params:
        idx=            dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam.idx",
        ContigMapping=  dir_out +"/co_assembly/rgi_prokka_cross/df_contig_rgi_taxa.xlsx",
        dir_out=        dir_out +"/{sample}/count/rgi_prokka_cross/",
        rscript=        dir_bin+"R/metagenomics/count_rgi_idx.R",
        sampleid=       "{sample}"
    threads: 1
    shell:
        '''
        Rscript {params.rscript} {params.idx} {params.ContigMapping} {params.dir_out} {params.sampleid} &> {log}
        '''

rule CountSummarise_rgi_cross:
    input:      expand(dir_out +"/{sample}/log/{sample}.count_rgi_cross.bmk",sample=samplelist)
    output:     dir_out +"/Summarise/log/CountSummarise_rgi_cross.bmk"
    benchmark:  dir_out +"/Summarise/log/CountSummarise_rgi_cross.bmk"
    log:        dir_out +"/Summarise/matrix/rgi_cross/CountSummarise_rgi_cross.log"
    params:
        prefix_in=   dir_out +"/*/count/rgi_prokka_cross/df_count_AMR.xlsx",
        file_list=  dir_out +"/Summarise/matrix/rgi_cross/rgi_AMR.list",
        prefix_out=    dir_out +"/Summarise/matrix/rgi_cross/matrix_rgiAMR",
        rscript=    dir_bin+"R/metagenomics/get_count_matrix.R",
    shell:
        '''
        ls {params.prefix_in} > {params.file_list}
        Rscript {params.rscript} {params.file_list} {params.prefix_out} 
        rm {params.file_list}
        '''

rule ContigMapping_diamond_cross:
    input:
        diamond=      dir_out +"/co_assembly/log/diamond_prokka_cross.{diamond_db}.bmk",
        ContigMapping_kraken=  dir_out +"/co_assembly/log/ContigMapping_kraken2_Prokka_cross.bmk"
    output:     dir_out +"/co_assembly/log/ContigMapping_diamond.{diamond_db}_cross.bmk"
    benchmark:  dir_out +"/co_assembly/log/ContigMapping_diamond.{diamond_db}_cross.bmk"
    log:        dir_out +"/co_assembly/diamond_prokka_cross/{diamond_db}/ContigMapping_diamond_cross.log"
    params:
        diamond=dir_out +"/co_assembly/diamond_prokka_cross/{diamond_db}/diamond_prokka_cross.{diamond_db}.txt",
        dir_out=dir_out +"/co_assembly/diamond_prokka_cross/{diamond_db}/",
        ContigMapping_kraken=dir_out +"/co_assembly/kraken2_prokka_cross/anno/df_contig_taxa.xlsx",
        rscript = dir_bin + "R/metagenomics/get_ContigMapping_diamond.R",
    shell:
        '''
        Rscript {params.rscript} {params.diamond} {params.ContigMapping_kraken} {params.dir_out}  &> {log}
        '''

rule count_diamond_cross:
    input:
        bowtieMapping=  dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross_rmDup.bmk",
        ContigMapping=  dir_out +"/co_assembly/log/ContigMapping_diamond.{diamond_db}_cross.bmk"
    output:             dir_out +"/{sample}/log/{sample}.count_diamond_cross.{diamond_db}.bmk"
    benchmark:          dir_out +"/{sample}/log/{sample}.count_diamond_cross.{diamond_db}.bmk"
    log:                dir_out +"/{sample}/count/diamond_prokka_cross.{diamond_db}/{sample}.count_diamond_cross.{diamond_db}.log"
    params:
        idx=            dir_out + "/{sample}/get_abundance/bowtieMapping_prokka_cross/{sample}.rmDup.bam.idx",
        ContigMapping=  dir_out +"/co_assembly/diamond_prokka_cross/{diamond_db}/df_contig_diamond.{diamond_db}_taxa.xlsx",
        dir_out=        dir_out +"/{sample}/count/diamond_prokka_cross.{diamond_db}/",
        rscript=        dir_bin+"R/metagenomics/count_diamond_idx.R",
        sampleid=       "{sample}"
    threads: 1
    shell:
        '''
        Rscript {params.rscript} {params.idx} {params.ContigMapping} {params.dir_out} {params.sampleid} &> {log}
        '''

rule CountSummarise_diamond_cross:
    input:      expand(dir_out +"/{sample}/log/{sample}.count_diamond_cross.{diamond_db}.bmk",sample=samplelist,diamond_db=diamond_db_list)
    output:     dir_out +"/Summarise/log/CountSummarise_diamond_cross.{diamond_db}.bmk"
    benchmark:  dir_out +"/Summarise/log/CountSummarise_diamond_cross.{diamond_db}.bmk"
    log:        dir_out +"/Summarise/matrix/diamond_cross/CountSummarise_diamond_cross.{diamond_db}.log"
    params:
        prefix_in=   dir_out +"/*/count/diamond_prokka_cross.{diamond_db}/df_count_diamond.{diamond_db}.xlsx",
        file_list=  dir_out +"/Summarise/matrix/diamond_cross/diamond.{diamond_db}.list",
        prefix_out=    dir_out +"/Summarise/matrix/diamond_cross/matrix_diamond.{diamond_db}",
        rscript=    dir_bin+"R/metagenomics/get_count_matrix.R",
    shell:
        '''
        ls {params.prefix_in} > {params.file_list}
        Rscript {params.rscript} {params.file_list} {params.prefix_out} 
        rm {params.file_list}
        '''


























