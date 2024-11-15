dir_ref="/scratch/zuhu/tmp/ref"

#software --------------------------------------------------------------------------------------------------------------
centrifuge="/home/zuhu/bin/centrifuge/bin/centrifuge"


#reference database ----------------------------------------------------------------------------------------------------
dict_diamond_ref_db={"nr":      dir_ref+"diamond/nr.dmnd",
                     "eggnog":  dir_ref+"diamond/nr.dmnd",
                     "kegg":    dir_ref+"diamond/keggdb.dmnd",
                     "VFDB":    dir_ref+"diamond/VFDB.dmnd",
}

dict_diamond_ref_db_para={"nt":     "--id 40 --block-size 50",
                          "eggnog": "--id 30 --block-size 50",
                          "kegg":   "--id 30 --block-size 50",
                          "VFDB":   "--more-sensitive --max-target-seqs 1"
}

dict_centrifuge_ref={
    "virus":"/scratch/zuhu/tmp/xinrui/ref/virus/abv"
}

bowtie_ref_host="/ref_genomes/genome/human/GRCh38/bowtie2_index/index/genome"

kneaddata_ref="/scratch/zuhu/tmp/ref/metagenomics/kneaddata/human_genome/"
kneaddata_trimmomatic=" /home/zgu_labs/anaconda3/envs/ITD/share/trimmomatic_tmp/"
kneaddata_fastqc="/home/zgu_labs/anaconda3/opt/fastqc-0.11.9/"

#others ----------------------------------------------------------------------------------------------------------------

fa_for_diamond_build=dir_ref+"metagenomics/VFDB/VFDB_setB_pro.fas.gz"

kraken2_db=dir_ref+"kraken2/plusPF"
diamond_nr=dir_ref+"diamond/nr.dmnd"
diamond_eggnog=dir_ref+"diamond/eggnog.dmnd"
diamond_kegg=dir_ref+"diamond/keggdb.dmnd"


# key_file=dir_out +"/co_assembly/log/bowtieBuild_prokka_cross.bmk"
# key_file=dir_out +"/co_assembly/log/diamond_nr_cross_prokka.bmk"
# key_file=dir_out +"/co_assembly/log/diamond_eggnog_cross_prokka.bmk"
# key_file=dir_out +"/co_assembly/log/diamond_kegg_cross_prokka.bmk"
# key_file=dir_out +"/co_assembly/log/diamond_nr_cross_prokka.bmk"
# key_file=dir_out +"/{sample}/log/{sample}.bowtieMapping_prokka_cross.bmk"

# samplelist=list(pd.read_table("43.list",header=None)[0])

# hmmid_list=['dmd_tmd_align3','dmm_align3','mauA','mauB','tdm_align3','tmm_align3']

# dir_fq_in="fq_srs"
# dir_fq="fq_srs"

# fq1_all=",".join(dir_fq_in+"/"+x+".1.fastq" for x in samplelist)
# fq2_all=",".join(dir_fq_in+"/"+x+".2.fastq" for x in samplelist)
