import pyranges as pr

# convert a DNA sequence into a peptide sequence
# Example usage
# dna_sequence = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAG"
# peptide_sequence = dna_to_peptide(dna_sequence)
def dna_to_peptide(dna_seq):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }    
    # Transcribing DNA to mRNA (replacing T with U)
    mrna_seq = dna_seq.replace('T', 'U')    
    # Splitting the mRNA sequence into codons
    codons = [mrna_seq[i:i+3] for i in range(0, len(mrna_seq), 3)]    
    # Translating the codons into a peptide sequence
    peptide_seq = ''.join([genetic_code.get(codon, 'X') for codon in codons if len(codon) == 3])    
    return peptide_seq

#get annotation
def get_anno(chr,start,end,genes_pr):
    positions_data = {'Chromosome': [chr],'Start': [int(start)],'End': [int(end)]}
    positions_pr = pr.PyRanges(pd.DataFrame(positions_data))
    pr_nearest= positions_pr.nearest(genes_pr, suffix='_gene')
    return pr_nearest

# get the reverse complement of a DNA sequence
# Example usage
# sequence = "ATCG"; print(reverse_complement(sequence))
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N': 'N'}
    reverse_comp = ''.join(complement[nucleotide] for nucleotide in reversed(dna_sequence))
    return reverse_comp


# function to detect mutation 
def detectMutation(dna, cdna): 
    count = 0      
    # x and y take each character in dna and cdna 
    # for character by character comparison 
    for x, y in zip(dna, cdna):          
        # if the character at the same index match 
        # then the count is increased 
        if x == y: 
            count = count + 1          
        # incase of mismatch the loop is broken 
        else: 
            break
          
    # the count value points to the index before the  
    # position of mutation 
    return count 

#search sequence in a seqeuence
def do_match(seq,seqs_compare_list):
    map_to_ref=[]
    start_of_refMap=[]
    for i in range(len(seqs_compare_list)):    
        seq_search=re.search(seqs_compare_list[i],seq)
        if seq_search is not None:
            # print(seqs_compare[0][i])
            map_to_ref.append(True)
            start_of_refMap.append(str(seq_search.start()))
    if any(map_to_ref): 
        map="Yes"
        start_="|".join(start_of_refMap)
    else: 
        map="No"
        start_=""
    return(map,start_)
    # return(map)


def get_flankingSeqs_position(chr,start0,fafile,flanking=10):
    ref_seq=fafile.fetch(chr,start0-flanking,start0+flanking+1)
    print(start0-flanking,start0+flanking+1)
    print(ref_seq,len(ref_seq))    
    ref_seq_rc=reverse_complement(ref_seq)
    return([ref_seq,ref_seq_rc])

def get_flankingSeqs_alt(chr,start0,altA,fafile,flanking=10):
    var_seq=fafile.fetch(chr,start0-flanking,start0)+ altA + fafile.fetch(chr,start0+1,start0+flanking+1)
    print(var_seq,len(var_seq))
    var_seq_rc=reverse_complement(var_seq)
    return([var_seq,var_seq_rc])

def get_flanking_refandVar_seq(chr,start0,altA,fafile,flanking=10):
    seq_ref=get_flankingSeqs_position(chr,start0,fafile,flanking=flanking)
    seq_var=get_flankingSeqs_alt(chr,start0,altA,fafile,flanking=flanking)   
    return(seq_ref,seq_var)
    
def parse_cigar(cigar_str):
    """Parse CIGAR string into a list of operations and their counts."""
    return re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)

def get_dfCigar(seq,cigar,align_start,length):
    cigar_seqs=[]
    i_seq=0
    i_ref=0
    end_ref=align_start
    for count, op in parse_cigar(cigar):
        slide_raw=int(count)
        
        slide_seq=np.where(op in "=MXISH",slide_raw,0)
        slide_ref=np.where(op in "=MXDN",slide_raw,0)
        
        start_ref=end_ref
        end_ref=end_ref+slide_ref
        start_seq=i_seq
        end_seq=i_seq+slide_seq
        seq_sub=seq[i_seq:end_seq]
        
        cigar_seqs.append([
            op,
            slide_raw,
            slide_seq,
            slide_ref,
            start_seq,
            end_seq,
            start_ref,
            end_ref,
            seq_sub
            ])
        i_seq+=slide_seq
        i_ref+=slide_ref
    
    df_cigar=pd.DataFrame(cigar_seqs,columns=["op","slide_raw","slide_seq","slide_ref","start_seq","end_seq","align_start","align_end","seq_sub"])
    return(df_cigar)

def get_varSeq(seq,df_cigar,start0,length):
    df_cigar1=df_cigar[(df_cigar['align_start'] <= start0) & (df_cigar['align_end'] > start0)]
    # print(df_cigar1)
        
    pos_match=np.where(df_cigar1.iloc[0,0] in "=MXP",start0-df_cigar1["align_start"]+df_cigar1['start_seq'],
                            np.where(df_cigar1.iloc[0,0] in "ISH",df_cigar1["start_seq"],
                                    np.where(df_cigar1.iloc[0,0] in "DN",df_cigar1["end_seq"],0))
            )
    op=df_cigar1.iloc[0,0]
    pos_match=int(pos_match[0])
    # list([seq[pos_match:pos_match+length],str(pos_match)])
    return(list([seq[pos_match:pos_match+length],str(pos_match),op]))