#this function will convert a strand of RNA into its DNA counterpart
def RNAtoDNA(RNA_strand):
    dna_list = []
    for i in range(len(RNA_strand)):
        if RNA_strand[i] == 'U':
            dna_list.append('T')
        else:
            dna_list.append(RNA_strand[i])
    dna = ''.join(dna_list)
    print(dna)


#this function converts a strand of DNA into its RNA counterpart
def DNAtoRNA(DNA_strand):
    rna_list = []
    for j in range(len(DNA_strand)):
        if DNA_strand[j] == 'T':
            rna_list.append("U")
        else:
            rna_list.append(DNA_strand[i])
        rna = ''.join(rna_list)
        print(rna)  
