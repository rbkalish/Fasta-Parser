genome = open("genome.fasta")
transcript = open("transcripts.fasta", "w")
tdict = {}
for line in genome:  #for loop parses genome.fasta
	if line.startswith(">"):
		id = line[1:9].strip()
	else: #puts info into tdict, sequence id = key, sequence = value
		try:
			tdict[id] += line.strip()
		except:
			tdict[id] = line.strip()
genome.close()
gtf = open("annotations.gtf")
gtf_dict = {}
for line in gtf:
	if line.startswith("#!"):
		pass
	else:
		tabs = line.split("\t")
		if tabs[2] == "CDS":
			header = tabs[0]
			start = tabs[3]
			stop = tabs[4]
			strand = tabs[6]
			attribute = tabs[8]
			#protein_id = attribute.split(";")[9].strip(" protein_id ").strip("\"")
			try:
				gtf_dict[header] += [[int(start), int(stop), strand]]
			except:
				gtf_dict[header] = [[int(start), int(stop), strand]]
gtf.close()
seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
sub_dict = {}
for key, annotations in gtf_dict.iteritems():
	sequence = tdict[key]
	for list in annotations:
		sub_seq = sequence[list[0]:list[1]]
		if "-" in list:  #reverse complements strands
			new_str = ""
			for nucleotide in sub_seq[::-1]:
				new_str += seq_dict.get(nucleotide)
			try:
				sub_dict[key] += [new_str]
			except:
				sub_dict[key] = [new_str]
		else:
			try:
				sub_dict[key] += [sub_seq]
			except:
				sub_dict[key] = [sub_seq]
for key, value in sub_dict.iteritems():
	transcript.write(">"), transcript.write(key), transcript.write("\n")
	for x in value:
		transcript.write(x), transcript.write("\n")
transcript.close()
subseq = open("transcripts.fasta")
codons = {}
for line in subseq:
	if line.startswith(">"):
		id = line.strip(">")
	else:
		sub = []
		for position in range(0,len(line)-2,3):
			sub.append(line[position:position+3].strip())
		try:
			codons[id.strip()] += [sub]
		except:
			codons[id.strip()] = [sub]		#random dinucleotides forming
aa_map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
proteins = {}
for key, sequences in codons.iteritems():
	for sequence in sequences:
		codonlist = []
		for codon in sequence:
			codonlist.append([aa_map.get(codon)])
		try:
			proteins[key] += [codonlist]
		except:
			proteins[key] = [codonlist]
test = open("proteins.fasta", "w")
for key, pros in proteins.iteritems():
	test.write(">"), test.write(key), test.write("\n")
	for protein in pros:
		for aa in protein:
			test.write(str(aa))
		test.write("\n")
test.close()
matrix = []
hydro_index = {'A' : 1.800,
	'R' : -4.500,
	'N' : -3.500,
	'D' : -3.500,
	'C' : 2.500,
	'E' : -3.500,
	'Q' : -3.500,
	'G' : -0.400,
	'H' : -3.200,
	'I' : 4.500,
	'L' : 3.800,
	'K' : -3.900,
	'M' : 1.900,
	'F' : 2.800,
	'P' : -1.600,
	'S' : -0.800,
	'T' : -0.700,
	'W' : -0.900,
	'Y' : -1.300,
	'V' : 4.200
}
for key, pros in proteins.iteritems():
	for protein in pros:
		length = len(protein)
		index = 0
		for aa in protein:
			if "STOP" in aa or None in aa:
				continue
			try:
				index += hydro_index[aa[0]] #key error when reaching stop codon
			except:
				index = hydro_index[aa[0]]
		matrix.append([key, length, (float(index)/length)])
gravy = open("gravy_score.tsv", "w")
gravy.write("".join(["sequence_id","	","length","	","gravy_score"])), gravy.write("\n")
for x in matrix:
	gravy.write("".join([x[0],"	",str(x[1]),"	",str(x[2])])), gravy.write("\n")
gravy.close()
