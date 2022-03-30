#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 15:33:25 2020
@author: christospapadopoulos

This tool extracts the transcriptome of a given genome and ellongates towards the
5`and 3`UTRs of all the genes for 50 nucletides. This tool is to be used when no
5`and 3`UTRs have been annotated on the GFF file of the organism.

This tool is used in ORFribo in order to map correctly all the reads on the transcriptome
and do not loose reads at the begining and the end of the genes.

Ex. of use:
ORFelongate -fna ${genome}.fasta -gff ${annotation}.gff -features_include CDS -type nucl -o transcriptome_ellongated
"""

import argparse,os,random
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio import SeqIO
import re
import linecache


def get_args():
    """
    Returns:
        Parameters
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Calculation')
    parser.add_argument("-fna",
                        type=str,
                        #action='store',
                        required=True,
                        nargs="?",
                        help="Genomic fasta file (.fna) ")
    parser.add_argument("-gff",
                        type=str,
                        #action='store',
                        required=True,
                        nargs="?",
                        help="Annotation file (.gff) ")
    parser.add_argument("-o",
                        type=str,
                        #action='store',
                        required=True,
                        nargs="?",
                        help="Output name of fasta file(s)")
    parser.add_argument("-features_include",
                        type=str,
                        action='store',
                        required=False,
                        nargs="*",
                        default = ['all'],
                        help="Annotation features to be considered (By definition is all)")
    parser.add_argument("-features_exclude",
                        type=str,
                        action='store',
                        required=False,
                        nargs="*",
                        default = ["None"],
                        help="Annotation features not to be considered (By definition is None)")
    parser.add_argument("-chr_exclude",
                        type=str,
                        action='store',
                        required=False,
                        nargs="?",
                        default = [],
                        help="Chromosomes to be excluded")
    parser.add_argument("-N",
                        type=int,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = False,
                        help="Size of the sample to be generated")
    parser.add_argument("-type",
                        type=str,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = "prot",
                        help="Sequences type of output [prot (default) / nucl / both];")

    args = parser.parse_args()
    return args


def read_genome(genome_fasta):
    with open(genome_fasta, 'r') as fasta_file:
        temp = ''.join([line.rstrip() if not line.startswith('>') else '\n' + line for line in fasta_file.readlines()]).split('\n')
        genome = {temp[i].replace('>', '').split()[0] : temp[i+1] for i in range(1, len(temp), 2)}
    return(genome)



def read_gff_info(gff,elements_in,elements_out):
    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    dico_info = {}
    with open(gff,'r') as f:
        print('\n')
        count = 0
        for x,line in enumerate(f):
            if line.startswith('#') == False:
                # And now read the line and extract the sequence
                element= line.split()[2].rstrip()

                if re.match(elements_out,element) and not re.match(elements_in,element):
                    continue

                if re.match(elements_out,element) and re.match(elements_in,element):
                    print("{} include and exclude at the same time!".format(element))
                    exit()

                if not re.match(elements_in,element) and elements_in != "(all)":
                    continue

                chrom  = line.split()[0].rstrip()
                if chrom in parameters.chr_exclude:
                    continue
                strand = line.split()[6].rstrip()
                start  = int(line.split()[3])
                stop   = int(line.split()[4])
                info   = line.split()[8]
                frame  = line.split()[7].rstrip()
                # In case CDS has no "ID" or "Name" attribute
                if info.split(';')[0].split('=')[0].rstrip() == "ID" or info.split(';')[0].split('=')[0].rstrip() == "Name":
                    gene   = info.split(';')[0].split('=')[1].rstrip()
                else:
                    continue

                if gene not in dico_info.keys():
                    count += 1
                    print('\r\t' + 'GFF' + '\t:\t' + str(count)+'\tsequences read', end = '')
                    dico_info[gene] = {}
                    dico_info[gene]['positions']=[]
                    dico_info[gene]['positions'].append(start)
                    dico_info[gene]['positions'].append(stop)
                    dico_info[gene]['DNA_seq'] = ''
                    dico_info[gene]["DNA_parts"] = []
                    gene_seen = False
                else:
                    gene_seen = True

                dico_info[gene]['strand'] = strand
                dico_info[gene]['chrom']  = chrom
                dico_info[gene]['info']   = info
                dico_info[gene]['frame']  = frame


                if strand == '+':

                    seq_dna = genome[chrom][start-1:stop]
                    my_seq_dna = Seq(seq_dna.upper())


                    if gene_seen == False:
                        dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                    elif gene_seen == True:

                        if stop < min(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq']

                        elif start > max(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                        else:
                            print(gene)
                            continue


                    dico_info[gene]['positions'].append(start)
                    dico_info[gene]['positions'].append(stop)

                    try:
                        dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')

                    except:
                        dico_info[gene]['AA_seq']  = ''

                elif strand == '-':
                    seq_dna = genome[chrom][start-1:stop]
                    my_seq_dna = Seq(seq_dna.upper())
                    my_seq_dna = my_seq_dna.reverse_complement()
                    dico_info[gene]["DNA_parts"].append(my_seq_dna)

                    if gene_seen == False:
                        dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                    elif gene_seen == True:
                        if stop < min(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)

                        elif start > max(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq']

                        else:
                            print(gene)
                            continue


                        dico_info[gene]['positions'].append(start)
                        dico_info[gene]['positions'].append(stop)


                    try:
                        dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')

                    except:
                        dico_info[gene]['AA_seq']  = ''



    return(dico_info)


def write_multifastas(dico_info,outname):
    names = []
    if parameters.type == "nucl" or parameters.type == "both":
        fwn = open(outname + '.nfasta', 'w')
        fwg = open(outname + '.gff','w')

    for gene in dico_info:
        try:
            if dico_info[gene]['trans_dna_seq'] != '':
                fwn.write('>'+gene.replace("CDS","mRNA")+'\n')
                fwn.write(dico_info[gene]['trans_dna_seq']+'\n')

                fwg.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene.replace("CDS","mRNA"),"ELL","gene","1", str((len(dico_info[gene]['trans_dna_seq']))),".","+",".","ID=" + gene.replace("CDS","") + ";"))
                fwg.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene.replace("CDS","mRNA"),"ELL","mRNA","1", str((len(dico_info[gene]['trans_dna_seq']))),".","+",".","ID=" + gene.replace("CDS","mRNA") + ";Parent=" + gene.replace("CDS","")))
                fwg.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene.replace("CDS","mRNA"),"ELL","CDS","51", str((len(dico_info[gene]['trans_dna_seq']))-50),".","+","0","ID=" + gene + ";Parent=" + gene.replace("CDS","mRNA")))
                fwg.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene.replace("CDS","mRNA"),"ELL","five_prime_UTR","1","50",".","+","0","ID=" + gene.replace("CDS","5UTR") + ";Parent=" + gene.replace("CDS","mRNA")))
                fwg.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene.replace("CDS","mRNA"),"ELL","three_prime_UTR",str((len(dico_info[gene]['trans_dna_seq']))-49),str((len(dico_info[gene]['trans_dna_seq']))),".","+","0","ID=" + gene.replace("CDS","3UTR") + ";Parent=" + gene.replace("CDS","mRNA")))

        except:
            pass

    try:
        fwn.close()
    except:
        pass
    try:
        fwg.close()
    except:
        pass

def find_matches(gff , elements_in , elements_out):
    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    matches = []
    with open(gff,'r') as f:
        for x,line in enumerate(f):
            if line.startswith('#') == False:
                element= line.split()[2].rstrip()

                if re.match(elements_out,element) and not re.match(elements_in,element):
                    continue

                if re.match(elements_out,element) and re.match(elements_in,element):
                    print("{} include and exclude at the same time!".format(element))
                    exit()

                if not re.match(elements_in,element) and elements_in != "(all)":
                    continue

                matches.append(x)

    return matches



def read_matched_gff_lines(matches , gff_file , genome):
    count = 0
    dico_info = {}
    for line in matches:
        my_line = linecache.getline(gff_file,line)
        my_line = my_line.split("\t")
        strand = my_line[6].rstrip()
        start  = int(my_line[3])
        stop   = int(my_line[4])
        chrom  = my_line[0].rstrip()
        info   = my_line[8]
        gene   = info.split(';')[0].split('=')[1].rstrip()

        if gene not in dico_info.keys():
            count += 1
            print('\r\t' + 'GFF' + '\t:\t' + str(count)+'\tsequences read', end = '')
            dico_info[gene] = {}
            dico_info[gene]['positions']=[]
            dico_info[gene]['positions'].append(start)
            dico_info[gene]['positions'].append(stop)
            dico_info[gene]['DNA_seq'] = ''
            dico_info[gene]["DNA_parts"] = []
            gene_seen = False
        else:
            gene_seen = True

        dico_info[gene]['strand'] = strand
        dico_info[gene]['chrom']  = chrom
        dico_info[gene]['info']   = info


        if strand == '+':

            seq_dna = genome[chrom][start-1:stop]
            my_seq_dna = Seq(seq_dna.upper())


            if gene_seen == False:
                dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
            elif gene_seen == True:

                if stop < min(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq']

                elif start > max(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                else:
                    print(gene)
                    continue


            dico_info[gene]['positions'].append(start)
            dico_info[gene]['positions'].append(stop)

            try:
                dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')

            except:
                dico_info[gene]['AA_seq']  = ''

        elif strand == '-':
            seq_dna = genome[chrom][start-1:stop]
            my_seq_dna = Seq(seq_dna.upper())
            my_seq_dna = my_seq_dna.reverse_complement()
            dico_info[gene]["DNA_parts"].append(my_seq_dna)

            if gene_seen == False:
                dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
            elif gene_seen == True:
                if stop < min(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)

                elif start > max(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq']

                else:
                    print(gene)
                    continue


                dico_info[gene]['positions'].append(start)
                dico_info[gene]['positions'].append(stop)


            try:
                dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')

            except:
                dico_info[gene]['AA_seq']  = ''

    return(dico_info)


def ellongation_3_5_UTR(dico):
    for i in dico:
        dico[i]['trans_positions'] = []
        dico[i]['trans_positions'].append(min(dico[i]['positions'])-50)
        dico[i]['trans_positions'].append(max(dico[i]['positions'])+50)

        if dico[i]['strand'] == "+":
            trans_seq = genome[dico[i]['chrom']][min(dico[i]['trans_positions']):min(dico[i]['positions'])].upper() + dico[i]['DNA_seq'] + genome[dico[i]['chrom']][max(dico[i]['positions']):max(dico[i]['trans_positions'])].upper()
            dico[i]["trans_dna_seq"] = trans_seq

        if dico[i]["strand"] == "-":
            trans_seq = str(Seq(genome[dico[i]['chrom']][max(dico[i]['positions']):max(dico[i]['trans_positions'])].upper()).reverse_complement()) + \
                 dico[i]['DNA_seq'] + \
                 str(Seq(genome[dico[i]['chrom']][min(dico[i]['trans_positions']):min(dico[i]['positions'])].upper()).reverse_complement())
            dico[i]["trans_dna_seq"] = trans_seq
    return dico





#def main():
global parameters
parameters = get_args()

print('''
      Genome file              : {}
      GFF file                 : {}
      Features to keep         : {}
      Features to exclude      : {}
      Chromosome to exclude    : {}
      Sample of population     : {}
      ''' . format(parameters.fna , parameters.gff , parameters.features_include , parameters.features_exclude , parameters.chr_exclude ,parameters.N))

genome_file = parameters.fna
global genome
genome  = read_genome(genome_fasta = genome_file)

elements_in = parameters.features_include

chomosomes_exclude = parameters.chr_exclude

gff_file = parameters.gff

if parameters.N != False:
    #print("Reading sequences and generating sample of {}".format(parameters.N))
    matches = find_matches(gff = gff_file,
                           elements_in =parameters.features_include,
                           elements_out=parameters.features_exclude )

    matches = sorted(random.sample(k=parameters.N , population= matches))

    seqs = read_matched_gff_lines(matches , gff_file = gff_file , genome = genome)

else:
    seqs = read_gff_info(gff=gff_file,
                         elements_in=parameters.features_include,
                         elements_out=parameters.features_exclude)



    seqs_ellongated = ellongation_3_5_UTR(dico = seqs)

#outname = os.path.basename(gff_file)
#outname = os.path.splitext(outname)[0]
#outname = outname + "_" + type_of_data
#outname  = parameters.o


write_multifastas(dico_info = seqs_ellongated , outname=parameters.o)

print()
