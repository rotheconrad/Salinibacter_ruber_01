#!/usr/bin/env python

'''Calculates coverage per gene from tab blast and prodigal fasta file.

This tool takes the following input parameters:

    * Tabular Blast file containing results for 1 genome and 1 metagenome
    * Genome fasta file used as reference for blast search.
    * Prodigal fasta file of predicted gene nucleotide sequence:
      Gene name format is xxx_xxx_xxx_xxx_Contig#_Gene# # start # end # x
      This script uses the start and end positions provided by prodigal
      and the contig name/number to calculate the coverage.

This script returns the following files:

    * 2 col tsv output of gene_name \t coverage

This script requires the following packages:

    * argparse
    * collection.defaultdict
    * itertools

This file can also be imported as a module and contains the follwing 
functions:

    * read_fasta - parses fasta file format into name, seq
    * retrieve_gene_coverage - using gene position
    * calc_genome_coverage - using genome length + tabblast file
    * read_genome_lengths - reads genome length file returns dict
    * operator - coordinates the other three functions
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Wednesday, August 14th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import itertools
from collections import defaultdict


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def read_genome_lengths(rgf):
    """ Reads genome lengths file returns dict genome_name: length """

    rgf_tad = defaultdict(dict) # initialize dicts
    rgf_ani = defaultdict(dict)
    rgf_len = {}
    wg_len = 0

    # read through genome fasta and build dictionary of dictionary
    # Containing base pair position for length of each contig.
    with open(gf, 'r') as f:
        for name, seq in read_fasta(f):

            contig_name = name.split(' ')[0].split('_')[-1] 
            length = len(seq) # calculate length of contig

            rgf_len[contig_name] = length
            wg_len += length

            # This populates the dictionary with value of zero for each
            # base pair position for each contig in the genome fasta
            for i in range(1, length+1, 1):
                rgf_tad[contig_name][i] = 0
                rgf_ani[contig_name][i] = []

    return rgf_tad, rgf_ani, rgf_len, wg_len


def calc_genome_coverage(tbf, rgf_tad, rgf_ani, thrshld):
    """ Reads tabblast file and adds coverage by genome position """

    with open(tbf, 'r') as f:
        for l in f:
            # split each line and define variables of interest
            X = l.rstrip().split('\t')
            pident = float(X[2])
            contig_name = X[1].split('_')[-1]
            strt = min(int(X[8]), int(X[9]))
            stp = max(int(X[8]), int(X[9]))
            
            # for each read above the user specified threshold, add
            # coverage of +1 to each basepair position for length of
            # read alignment along the subject sequence.

            if pident >= thrshld:
                for i in range(strt, stp+1, 1):
                    rgf_tad[contig_name][i] += 1
                    rgf_ani[contig_name][i].append(pident)

    return rgf_tad, rgf_ani


def write_genome_cov_by_bp(rgf_tad, outpre):
    """ writes the tad coverage per bp as 2 col tsv file pos: value """

    counter = 1
    with open(f'{outpre}_genome_by_bp.tsv', 'w') as o:
        o.write('Position\tDepth\n')
        for k, v in rgf_tad.items():
            for i in v.values():
                o.write(f'{counter}\t{i}\n')
                counter += 1


def retrieve_gene_coverage(pgf, rgf_tad, rgf_ani):
    """ Retrieves list of depths for each bp position of gene length """

    gn_tad = defaultdict(list) # initialize dictionary
    gn_ani = defaultdict(list)
    gn_len = {}

    with open(pgf, 'r') as f:
        for name, seq in read_fasta(f):
            X = name.split(' # ')
            gene_name = X[0][1:]
            contig_name = gene_name.split('_')[-2]
            strt = min(int(X[1]), int(X[2]))
            stp = max(int(X[1]), int(X[2]))

            gn_len[gene_name] = len(seq)

            for i in range(strt, stp+1, 1):
                gn_tad[gene_name].append(rgf_tad[contig_name][i])
                gn_ani[gene_name].extend(rgf_ani[contig_name][i])

    return gn_tad, gn_ani, gn_len


def truncate(x, tad):
    """ returns tad range of a list/array """

    xsorted = sorted(x)
    xlen = len(x) # get length of the list
    inverse_tad = round((1.0 - tad) / 2.0, 2) # to get top and bottom 
    q = int(xlen * inverse_tad) # to get top and bottom
    bottom = q
    top = xlen - q
    tad_range = xsorted[bottom:top] # slice list

    #print(xlen, bottom, top, len(tad_range), sum(xsorted[:bottom+1]))

    return tad_range


def get_contig_tad(rgf_tad, tad):
    """ reads through rgf_tad and returns dict of tads by contig """

    contig_tad = {}
    wg_tad = []

    for k, v in rgf_tad.items():
        values = list(v.values())
        coverage = get_average(values, tad)
        contig_tad[k] = coverage
        wg_tad.extend(values)

    return contig_tad, wg_tad


def get_gene_tad(gn_tad, tad):
    """ reads through gn_tad and returns dict of tads by gene """

    gene_tad = {}

    for k, v in gn_tad.items():
        gene_tad[k] = get_average(v, tad)

    return gene_tad
    

def get_contig_anir(rgf_ani, tad):
    """ loops through in_d and calculates tad/ani for each key """

    contig_ani = {}
    wg_ani = []

    for k, v in rgf_ani.items():
        values = list(itertools.chain.from_iterable(list(v.values())))
        average = get_average(values, tad)
        if average > 0: contig_ani[k] = average
        wg_ani.extend(values)

    return contig_ani, wg_ani


def get_gene_anir(gn_ani, tad):
    """ loops through in_d and calculates tad/ani for each key """

    gene_ani = {}

    for k, v in gn_ani.items():
        average = get_average(v, tad)
        if average > 0: gene_ani[k] = average

    return gene_ani


def get_average(l, tad):
    """ returns anir from truncated list """

    trunc_val = truncate(l, tad)

    if sum(trunc_val) == 0:
        average = 0
    else:
        average = sum(trunc_val) / len(trunc_val)

    return average


def get_relative_abundance(wg_tad, mtg):
    """ calculates and returns relative abundance from wg_TAD """

    total_metagenome_bp = 0

    with open(mtg, 'r') as f:
        for name, seq in read_fasta(f):
            total_metagenome_bp += len(seq)

    relabndc = (sum(wg_tad) / total_metagenome_bp) * 100

    return relabndc, total_metagenome_bp


def write_file(in_d, len_d, outpre, outpost, precision):
    """ writes dictionary to file """
    
    outfile = outpre + outpost

    with open(outfile, 'w') as o:
        o.write('Name\tValue\tLength\n')
        for k, v in in_d.items():
            o.write(f'{k}\t{v:.{precision}f}\t{len_d[k]}\n')


def calc_tad_anir_relabndc(
                            mtg,
                            wglen,
                            rgf_tad,
                            rgf_ani,
                            rgf_len,
                            gn_tad,
                            gn_ani,
                            gn_len,
                            tad,
                            outpre,
                            precision
                            ):

    """ Calculate tad and anir for whole genome, contig, and gene """

    print('... Calculating TADs for Contigs')
    contig_tad, wg_tad = get_contig_tad(rgf_tad, tad)

    print('... Calculating TADs for Genes')
    gene_tad = get_gene_tad(gn_tad, tad)

    print('... Calculating TAD for Genome')
    wgtad = get_average(wg_tad, tad)

    print('... Calculating Total Metagenome Size & Relative Abundance')
    relabndc, total_metagenome_bp = get_relative_abundance(wg_tad, mtg)

    print('... Calculating ANI for Contigs')
    contig_ani, wg_ani = get_contig_anir(rgf_ani, tad)

    print('... Calculating ANI for Genes')
    gene_ani = get_gene_anir(gn_ani, tad)

    print('... Calculating ANI for Genome')
    wgani = get_average(wg_ani, tad)

    _ = write_file(contig_tad, rgf_len, outpre, '_contig_tad.tsv', precision)
    _ = write_file(contig_ani, rgf_len, outpre, '_contig_ani.tsv', precision)
    _ = write_file(gene_tad, gn_len, outpre, '_gene_tad.tsv', precision)
    _ = write_file(gene_ani, gn_len, outpre, '_gene_ani.tsv', precision)

    return wgtad, wgani, relabndc, total_metagenome_bp


def operator(mtg, rgf, tbf, pgf, thd, tad, outpre):
    """ Runs the different functions and writes out results """

    tadp = tad / 100
    precision = 2 # number of decimals places to keep.

    print(f'Using values {tad}% for TAD & {thd}% for ANIr')

    print('Preparing base pair array for each contig in genome.')
    rgf_tad, rgf_ani, rgf_len, wglen = read_genome_lengths(rgf)

    print('Calculating coverage for each base pair position in genome')
    rgf_tad, rgf_ani = calc_genome_coverage(tbf, rgf_tad, rgf_ani, thd)

    print('Writing Whole genome per base pair depth')
    _ = write_genome_cov_by_bp(rgf_tad, outpre)

    print('Retrieving coverage for each contig & gene')
    gn_tad, gn_ani, gn_len = retrieve_gene_coverage(pgf, rgf_tad, rgf_ani)

    print(f'Calculating {tad}% truncated average depth and {thd}% ANIr')
    (
        wgtad,
        wgani,
        relabndc,
        total_metagenome_bp
            ) = calc_tad_anir_relabndc(
                                        mtg,
                                        wglen,
                                        rgf_tad,
                                        rgf_ani,
                                        rgf_len,
                                        gn_tad,
                                        gn_ani,
                                        gn_len,
                                        tadp,
                                        outpre,
                                        precision
                                        )

    with open(f'{outpre}_genome.tsv', 'w') as o:
        o.write(
            f'Genome_Name\tTAD_{int(tad)}\tANIr_{int(thd)}\t'
            f'Relative_Abundance(%)\tGenome_Length(bp)\tMetagenome_Length(bp)\n'
            )
        o.write(
            f'{outpre}\t{wgtad:.{precision}f}\t{wgani:.{precision}f}%\t'
            f'{relabndc:.{precision}f}%\t{wglen}\t{total_metagenome_bp}\n'
            )

    print('\nScript seems to have finished successfully.\n')

    print('\nWhole Genome Values:\n')
    print(
        f'Genome_Name\tTAD_{int(tad)}\tANIr_{int(thd)}\t'
        f'Relative_Abundance(%)\tGenome_Length\tMetagenome_Length(bp)'
        )
    print(
        f"{outpre}\t{wgtad:.{precision}f}\t{wgani:.{precision}f}%\t"
        f"{relabndc:.{precision}f}%\t{wglen}\t{total_metagenome_bp}\n\n"
        )


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-mtg', '--metagenome_file',
        help='Please specify the query metagenome fasta file!',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-rgf', '--ref_genome_file',
        help='Please specify the genome fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-tbf', '--tabular_blast_file',
        help='Please specify the tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pgf', '--prodigal_gene_fasta_file',
        help='Please specify the prodigal gene fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-thd', '--pIdent_threshold_cutoff',
        help='Please specify pIdent threshold to use! (ie: 95)',
        metavar='',
        type=float,
        required=True
        )
    parser.add_argument(
        '-tad', '--truncated_avg_depth_value',
        help='Please specify TAD value! (ie: 80 or 90)',
        metavar='',
        type=float,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file_prefix',
        help='What do you like the output file prefix to be?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    operator(
            args['metagenome_file'],
            args['ref_genome_file'],
            args['tabular_blast_file'],
            args['prodigal_gene_fasta_file'],
            args['pIdent_threshold_cutoff'],
            args['truncated_avg_depth_value'],
            args['out_file_prefix']
            )


if __name__ == "__main__":
    main()
