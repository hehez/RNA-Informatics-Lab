#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 05/30/17 12:56AM
# @Author  : hongfei
# @File    : fr3d_search.py
# @Version : 0.3.0

"""Support following arguments.

None
    Output the whole sequence of all chains of each pdb file in PDBs folder.
pdbid
    print the whole sequence of all chains by a given certain pdbid.
pdbid chainid
    print the whole sequence
pdbid chainid index_start index_end
    print the sliced sequence
"""


# Import the modules needed to run the script.
import sys
import os
import numpy as np
import itertools
import gzip

from prody import *

## mapping dictionary
AAMAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
    'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O', 'XLE': 'J',
    'PTR': 'Y', 'TPO': 'T', 'SEP': 'S', 'CSO': 'C', 'HSD': 'H', 'HSP': 'H', 'HSE': 'H'
}

## update mapping dictionary
AAMAP.update({
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'U': 'U',
    'T': 'T',
    'DA': 'A',
    'DC': 'C',
    'DG': 'G',
    'DU': 'U',
    'DT': 'T',
})
## update mapping dictionary
with open('fr3d_mapping', 'r') as rf:
    AAMAP.update([line.split() for line in rf])

## update mapping dictionary
AAMAP.update({
    'MSE': 'M', 'NH2': 'X', 'UNK': 'X', 'PYY': 'N', 'PED': 'N', 'PHA': 'F', 'QUA': 'X', 'DHA': 'S', 'BB9': 'C', 'DBU': 'T', 'DCY': 'C', 'TS9': 'I', 'MH6': 'S',
    'N': 'N', 'MHW': 'X', 'DBB': 'T', 'MHU': 'F', 'MHV': 'X', '004': 'X', 'MHT': 'X', 'P1P': 'N', 'NF2': 'N', 'DPR': 'P', 'DAR': 'R', 'ZUK': 'X', 'S9L': 'X',
    '3GL': 'E', 'MLZ': 'K', '1DP': 'X', 'ACE': 'X', 'ILX': 'N', 'TRX': 'P', 'CSX': 'I', 'HYP': 'C', 'ACA': 'X', 'HFA': 'X', 'BTN': 'X', 'MEA': 'F', 'TSE': 'X',
    'PAE': 'X', '2OP': 'X', 'PO2': 'X', 'CMT': 'C', '3DR': 'N',
    ## 3J5T
    'KBE': '',
    'DPP': '',
    'UAL': '',
    '5OH': '',
    'MYN': '',
    'PO4': 'X',
})

PDB_DB = '/home/common/PDBs/'
DEBUG = 0
## True if wants to parse PDB file manually, otherwise uses prody package
MANUAL_PARSE_ON = 1

def extract_nuclotide_seq_indices(pdbid, chainid, s, e):
    '''This method is used to get a sequence with pdb id, chain id and indices.

    Args:
        pdbid: Case insensitive
        chainid: Case sensitive
        s: Start position of sequence
        e: End position of sequence

    Returns:
        A substring sequence with given indices.

    Raises:
        No, only print error info
    '''
    try:
        if MANUAL_PARSE_ON:
            seq = extract_nuclotide_seq_manual_parse(pdbid, chainid, 'seq')
        else:
            seq = extract_nuclotide_seq(pdbid, chainid)

        if s == 0 and e == 0:
            ## full sequence
            return seq
        else:
            ## either s or e is 0
            if s == 0:
                ## start 0 to end index
                return seq[:e]
            if e == 0:
                ## start index to end
                return seq[s-1:]
        return seq[s-1:e]
    except Exception, e:
        print 'extract_nuclotide_seq_indices >>>', e


def extract_nuclotide_seq_manual_parse(pdbid, chain=None, retMode='chain'):
    '''Manually parse the SEQRES entries in a pdb file. If this fails, use the ATOM
    entries(unimplemented). Return either a list of sequences or single sequence.

    Args:
        pdbid: Case insensitive
        chain: Case sensitive, None return a list of sequences, otherwise a sequence with given chain id
        retMode:
                'chain' as default, return sequence along with its chain number, i.e A   DGDADGDADGDGDADADGDADGDADA
                'seq', return only sequence
    Returns:
        A list of sequences: if chain id is not given, the list contains sequences of all chain
        single sequence: if chain id is given

    Raises:
        No, only print error info
    '''

    try:
        seq_type = 'SEQRES'
        f_pdb = os.path.join(PDB_DB, pdbid + '.pdb.gz')
        pdbIsExist = os.path.isfile(f_pdb)
        if not pdbIsExist:
            fetchPDB(pdbid)

        f = gzip.open(f_pdb, 'r')
        pdb = f.readlines()
        f.close()
        seq = [l for l in pdb if l[0:6] == seq_type]
        chain_dict = dict([(l[11],[]) for l in seq])
        if chain:
            if chain_dict.get(chain) is None: return ''
            ## single chain
            chain_seq = [l[19:70].split() for l in seq if l[11] == chain]
            for fragment in chain_seq:
                chain_dict[chain].extend([AAMAP.get(e) for e in fragment])
            if 'chain' == retMode:
                return'%s\t%s' %(chain, ''.join(chain_dict.get(chain)))
            else:
                return ''.join(chain_dict.get(chain))
        else:
            ## all chain
            result = []
            for c in chain_dict.keys():
                chain_seq = [l[19:70].split() for l in seq if l[11] == c]
                for fragment in chain_seq:
                    chain_dict[c].extend([AAMAP.get(e) for e in fragment])

                if 'chain' == retMode:
                    result.append('%s\t%s' %(c, ''.join(chain_dict.get(c))))
                else:
                    result.append(''.join(chain_dict.get(c)))
            return result
    except Exception, e:
        print 'extract_nuclotide_seq_manual_parse >>>', e


def extract_nuclotide_seq(pdbid, chainid, atomgroup=None, multi_ret=False):
    '''Using prody package to parse a pdb file.

    PS: pdbter will find the final index, in most case it gets a full of sequence
    unless the indices numbers in reverse order, so see method

        extract_nuclotide_seq_manual_parse(pdb_file_path)

    Args:
        pdbid: Case insensitive
        chain: Case sensitive
        atomgroup: a class after parsing pdb file, None as default
        multi_ret: False as default, True is return both sequence and a list of
                    each element of this sequence. For debug usage

    Returns:
        If multi_ret is true, return a sequence along with a list of each element of
        this sequence. Otherwise, only return sequence.

    Raises:
        No, only print error info
    '''

    try:
        pdb = atomgroup[chainid] if atomgroup else parsePDB(pdbid, chain=chainid)
        result = ''
        tmp = []
        if pdb is not None:
            if pdb.select('pdbter'):
                ter = pdb.select('pdbter').getResnums()[0]
                # ter = pdb[chainid].select('pdbter').getResnums()[0]
                indices = np.unique(pdb.getResnums())
                nindices = indices[indices <= ter]
                params = ['resnum %s' % ('`%s`' % index if index < 0 else index) for index in nindices]
                # result = ''.join([np.unique(pdb.select(p).getResnames()[0][-1:]).astype('|S1').tostring() for p in params])
                # for c in pdb.iterChains():
                #     print '===\t\t',c.getSequence()

                for p in params:
                    syb = np.unique(pdb.select(p).getResnames())[0].tostring()
                    if syb in AAMAP.keys():
                        tmp.append(AAMAP.get(syb))
                    else:
                        tmp.append(syb)
                result = ''.join(tmp)
        else:
            print 'NO result for case: %s %s' % (pdbid, chainid)
        if multi_ret:
            return result, tmp
        return result
    except Exception, e:
        print 'extract_nuclotide_seq >>>', e


def isExist(pdbid):
    '''Check if a pdb file exists in the folder.

    Args:
        pdbid: Case insensitive

    Returns:
        True if a pdb file exist, otherwise false.
    '''

    return os.path.isfile('%s%s.pdb.gz' % (PDB_DB, pdbid)) or os.path.isfile('%s%s.pdb' % (PDB_DB, pdbid))


def extract_nuclotide_all_seqs_frm_one(pdbid):
    '''This method is used to get all sequence from one pdb id.

    Args:
        pdbid: Case insensitive

    Returns:
        A list of sequences.
    '''

    result = []
    if MANUAL_PARSE_ON:
        seqs = extract_nuclotide_seq_manual_parse(pdbid, None, 'chain')
        if DEBUG:
            pass
        else:
            result = ['%s\t%s\n' % (pdbid, chain_plus_seq) for chain_plus_seq in seqs]
    else:
        pdb = parsePDB(pdbid)
        if DEBUG:
            for chain in pdb.iterChains():
                print '%s\t%s\t%s\n' % (pdbid, chain.getChid(), extract_nuclotide_seq(pdbid, chain.getChid(), pdb))
        else:
            result = ['%s\t%s\t%s\n' % (pdbid, chain.getChid(), extract_nuclotide_seq(pdbid, chain.getChid())) for chain in pdb.iterChains()]
    return result


def extract_nuclotide_all_seqs_frm_directory(directory=PDB_DB):
    '''This method is used to get all sequences of each pdb file in a folder.

    Args:
        directory: PDBs folder as default, given the directory where all pdb files are.

    Returns:
        No, write results into file.
    '''

    pdbs = [f.split('.')[0] for f in os.listdir(directory)]
    pdbs.sort()
    # return list(itertools.chain.from_iterable(map(extract_nuclotide_all_seqs_frm_one, pdbs)))
    with open('all_pdb_seqs.txt', 'w') as wf:
        for pdbid in pdbs:
            print pdbid
            # result.extend(extract_nuclotide_all_seqs_frm_one(pdbid))
            for line in extract_nuclotide_all_seqs_frm_one(pdbid):
                wf.write(line)


def extract_nuclotide_all_seqs_frm_directory_debug(directory=PDB_DB):
    '''This method is for debugging, listing distinct nucleotide letters of each sequence.

    Args:
        directory: PDBs folder as default, given the directory where all pdb files are.

    Returns:
        No, write results into file.

    Raises:
        No, only print error info
    '''

    try:
        pdbs = [f.split('.')[0] for f in os.listdir(directory)]
        pdbs.sort()
        with open('all_pdb_seqs_debug.txt', 'w') as wf:
            for pdbid in pdbs:
                print pdbid
                pdb = parsePDB(pdbid)
                for chain in pdb.iterChains():
                    s, l = extract_nuclotide_seq(pdbid, chain.getChid(), pdb, 1)
                    # wf.write('%s\t%s\t%s\n' % (pdbid, chain.getChid(), s))
                    wf.write('%s\t%s\t%s\n' % (pdbid, chain.getChid(), ','.join(set(l))))
    except Exception, e:
        print 'extract_nuclotide_all_seqs_frm_directory_debugV >>>', e


def extract_c1_atom_coordinates(pdbid, chainid, s, e, model=1):
    '''This method is extracting c1 coordinates for each nucleotide

    Args:
        pdbid: Case insensitive
        chainid: Case sensitive
        s: Start position of sequence
        e: End position of sequence

    Returns:
        A substring sequence with given indices.

    Raises:
        No, only print error info
    '''

    try:
        pdb = parsePDB(pdbid, chain=chainid, model=model)
        if not pdb: return
        c1s = pdb.select('name C1\'')
        indices = c1s.getResnums()
        # indices = indices[s:e]
        if s == 0 and e == 0:
            ## full sequence
            # return seq
            pass
        else:
            ## either s or e is 0
            if s == 0:
                ## start 0 to end index
                # return seq[:e]
                indices = indices[:e]
            elif e == 0:
                ## start index to end
                # return seq[s-1:]
                indices = indices[s-1:]
            else:
                indices = indices[s-1:e]

        for index in indices:
            a = c1s.select('resnum `%s`' % (index,)).getCoords()
            if a.size == 0:
                print 'NA\tNA\tNA\n'
            else:
                print '\n'.join('\t'.join('%0.4f' %x for x in y) for y in a)


    except Exception, e:
        print 'extract_nuclotide_seq >>>', e


def init():
    '''Prody package setting, ignore debugging info from package.'''
    confProDy(verbosity='none')
    pathPDBFolder(PDB_DB)


def main():
    '''
    Change to command parse later
    '''
    try:
        if len(sys.argv) == 5:
            '''Providing all arguments, pdbid, chainid, index start number (from 1), index end.

            If any index number is given 0, it will extract the whole sequence.
            '''
            pdbid = sys.argv[1].lower()
            # upper case and lower case could be in the same pdb file
            chainid = sys.argv[2]
            start = int(sys.argv[3])
            end = int(sys.argv[4])
            print extract_nuclotide_seq_indices(pdbid, chainid, start, end)
            extract_c1_atom_coordinates(pdbid, chainid, start, end)

        elif len(sys.argv) == 3:
            '''Providing only pdb ID and chain ID.

            Indices will be set 0 as default, extracting the whole sequence.
            '''
            pdbid = sys.argv[1].lower()
            chainid = sys.argv[2]
            print extract_nuclotide_seq_indices(pdbid, chainid, 0, 0)
            extract_c1_atom_coordinates(pdbid, chainid, 0, 0)

        elif len(sys.argv) == 2:
            '''Providing only pdb ID, extracting the whole sequence for all chains. '''
            pdbid = sys.argv[1].lower()
            for line in extract_nuclotide_all_seqs_frm_one(pdbid):
                print line.rstrip()

        elif len(sys.argv) == 1:
            '''Extracting all the whole sequences of every pdb in PDBs folder. '''
            # for line in extract_nuclotide_all_seqs_frm_directory('/home/common/PD'):
            #     print line.rstrip()

            # extract_nuclotide_all_seqs_frm_directory_debug()
            extract_nuclotide_all_seqs_frm_directory()

            # with open('all_pdb_seqs.txt', 'w') as wf:
            #     # map(wf.write, extract_nuclotide_all_seqs_frm_directory())
            #     # map(wf.write, extract_nuclotide_all_seqs_frm_directory('/home/common/PD'))
            #     for line in extract_nuclotide_all_seqs_frm_directory():
            #         print line
            #         wf.write(line)

        else:
            print '''arguments invalid!
            \'pdbid\'\t\t\t\t\textracting all chain's whole sequence.
            \'pdbid chainid\'\t\t\t\textracting the whole sequence.
            \'pdbid chainid index_start index_end\'\textracting a partition sequence. '''
            sys.exit(1)
    except Exception, e:
        if DEBUG:
            print 'e >>> ', e


if __name__ == '__main__':
    init()
    main()
