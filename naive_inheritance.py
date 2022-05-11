#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:32:53 2022

@author: pazhao
"""
#adapted from Eichler Lab's existing naive_inheritance_trio.py
from __future__ import print_function
import sys
import gzip
import argparse
import pandas as pd
import numpy as np
parser = argparse.ArgumentParser(description='adding denovo or inheritance information as an extra info field to vcf output from variant callers such as freebayes, GATK')
parser.add_argument('-i', action = 'store', dest = 'inputfile', required = True, type=str, help='input vcf file that has been reheaded, i.e. sample id to family id')
parser.add_argument('-o',action='store',dest='output',required=True,type=str,help='')
parser.add_argument('-m', action = 'store', dest = 'fmManifest', required = True, type=str, help='manifest file for the family,which provides col indices for fa, mo, probands, sex of probands, sibling and sibling sex if applicable')

#father_index, mother_index, prob_index, list_of_prob_sex, sib_index=None,list_of_sib_sex=None)
args = parser.parse_args()
infile=str(args.inputfile)
outfile=str(args.output)
fmManifest=pd.read_table(args.fmManifest,header=None)
father_index=int(np.where(fmManifest.iloc[:,1].str.contains('\.fa'))[0][0]) + 9
mother_index=int(np.where(fmManifest.iloc[:,1].str.contains('\.mo'))[0][0]) + 9
proband_index=[ i  + 9  for i in np.where(fmManifest.iloc[:,1].str.contains('\.p'))[0].tolist() ]
if np.where(fmManifest.iloc[:,1].str.contains('\.s'))[0].size > 0:
    sib_index=[ i  + 9 for i in np.where(fmManifest.iloc[:,1].str.contains('s'))[0].tolist()]
    list_of_sib_sex=[ fmManifest.iloc[i-9,1] . split('.')[1] for i in sib_index]
else:
    sib_index=None
    list_of_sib_sex=None
list_of_prob_sex=[fmManifest.iloc[i-9,1] . split('.')[1] for i in proband_index]
print (list_of_prob_sex)

# INFO header lines added
INFO_HEADER_LINES = [
    '##INFO=<ID=INH,Number=A,Type=String,'
    'Description="inheritance pattern treating autosome and sex chromosome separately">\n'
]

class VariantReader(object):
        def __init__(self, filename):
                self.fh = gzip.GzipFile(filename, 'r')
        def __iter__(self):
                return self
        def next(self):
                while True:
                        line = self.fh.readline()
                        if line.startswith( '#' ):
                                continue
                        if line == "":
                                self.fh.close()
                                raise StopIteration
                        line = line[:-1]
                        return Variant(line.split('\t'))
def loop_one_parent(list_of_offspring_entries, offspring_name='pro',het='fa'):
    #this function is to query offspring when the parents are (0/0 ,0/1) or (0/0 ,1/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0':
            prob_info_str=het+'_to_'+offspring_name+str(i+1)+'_'
        elif offspring.split(":")[0] == '0/0':
            prob_info_str=''
        else:
            prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
        if i==0:
            full_prob_info_str=prob_info_str
        else:
            full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)
def loop_two_parents(list_of_offspring_entries, offspring_name='pro'):
    #this function is to query offspring when the parents are (0/1 ,0/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0' or offspring.split(":")[0] == '1/1' :
            prob_info_str='mo_fa_to_'+offspring_name+str(i+1)+'_'
        elif offspring.split(":")[0] == '0/0':
            prob_info_str=''
        else  :
            prob_info_str=offspring_name+str(i+1)+'_NA_'
        if i==0:
            full_prob_info_str=prob_info_str
        else:
            full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)
def loop_two_parents_one_hom(list_of_offspring_entries, offspring_name='pro',hom='fa'):
    #this function is to query offspring when the parents are (0/1 ,1/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0' :
            prob_info_str=hom+'_to_'+offspring_name+str(i+1)+'_'
        elif offspring.split(":")[0] == '0/0':
            prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
        elif offspring.split(":")[0] == '1/1' :
            prob_info_str='fa_mo_to_'+offspring_name+str(i+1)+'_'
        else:
            prob_info_str=offspring_name+str(i+1)+'_NA_'
        if i==0:
            full_prob_info_str=prob_info_str
        else:
            full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)
def loop_two_parents_two_hom(list_of_offspring_entries, offspring_name='pro'):
    #this function is to query offspring when the parents are  (1/1,1/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if '0' in offspring.split(":")[0]  :
            prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
        elif offspring.split(":")[0] == '1/1' :
            prob_info_str='mo_fa_to_'+offspring_name+str(i+1)+'_'
        else:
            prob_info_str=offspring_name+str(i+1)+'_NA_'
        if i==0:
            full_prob_info_str=prob_info_str
        else:
            full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)

def loop_homreffa(list_of_offspring_entries,list_of_offspring_sex, offspring_name='pro'):
    #this function is to query offspring when the parents are (0/0 ,0/1) or (0/0 ,1/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if list_of_offspring_sex[i]=='F':
            if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0':
                prob_info_str='mo_to_'+offspring_name+str(i+1)+'_'
            elif offspring.split(":")[0] == '0/0':
                prob_info_str=''
            else:
                prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
            if i==0:
                full_prob_info_str=prob_info_str
            else:
                full_prob_info_str=full_prob_info_str+prob_info_str
        elif list_of_offspring_sex[i]=='M':
            if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0':
                prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
            elif offspring.split(":")[0] == '0/0':
                prob_info_str=''
            elif offspring.split(":")[0] == '1/1':
                prob_info_str='mo_to_'+offspring_name+str(i+1)+'_'
            else:
                prob_info_str=offspring_name+str(i+1)+'_NA_'
            if i==0:
                full_prob_info_str=prob_info_str
            else:
                full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)

def loop_hetmo_homvarfa(list_of_offspring_entries, list_of_offspring_sex, offspring_name='pro',hom='fa'):
    #this function is to query chrX offspring when the parents are (0/1 ,1/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if list_of_offspring_sex[i]=='F':
            if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0' :
                prob_info_str=hom+'_to_'+offspring_name+str(i+1)+'_'
            elif offspring.split(":")[0] == '0/0':
                prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
            elif offspring.split(":")[0] == '1/1' :
                prob_info_str='fa_mo_to_'+offspring_name+str(i+1)+'_'
            else:
                prob_info_str=offspring_name+str(i+1)+'_NA_'
            if i==0:
                full_prob_info_str=prob_info_str
            else:
                full_prob_info_str=full_prob_info_str+prob_info_str
        elif  list_of_offspring_sex[i]=='M':
            if offspring.split(":")[0] == '0/1'or offspring.split(":")[0] == '1/0' :
                prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
            elif offspring.split(":")[0] == '0/0':
                prob_info_str=''
            elif offspring.split(":")[0] == '1/1':
                prob_info_str='mo_to_'+offspring_name+str(i+1)+'_'
            else:
                prob_info_str=offspring_name+str(i+1)+'_NA_'
            if i==0:
                full_prob_info_str=prob_info_str
            else:
                full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)

def loop_homvarmo_homvarfa(list_of_offspring_entries, list_of_offspring_sex,offspring_name='pro'):
    #this function is to query chrX offspring when the parents are (0/1 ,1/1)
    for i,offspring in enumerate(list_of_offspring_entries):
        if list_of_offspring_sex[i]=='F':
            if offspring.split(":")[0] == '1/1' :
                prob_info_str='fa_mo_to_'+offspring_name+str(i+1)+'_'
            elif '0' in offspring.split(":")[0]  :
                prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
            else:
                prob_info_str=offspring_name+str(i+1)+'_NA_'
            if i==0:
                full_prob_info_str=prob_info_str
            else:
                full_prob_info_str=full_prob_info_str+prob_info_str
        elif  list_of_offspring_sex[i]=='M':
            if offspring.split(":")[0] == '1/1':
                prob_info_str='mo_to_'+offspring_name+str(i+1)+'_'
            elif '0' in offspring.split(":")[0]  :
                prob_info_str=offspring_name+str(i+1)+'unMendelian'+'_'
            else:
                prob_info_str=offspring_name+str(i+1)+'_NA_'
            if i==0:
                full_prob_info_str=prob_info_str
            else:
                full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)

def  loop_chrYvariants(list_of_offspring_entries,list_of_offspring_sex,offspring_name='pro'):
    for i,offspring in enumerate(list_of_offspring_entries):
        if list_of_offspring_sex[i] =='M':
            if '1' in offspring.split(":")[0]  :
                prob_info_str='fa_to_'+offspring_name+str(i+1)+'_'
            else:
                prob_info_str=''
        else:
            prob_info_str=''
        if i==0:
            full_prob_info_str=prob_info_str
        else:
            full_prob_info_str=full_prob_info_str+prob_info_str
    return(full_prob_info_str)


class Variant(object):
        def __init__(self, row,father_index, mother_index, prob_index, list_of_prob_sex, sib_index=None,list_of_sib_sex=None):
                self.chrom = row[0]
                self.end = row[1]
                self.id = row[2]
                self.ref = row[3]
                self.alt = row[4]
                self.qual = row[5]
                self.filter = row[6]
                self.info = row[7]
                self.format = row[8]
                self.father = row[father_index]
                self.mother = row[mother_index]
                self.proband=[]
                self.proband_sex=[]
                for j,i in enumerate(prob_index):
                    self.proband.append(row[i])
                    self.proband_sex.append(list_of_prob_sex[j])
                if sib_index:
                    self.sibling = []
                    self.sibling_sex=[]
                    for j,i in enumerate(sib_index):
                        self.sibling.append (row[i])
                        self.sibling_sex.append(list_of_sib_sex[j])
#                if not self.chrom.startswith("chr"):
#                        self.chrom = 'chr' + self.chrom
                else:
                    full_sibling_info_str = ''
                if self.chrom!='X' and self.chrom!='Y':
                    #treat sex chromosomes and autosomes separately
                    if self.father.split(':')[0] == '0/0' and (self.mother.split(':')[0] == '0/1' or self.mother.split(':')[0] == '1/1'or self.mother.split(':')[0] == '1/0' ) :
                        #get mo_transmission, one het, one homref
                        full_prob_info_str=loop_one_parent(self.proband,het='mo')
                        if sib_index:
                            full_sibling_info_str = loop_one_parent(self.sibling,offspring_name='sib',het='mo')
                        self.info=self.info + ';INH='+    full_prob_info_str +      full_sibling_info_str
                    elif self.mother.split(':')[0] == '0/0' and (self.father.split(':')[0] == '0/1' or self.father.split(':')[0] == '1/1'or self.father.split(':')[0] == '1/0' ) :
                        #get_fa_transmission one het, one homref
                        full_prob_info_str=loop_one_parent(self.proband,het='fa')
                        if sib_index:
                            full_sibling_info_str=loop_one_parent(self.sibling,offspring_name='sib',het='fa')
                        self.info=self.info + ';INH='+     full_prob_info_str +      full_sibling_info_str
                    elif (self.mother.split(':')[0] == '0/1' or self.mother.split(':')[0] == '1/0' )and ( self.father.split(':')[0] == '0/1' or  self.father.split(':')[0] == '1/0'):
                        #get transmission, albeit unknown origin, both parents het
                        full_prob_info_str=loop_two_parents(self.proband)
                        if sib_index:
                            full_sibling_info_str=loop_two_parents(self.sibling,offspring_name='sib')
                        self.info=self.info + ';INH='+     full_prob_info_str +      full_sibling_info_str
                    elif (self.mother.split(':')[0] == '0/1' or self.mother.split(':')[0] == '1/0' ) and self.father.split(':')[0] == '1/1':
                        #get transmission  one het, one homvar
                         full_prob_info_str=loop_two_parents_one_hom(self.proband)
                         if sib_index:
                             full_sibling_info_str=loop_two_parents_one_hom(self.sibling,offspring_name='sib')
                         self.info=self.info + ';INH='+  full_prob_info_str +      full_sibling_info_str
                    elif (self.father.split(':')[0] == '0/1' or self.father.split(':')[0] == '1/0' ) and self.mother.split(':')[0] == '1/1':
                        ##get transmission , one het, one homvar
                         full_prob_info_str=loop_two_parents_one_hom(self.proband, hom='mo')
                         if sib_index:
                             full_sibling_info_str=loop_two_parents_one_hom(self.sibling,offspring_name='sib',hom='mo')
                         self.info=self.info + ';INH='+  full_prob_info_str +      full_sibling_info_str
                    elif self.father.split(':')[0] == '1/1' and self.mother.split(':')[0] == '1/1' :
                        ##get transmission , two homvar
                         full_prob_info_str=loop_two_parents_two_hom(self.proband)
                         if sib_index:
                             full_sibling_info_str=loop_two_parents_two_hom(self.sibling,offspring_name='sib')
                         self.info=self.info + ';INH='+  full_prob_info_str +      full_sibling_info_str
                    #get denovo
                    elif self.father.split(':')[0] == '0/0' and self.mother.split(':')[0] == '0/0' :
                        for i,proband in enumerate(self.proband):
                            if '1' in proband.split(":")[0] :
                                prob_info_str='pro'+str(i+1)+'_'
                            else:
                                prob_info_str=''
                            if i==0:
                                full_prob_info_str=prob_info_str
                            else:
                                full_prob_info_str=full_prob_info_str+prob_info_str
                        if sib_index:
                            for i, sibling in enumerate(self.sibling):
                                if '1' in sibling.split(":")[0] :
                                    sibling_info_str='sib'+str(i+1)+'_'
                                else:
                                    sibling_info_str=''
                                if i==0:
                                    full_sibling_info_str=sibling_info_str
                                else:
                                    full_sibling_info_str=full_sibling_info_str+sibling_info_str

                        self.info=self.info + ';INH=denovo_'+     full_prob_info_str +      full_sibling_info_str
                    else:
                        self.info=self.info+';INH=NA'
                elif self.chrom=='X' :
                    if  self.father.split(':')[0] == '0/0' and (self.mother.split(':')[0] == '0/1' or self.mother.split(':')[0] == '1/1'or self.mother.split(':')[0] == '1/0' ) :
                        full_prob_info_str=loop_homreffa(self.proband,self.proband_sex)
                        if sib_index:
                            full_sibling_info_str=loop_homreffa(self.sibling,self.sibling_sex,offspring_name='sib')
                        self.info=self.info + ';INH='+     full_prob_info_str +      full_sibling_info_str
                    elif (self.mother.split(':')[0] == '0/1' or self.mother.split(':')[0] == '1/0' ) and self.father.split(':')[0] == '1/1' :
                        #get transmission , albeit unknown origin, one het, one homvar
                         full_prob_info_str=loop_hetmo_homvarfa(self.proband,self.proband_sex)
                         if sib_index:
                             full_sibling_info_str=loop_hetmo_homvarfa(self.sibling,self.sibling_sex,offspring_name='sib')
                         self.info=self.info + ';INH='+  full_prob_info_str +      full_sibling_info_str
                    elif self.father.split(':')[0] == '1/1' and self.mother.split(':')[0] == '1/1' :
                        full_prob_info_str=loop_homvarmo_homvarfa(self.proband, self.proband_sex)
                        if sib_index:
                            full_sibling_info_str=loop_homvarmo_homvarfa(self.sibling, self.sibling_sex,offspring_name='sib')
                        self.info=self.info + ';INH=' + full_prob_info_str +      full_sibling_info_str
                    #get denovo
                    elif self.father.split(':')[0] == '0/0' and self.mother.split(':')[0] == '0/0' :
                        for i,proband in enumerate(self.proband):
                            if '1' in proband.split(":")[0] :
                                prob_info_str='pro'+str(i+1)+'_'
                            else:
                                prob_info_str=''
                            if i==0:
                                full_prob_info_str=prob_info_str
                            else:
                                full_prob_info_str=full_prob_info_str+prob_info_str
                        if sib_index:
                            for i, sibling in enumerate(self.sibling):
                                if '1' in sibling.split(":")[0] :
                                    sibling_info_str='sib'+str(i+1)+'_'
                                else:
                                    sibling_info_str=''
                                if i==0:
                                    full_sibling_info_str=sibling_info_str
                                else:
                                    full_sibling_info_str=full_sibling_info_str+sibling_info_str

                            self.info=self.info + ';INH=denovo_'+     full_prob_info_str +      full_sibling_info_str
                    else:
                        self.info=self.info+';INH=NA'
                elif self.chrom=='Y':
                    if '1' in self.father.split(':')[0]:
                        full_prob_info_str=loop_chrYvariants(self.proband, self.proband_sex)
                        if sib_index:
                            full_sibling_info_str=loop_chrYvariants(self.sibling, self.sibling_sex,offspring_name='sib')
                        self.info=self.info+';INH='+ full_prob_info_str +      full_sibling_info_str
                        #get denovo
                    elif self.father.split(':')[0] == '0/0'  :
                        for i,proband in enumerate(self.proband):
                            if '1' in proband.split(":")[0] and self.proband_sex[i] =='M' :
                                prob_info_str='pro'+str(i+1)+'_'
                            else:
                                prob_info_str=''
                            if i==0:
                                full_prob_info_str=prob_info_str
                            else:
                                full_prob_info_str=full_prob_info_str+prob_info_str
                        if sib_index:
                            for i, sibling in enumerate(self.sibling):
                                if '1' in sibling.split(":")[0] and self.sibling_sex[i] =='M' :
                                    sibling_info_str='sib'+str(i+1)+'_'
                                else:
                                    sibling_info_str=''
                                if i==0:
                                    full_sibling_info_str=sibling_info_str
                                else:
                                    full_sibling_info_str=full_sibling_info_str+sibling_info_str
                        if full_prob_info_str!='' or full_sibling_info_str !='':
                            self.info=self.info + ';INH=denovo_'+     full_prob_info_str +      full_sibling_info_str
                    else:
                        self.info=self.info+';INH=NA'
        def get_sib(self):
            if hasattr(self,'sibling'):
                return(True)
            else:
                return(False)
        def __repr__(self):
            if Variant.get_sib(self):
                return "\t".join([self.chrom, str(self.end), self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format, self.father, self.mother] + [i for i in  self.proband] +  [i for i in self.sibling])
            else:
                return "\t".join([self.chrom, str(self.end), self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format, self.father, self.mother] + [i for i in  self.proband] )

def is_gzip_file_name(file_name):
    """
    Determine if a filename refers to a gzipped file.

    :param file_name: File name.

    :return: `True` if `file_name` is a gzipped file. `False` if the file name is `None` or does not end with ".gz".
    """

    if file_name is None:
        return False

    file_name = file_name.strip()

    return file_name.endswith('.gz')
class ContextFile:
    """
    Open a file as a regular file or a gzipped file depending on the file name. This class is a context-guard allowing
    the file to be opened in "with" statements.
    """

    def __init__(self, file_name, mode):
        """
        Create a new context-guarded file. Automatically opens as GZIP or not depending on the file extension.

        :param file_name: File name to open.
        :param mode: "r" (read) or "w" (write).
        """

        # Check arguments
        if file_name is not None:
            file_name = file_name.strip()

        if not file_name:  # No file name
            file_name = None

        if mode not in ('r', 'w'):
            raise RuntimeError('File mode must be "r" or "w": {}'.format(mode))

        # Set fields
        self.file_name = file_name
        self.mode = mode

        # Set state
        self.file = None
    def __enter__(self):
        """
        Open the file.

        :return: Open file handle.
        """

        # Get default file if no file name
        if self.file_name is None:
            if self.mode == 'r':
                return sys.stdin
            else:
                return sys.stdout

        # Open file
        if self.file_name.endswith('.gz'):
            self.file = gzip.open(self.file_name, self.mode)
        else:
            self.file = open(self.file_name, self.mode)

        return self.file

    def __exit__(self, type, value, traceback):
        """
        Close the open file.
        """

        if self.file is None:
            return  # Do not close stdin or stdout

        self.file.close()
        self.file = None

def write_info_header_lines(out_file, is_gz):
    """
    Write VCF header lines added by this script.

    :param out_file: Open output file.
    :param is_gz: Lines are encoded to bytes before writing if `True` (for gzip output).
    """

    for line in INFO_HEADER_LINES:
        if is_gz:
            out_file.write(line.encode())
        else:
            out_file.write(line)
if __name__ == '__main__':
    out_is_gz = is_gzip_file_name(outfile)
    vcfHeaderList = list()
with ContextFile(outfile, 'w') as outVcfFile:
    with ContextFile(infile, 'r') as inVcfFile:
        formatFound = False
        line = ''
        for line in inVcfFile:
            if isinstance(line, bytes):
                line = line.decode('utf-8')
            if not line:
                continue
            if line.startswith ("#CHROM"):
                if sib_index:
                    line=str('\t'.join([ '#CHROM' ,'POS',     'ID',      'REF',     'ALT',     'QUAL',    'FILTER'  ,'INFO',    'FORMAT','fa','mo'] +['p'+str(i) for i in range(1,len(proband_index)+1)] + ['s'+str(i)  for i in range(1,len(sib_index)+1)])) + '\n'
                else:
                    line=str('\t'.join([ '#CHROM' ,'POS',     'ID',      'REF',     'ALT',     'QUAL',    'FILTER'  ,'INFO',    'FORMAT','fa','mo'] +['p'+str(i) for i in range(1,len(proband_index)+1)] )) + '\n'
            if not line.startswith('#'):
                break
            if line.startswith('##FORMAT') and not formatFound:
                write_info_header_lines(outVcfFile,out_is_gz)
                formatFound=True
            if not line.startswith("##"):
                if not formatFound:
                    write_info_header_lines(outVcfFile, out_is_gz)
            if out_is_gz:
                line=line.encode()
            outVcfFile.write(line)
        while True:
            # Prepare line
            if isinstance(line, bytes):  # Line may be string or bytes depending on gzip file or not
                line = line.decode('utf-8')
            line = line.strip('\n')
            line=line.split('\t')
            if not line:
                break
            if sib_index:
                record_line = str(Variant(line, father_index, mother_index, proband_index,list_of_prob_sex,sib_index,list_of_sib_sex)) + '\n'
            else:
                record_line = str(Variant(line, father_index, mother_index, proband_index,list_of_prob_sex)) + '\n'
            if out_is_gz:
                record_line = record_line.encode()
            outVcfFile.write(record_line)
            try:
                line = next(inVcfFile)
            except StopIteration:
                break
