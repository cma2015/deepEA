# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 15:35:41 2018

@author: Administrator
liuze
"""
from __future__ import division
import sys
from functools import reduce
import re
import operator
from math import log
import numpy as np
def kmer(seq):
    mer2={}
    mer3={}
    mer4={}
    for n1 in 'ATCG':
        for n2 in 'ATCG':
            mer2[n1+n2]=0
            for n3 in 'ATCG':
                mer3[n1+n2+n3]=0
                for n4 in 'ATCG':
                    mer4[n1+n2+n3+n4]=0
    seq=seq.replace('N','')
    seq_len=len(seq)
    for p in range(0,seq_len-3):
        mer2[seq[p:p+2]]+=1
        mer3[seq[p:p+3]]+=1
        mer4[seq[p:p+4]]+=1
    mer2[seq[p+1:p+3]]+=1
    mer2[seq[p+2:p+4]]+=1
    mer3[seq[p+1:p+4]]+=1
    v2=[]
    v3=[]
    v4=[]
    for n1 in 'ACGT':
        for n2 in 'ACGT':
            v2.append(mer2[n1+n2])
            for n3 in 'ACGT':
                v3.append(mer3[n1+n2+n3])
                for n4 in 'ACGT':
                    v4.append(mer4[n1+n2+n3+n4])
    v=v2+v3+v4
    return v

def ksnpf(seq):
    kn=5
    freq=[]
    v=[]
    for i in range(0,kn):
        freq.append({})
        for n1 in 'ATCGN':
            freq[i][n1]={}
            for n2 in 'ATCGN':
                freq[i][n1][n2]=0
    seq=seq.strip('N')
    seq_len=len(seq)
    for k in range(0,kn):
        for i in range(seq_len-k-1):
            n1=seq[i]
            n2=seq[i+k+1]
            freq[k][n1][n2]+=1
    for i in range(0,kn):
        for n1 in 'ATCG':
            for n2 in 'ATCG':
                v.append(freq[i][n1][n2])
    return v
def binary_code(seq):
    binary_dictionary={'A':[1,1,1],'T':[0,1,0],'G':[1,0,0],'C':[0,0,1],'N':[0,0,0]}
    cnt=[]
    for i in seq:
        cnt.append(binary_dictionary[i])
    return reduce(operator.add,cnt)
def nucleic_shift(seq):
    nucleic_dictionary={'A':0,'T':0,'C':0,'G':0,'N':0}
    p=0
    v=[]
    for n in seq:
        p=p+1
        nucleic_dictionary[n]+=1
        v.append(nucleic_dictionary[n]/p)
    return v

#v1=nucleic_shift('ACGTTTT')
#v=ksnpf('ACGTACGT')                    