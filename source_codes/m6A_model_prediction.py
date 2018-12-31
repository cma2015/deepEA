#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 11:11:20 2018

@author: liuze
"""
import sys
from sklearn.externals import joblib
from feature_selection import *
from sklearn import preprocessing
import misvm
from create_bag_dict import file_to_features
from sklearn import metrics

def exit_with_help(argv):
	print("""\
Usage: python {0} dataset model_file scale_file predict_res

This script will generate scores to relabel the candidate m6As.""".format(argv[0]))
	exit(1)

def process_options(argv):
    argc=len(argv)
    if argc!=5:
        exit_with_help(argv)
    dataset=open(argv[1],'r')
    model_file=joblib.load(argv[2])
    scale_file=joblib.load(argv[3])
    predic_res=open(argv[4],'w+')
    return dataset, model_file, scale_file, predic_res
        
def main(argv=sys.argv):
    data_file, model_file, scale_file, predict_file=process_options(argv)
######test feature extraction
    feature_matrix=[]
    info_vector=[]
    for line in data_file:
        feature_vector=[]
        sequence_infor=line.split(',')
        header=sequence_infor[0]
        info_vector.append(header)
        sequence=sequence_infor[1].strip()
        feature_vector.append(header)
        feature_vector.extend(kmer(sequence)+ksnpf(sequence)+nucleic_shift(sequence))
        #feature_vector.extend(ksnpf(sequence))
        feature_vector.append('1')
        feature_matrix.append(feature_vector)
    data_file.close()
#        
    feature_vector = np.array([b[1:-1] for b in feature_matrix],dtype=np.float32)
    feature_scaled= scale_file.transform(feature_vector) 
    bag_predictions,inst_predictions=model_file.predict(feature_scaled,1)
    
    predict_file.write("name\tscore\tlabel\n")
    for i in range(0,len(inst_predictions)):
        if inst_predictions[i]>=0:
            predict_file.write(info_vector[i]+"\t"+str(inst_predictions[i])+"\t"+"1\n")
        else:
            predict_file.write(info_vector[i]+"\t"+str(inst_predictions[i])+"\t"+"-1\n")
    predict_file.close()
if __name__=='__main__':
    main(sys.argv)