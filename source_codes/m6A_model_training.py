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
Usage: python {0} positive_dataset negative_dataset method model_file scale_file score_file reserved_file

This script trains a miSVM/MISVM model.

method : method of selection
     0 -- miSVM
     1 -- MISVM
A classifier model file will be created in current path.""".format(argv[0]))
	exit(1)

def process_options(argv):
    argc=len(argv)
    if argc!=8:
        exit_with_help(argv)
    posCV=open(argv[1],'r')
    negCV=open(argv[2],'r')
    method=int(argv[3])
    model=argv[4]
    scale=argv[5]
    score=open(argv[6],'w+')
    res=open(argv[7],'w+')
    if method not in [0,1]:
        print ("Unknown selection method {0}".format(method))
        exit_with_help(argv)
    return posCV, negCV, method, model, scale, score, res
        
def main(argv=sys.argv):
    pos_train_file, neg_train_file, method, model_file, scale_file, score_file, res_file=process_options(argv)
######train feature extraction
    feature_matrix=[]
    for line in pos_train_file:
        feature_vector=[]
        sequence_infor=line.split(',')
        header=sequence_infor[0]
        bag=sequence_infor[1]
        sequence=sequence_infor[3].strip()
        feature_vector.append(header)
        feature_vector.append(bag)
        feature_vector.extend(kmer(sequence)+ksnpf(sequence)+nucleic_shift(sequence))
        #feature_vector.extend(ksnpf(sequence))
        feature_vector.append('1')
        feature_matrix.append(feature_vector)
    pos_train_file.close()
    
    for line in neg_train_file:
        feature_vector=[]
        sequence_infor=line.split(',')
        header=sequence_infor[0]
        bag=sequence_infor[1]
        sequence=sequence_infor[3].strip()
        feature_vector.append(header)
        feature_vector.append(bag)
        feature_vector.extend(kmer(sequence)+ksnpf(sequence)+nucleic_shift(sequence))
        #feature_vector.extend(ksnpf(sequence))
        feature_vector.append('-1')
        feature_matrix.append(feature_vector)
    feature_array = np.array([b[2:-1] for b in feature_matrix],dtype=np.float32)
    min_max_scaler = preprocessing.MinMaxScaler(copy=True, feature_range=(-1, 1))
    feature_scaled= min_max_scaler.fit_transform(feature_array)
    
    feature_matrix_T=map(list,zip(*feature_matrix))
    feature_scaled_T=map(list,zip(*feature_scaled))
    k=0
    train_feature_matrix_T=[]
    train_feature_matrix_T.append(feature_matrix_T[0])
    train_feature_matrix_T.append(feature_matrix_T[1])
    for i in range(len(feature_scaled_T)):
        train_feature_matrix_T.append(feature_scaled_T[k])
        k=k+1
    train_feature_matrix_T.append(feature_matrix_T[-1])
    train_feature_matrix=map(list,zip(*train_feature_matrix_T))
    neg_train_file.close()
    np.savetxt("train_features.txt",train_feature_matrix, fmt='%s',delimiter=',')    
######put samples into bags
    train_file_path='./train_features.txt'
    [train_bag_targets,train_bag_samples,train_bag_instance_targets,sample_info]=file_to_features(train_file_path)    
    if method==0:
        svc=misvm.miSVM(kernel='quadratic',C=5.4,max_iters=10)
    elif method==1:
        svc=misvm.MISVM(kernel='quadratic',C=5.4,max_iters=10)
    svc.fit(train_bag_samples,train_bag_targets)
    #joblib.dump(svc,'./svc.pkl')
    joblib.dump(svc,model_file)
    joblib.dump(min_max_scaler,scale_file)
    bag_predictions,inst_predictions=svc.predict(train_bag_samples,1)
    #score_file=open('training_score.txt','w+')
    score_file.write("name\tscore\n")
    for i in range(0,len(sample_info)):
        score_file.write(sample_info[i]+"\t"+str(inst_predictions[i])+"\n")
    score_file.close()
    #reserved_samples=open('reserved_samples.txt','w+')
    res_file.write("name\n")
    for j in range(0,len(inst_predictions)):
        if inst_predictions[j]>=0:
            res_file.write(sample_info[j]+"\n")
    res_file.close()
    
if __name__=='__main__':
    main(sys.argv)



