# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:03:20 2018

@author: Administrator
liuze
"""
from sklearn import preprocessing
import numpy as np
import csv
def file_to_features(filepath):
    bag_dic={}
    info_dic={}
    sample_instances=[]
    j=0
    #filepath='features.txt'
    filein=open(filepath,'r')
    for line in filein:
        head_info=line.split(',')[0]
        bag_key=int(line.split(',')[1])
        sample_instance=','.join(line.split(',')[2:])
        #sample_instance=line.split(',')[1:]
        #bag_dic[bag_key].append(sample_instance)
        if bag_key in bag_dic:
            bag_dic[bag_key].append(sample_instance)
            info_dic[bag_key].append(head_info)
        else:
            bag_dic[bag_key]=[]
            bag_dic[bag_key].append(sample_instance)
            info_dic[bag_key]=[]
            info_dic[bag_key].append(head_info)
    bag_index=bag_dic.keys()
    sample_list=bag_dic.values()
    info_index=info_dic.keys()
    info_list=info_dic.values()
    bag_instances=[[] for i in range(0,len(sample_list))]
    label_in_bag=[[] for i in range(0,len(sample_list))]
    np.savetxt("bag_index",bag_index, fmt='%s',delimiter='\n')
    
    sample_info=[]
    for info in info_list:
        for info_sample in info:
            sample_info.append(info_sample)
    
    for bag in sample_list:
        for instance in bag:
    #for i in range(0,len(sample_list)):
     #   for j in range(0,len(sample_list[i])):
            #sample_in_bag=str(sample_list[j]).split(',')
            sample_in_bag=str(instance).split(',')
            label_in_bag[j].append(int(str(sample_in_bag.pop()).strip()))
            #bag_instances[j].append(sample_in_bag)
            feature_in_bag=[]
            for i in range(0,len(sample_in_bag)):
                instFeature=float(sample_in_bag[i])
                feature_in_bag.append(instFeature)
            #bag_instances=np.append(bag_instances[j],feature_in_bag)
            bag_instances[j].append(np.array(feature_in_bag))
        j=j+1
    bag_targets=[]
    for this_bag_instance_targets in label_in_bag:
        if 1 in this_bag_instance_targets:
            bag_targets.append(1)
        else:
            bag_targets.append(-1)
    for i in range(0,len(bag_instances)):
        #bag_instances[i]=np.array(bag_instances[i], dtype=np.float32)
        bag_instances[i]=bag_instances[i]
    return [bag_targets, bag_instances, label_in_bag, sample_info]
 
#[instance_bag_ids, instance_targets, instance_samples]=file_to_features('features.txt')
#np.savetxt("bag_name.txt",instance_bag_ids, fmt='%s',delimiter=',')
#np.savetxt("bag_instance_targets.txt",instance_targets, fmt='%s',delimiter=',')
#np.savetxt("bag_instance_samples.txt",instance_samples, fmt='%s',delimiter=',')
    
#with open('dict.csv', 'wb') as csv_file:
#    writer = csv.writer(csv_file)
#    for key, value in bag_dic.items():
#       writer.writerow([key, value])
#filein.close() 
 
    