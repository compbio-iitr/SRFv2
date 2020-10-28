#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 01:39:51 2020

@author: backpropagator
"""

import json
import numpy as np
import math

def matching_score(seq1,seq2):
    length = len(seq1)
    score = 0.0
    for i in range(length):
        if seq1[i] == seq2[i]:
            score = score + 1
    
    percent = score/length
    return percent
    
def score_match(short_seq,long_seq):
    short_length = len(short_seq)
    long_length = len(long_seq)
    
    diff = long_length - short_length
    
    percent_score = 0
    
    seq1 = short_seq
    
    for i in range(diff + 1):
        seq2 = long_seq[i:short_length + i]
        percent_match = matching_score(seq1,seq2)
        if percent_match > percent_score:
            percent_score = percent_match
            
    return percent_score

#data_json = './data1.json'
#modified_json = './modified1.json'
#end = 2000

def modification_json(data_json,modified_json,end):

    if end <= 500:
        constraint = 0.2*end
    elif end <= 1000:
        constraint = 100
    else:
        constraint = 0.1*end
    modified_file = {}
    modified_data = {}
    unique_mer_indexes = {}
        
    with open(data_json) as f:
        file = json.load(f)
    
    m = file['mers']
    data = file['data']
    
    mers = []
    
    for mer_length in m:
        mers.append(int(mer_length))
    
    #mers = [23]
    mers = np.array(mers)
    mers.sort()
    #
    for mer in mers:
        data_mer = data[str(mer)]
        
        data_mer_length = len(data_mer)
        
        unique_mer_indexes[str(mer)] = []
    #    
        
        if data_mer_length == 1:
            data_mer[0]['modified'] = "true"
            unique_mer_indexes[str(mer)].append(0)
    #    
        if data_mer_length > 1:
            details = {}
            ncopies = []
            for i in range(data_mer_length):
                details[str(i)] = []
                data_mer[i]['modified'] = "false"
                ncopies.append(data_mer[i]['patterns'][0]['number_of_copies'])
                details[str(i)].append(data_mer[i]['from'])
                details[str(i)].append(data_mer[i]['upto'])
                details[str(i)].append(data_mer[i]['patterns'][0]['pattern_searched'])
                details[str(i)].append(data_mer[i]['patterns'][0]['consensus'])
            ncopies = np.array(ncopies)
            indexes = ncopies.argsort()[::-1]
            data_mer[indexes[0]]['modified'] = 'true'
            
            for i in range(indexes.shape[0]-1):
                index1 = indexes[i]
                index2 = indexes[i+1]
                start1 = details[str(index1)][0]
                end1 = details[str(index1)][1]
                pattern1 = details[str(index1)][2]
                consensus1 = details[str(index1)][3]
                start2 = details[str(index2)][0]
                end2 = details[str(index2)][1]
                pattern2 = details[str(index2)][2]
                consensus2 = details[str(index2)][3]
                
                pattern_score = matching_score(pattern1,pattern2)
                consensus_score = matching_score(consensus1,consensus2)
                
                if (abs(start1-start2) < constraint and abs(end1-end2) < constraint):
                    data_mer[index2]['modified'] = 'false'
                else:
                    data_mer[index2]['modified'] = 'true'
                    data_mer[index1]['modified'] = 'true'
            for i in range(data_mer_length):
                if data_mer[i]['modified'] == "true":
                    unique_mer_indexes[str(mer)].append(i)
        modified_data[str(mer)] = data_mer
    a = mers.shape[0] - 1
    ##
    while a >= 0:
        match = {}
        b = True
    #    
        mer_length = mers[a]
        match[str(mer_length)] = []
        unique_mer_index = unique_mer_indexes[str(mer_length)]
        data_mer = data[str(mer_length)]
        length = mer_length
        while b:
            c = a - 1
            if c >= 0:
                mer_length1 = mers[c]
                match[str(mer_length1)] = []
                unique_mer_index1 = unique_mer_indexes[str(mer_length1)]
                data_mer1 = data[str(mer_length1)]
                length1 = mer_length1
                for i in range(len(unique_mer_index)):
                    for j in range(len(unique_mer_index1)):
                        start = data[str(mer_length)][unique_mer_index[i]]['from']
                        end = data[str(mer_length)][unique_mer_index[i]]['upto']
                        consensus = data[str(mer_length)][unique_mer_index[i]]['patterns'][0]['consensus']
                        pattern = data[str(mer_length)][unique_mer_index[i]]['patterns'][0]['pattern_searched']
                        
                        start1 = data[str(mer_length1)][unique_mer_index1[j]]['from']
                        end1 = data[str(mer_length1)][unique_mer_index1[j]]['upto']
                        consensus1 = data[str(mer_length1)][unique_mer_index1[j]]['patterns'][0]['consensus']
                        pattern1 = data[str(mer_length1)][unique_mer_index1[j]]['patterns'][0]['pattern_searched']
                        
                        pattern_score = score_match(pattern1,pattern)
                        consensus_score = score_match(consensus1,consensus)
                        
                        if (consensus_score >= 0.9 or pattern_score >= 0.9) and (abs(start-start1) < constraint and abs(end-end1) < constraint):
                            match[str(mer_length)].append(unique_mer_index[i])
                            match[str(mer_length1)].append(unique_mer_index1[j])
                
                if len(match[str(mer_length1)]) == 0:
    #                print(match)
                    b = False
            else:
                b = False
    ##        
            a = a - 1
        numbers = 0
        for m in match.keys():
            match[m] = np.array(match[m])
            match[m] = np.unique(match[m])
            numbers = numbers + match[m].shape[0]
            
    ##    
        if numbers > 1:
            scores = []
            copies = []
            lengths = []
            copy = 0
            copy_mer = 0
            length = 0
            for m in match.keys():
                data_mer = data[m]
                for i in range(match[m].shape[0]):
                    copy_number = int(data_mer[match[m][i]]['patterns'][0]['number_of_copies'])
                    if copy_number > copy:
                        copy = copy_number
                        copy_mer = int(m)
                    
                    if int(m) > length:
                        length = int(m)
            maxi = max(copy_mer,length)
            mini = min(copy_mer,length)
            if copy_mer != length:
                if maxi <= max(1+mini,1.1*mini):
                    maxi = mini
            if maxi != mini:
                long_mer = match[str(maxi)]
                short_mer = match[str(mini)]
                
                for i in range(long_mer.shape[0]):
                    for j in range(short_mer.shape[0]):
                        short_data_mer = data[str(mini)][short_mer[j]]['patterns'][0]
                        long_data_mer = data[str(maxi)][long_mer[i]]['patterns'][0]
                        
                        short_pattern = short_data_mer['pattern_searched']
                        long_pattern = long_data_mer['pattern_searched']
                        
                        short_consensus = short_data_mer['consensus']
                        long_consensus = long_data_mer['consensus']
                        
                        short_length = mini
                        long_length = maxi
                        
                        
                        mul_factor = math.ceil(long_length/short_length)
                        
                        short_pattern2 = short_pattern + short_pattern
                        short_consensus2 = short_consensus + short_consensus
                
                        if mul_factor == 2:
                            if (long_pattern in short_pattern2 or long_consensus in short_consensus2):
                                maxi = mini
                                
                        if mul_factor > 2:
                            if (short_pattern2 in long_pattern or short_consensus2 in long_consensus):
                                maxi = mini
                                
            data_mer = data[str(mini)]
            data_mer[match[str(mini)][0]]['modified'] = 'true'
            modified_data[str(mini)] = data_mer
            
            data_mer = data[str(maxi)]
            data_mer[match[str(maxi)][0]]['modified'] = 'true'
            modified_data[str(maxi)] = data_mer
            
            for m in match.keys():
                if match[m].shape[0] > 0:
                    if not (int(m) == maxi or int(m) == mini):
                        data_mer = data[m]
                        data_mer[match[m][0]]['modified'] = 'false'
                        modified_data[m] = data_mer
                else:
                    data_mer = data[m]
                    data_mer[0]['modified'] = 'true'
                    modified_data[m] = data_mer
                    
    modified_file['data'] = modified_data
    modified_file['mers'] = file['mers']
    modified_file['region_info'] = file['region_info']
    #
    if end <= 300:
        modified_file['end'] = "true"
    else:
        modified_file['end'] = "false"
    
    
    
    with open(modified_json,'w') as f:
        json.dump(modified_file,f)    
            


        
        
            
            
        
    

