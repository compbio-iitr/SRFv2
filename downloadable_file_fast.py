#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 21:49:22 2020

@author: backpropagator
"""

import json

#json_file = './mod_75.json'
#pat_file = './pat_75.file'
#downloadable_file = './download.file'

def downloadable(json_file,pat_file,downloadable_file):
    with open(json_file) as f:
        file = json.load(f)
    
    data_file = file['data']
    chars = []
    with open(pat_file) as f:
        pattern = ''
        for line in f:
            chars.append(line.strip())
     
    with open(downloadable_file,'a') as f:
        f.write('OUTPUT OF SRFv2 - FAST\n\n')
        f.write('Format of the Results\n\n')
        f.write('LENGTH:    Repeat Length\n')
        f.write('DETAILS:   Pattern Searched | Region | Consensus | #Copies | Score\n')
        f.write('POSITIONS: Location of repeats detected\n')
        f.write('\n')
        for i in range(len(chars)):
            if chars[i][0] == '#':
                splits = chars[i].split(':')
                mer = splits[-1]
                beg = splits[2]
                end = splits[3]
                pattern_pat = chars[i+1].split(':')[1]
                ncopies = chars[i+1].split(':')[2]
                regions = chars[i+4].split(':')[1]
                length = len(data_file[mer])
                for j in range(length):
                    pattern = data_file[mer][j]['patterns'][0]['pattern_searched']
                    if pattern == pattern_pat:
                        final_pattern = data_file[mer][j]['patterns'][0]['pattern_searched']
                        final_consensus = data_file[mer][j]['patterns'][0]['consensus']
                        final_score = data_file[mer][j]['patterns'][0]['score']
                        redundant = data_file[mer][j]['modified']
                
                if redundant == 'true':    
                    f.write('LENGTH:    '+mer+'\n')
                else:
                    f.write('LENGTH:    '+mer+' (Redundant) \n')
                f.write('DETAILS:   '+final_pattern+' | '+beg+':'+end+' | ' + final_consensus+' | '+ncopies+' | '+final_score+'\n')
                f.write('POSITIONS: '+regions+'\n')
                f.write('\n')
        
        
        
        
        
        
       
            
            