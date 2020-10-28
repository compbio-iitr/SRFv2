#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 02:00:21 2020

@author: backpropagator
"""

import os
from glob import glob
import json
import numpy as np

def html(output_dir):    
    html_string = ''
    if os.path.exists(output_dir+'/modified1.json'):
        json_files = glob(output_dir+'/modified*.json')
        no_of_regions = len(json_files)
        
        if no_of_regions > 1:
            html_string = '<div class="mer-row">[Fast]<span style="float: right;">[<a href="#" class="redundant_button">Hide/Unhide Redundant</a>] </span></div>'
            html_string = html_string + '<div class="regions">'
            html_string = html_string + '<table class="table table-bordered"><tbody>'
            html_string = html_string + '<tr class="text-center success heading"><td colspan="2" class="success">Regions</td></tr>'
            
            for region_no in range(1,no_of_regions+1):
                json_string = output_dir + '/data' + str(region_no) + '.json'
                with open(json_string) as f:
                    data = json.load(f)
                region_begin = data['region_info']['region_begin']
                region_end = data['region_info']['region_end']
                r_range = str(region_begin) + ':' + str(region_end)
                html_string = html_string + '<tr class="text-center droid-sans"><td>Region '
                html_string = html_string + str(region_no)
                html_string = html_string + '</td><td>'
                html_string = html_string + '<a class="in-div-anchor" data-number="'
                html_string = html_string + str(region_no) 
                html_string = html_string + '" href="#region' 
                html_string = html_string + str(region_no) 
                html_string = html_string + '">' 
                html_string = html_string + str(r_range)
                html_string = html_string + '</a></td></tr>'
            
            html_string = html_string + '</tbody></table></div>'
        
        for region_no in range(1,no_of_regions+1):
            if no_of_regions > 1:
                html_string = html_string + '<div id="show-details' + str(region_no) + '" class="show-details mer-row primary">Region ' + str(region_no) + '</div>'
            else:
                html_string = ''
            
            json_string = output_dir + '/modified' + str(region_no) + '.json'
            with open(json_string) as f:
                data = json.load(f)
            
            pattern_files = []
            mers = data['mers']
            length_condition = data['end']
            
            if no_of_regions > 1:
                html_string = html_string + '<div class="mer-row" id="region'
                html_string = html_string + str(region_no)
                html_string = html_string + '">[View Full Sequence Spectrum: <a href="#" class="view-winspectra" data-mer="0" data-spectype="fft" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">FFT</a> | <a href="#" class="view-winspectra" data-mer="0" data-spectype="dct" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">DCT</a> | <a href="#" class="view-winspectra" data-mer="0" data-spectype="dst" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">DST</a> | <a href="#" class="view-winspectra" data-mer="0" data-spectype="strans" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">AST</a>]</div>';
                
            else:
                html_string = html_string + '<div class="mer-row" id="region'
                html_string = html_string + str(region_no)
                html_string = html_string + '">[Fast][View Full Sequence Spectrum: <a href="#" class="view-winspectra" data-mer="0" data-spectype="fft" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">FFT</a> | <a href="#" class="view-winspectra" data-mer="0" data-spectype="dct" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">DCT</a> | <a href="#" class="view-winspectra" data-mer="0" data-spectype="dst" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">DST</a> | <a href="#" class="view-winspectra" data-mer="0" data-spectype="strans" data-region="'
                html_string = html_string + str(region_no)
                html_string = html_string + '">AST</a>]<span style="float: right;">[<a href="#" class="redundant_button">Hide/Unhide Redundant</a>]</span>'
            
            html_string = html_string + '<table class="table table-bordered"><tbody>'
            html_string = html_string + '<tr class="text-center success heading">'
            html_string = html_string + '<td class="success">Pattern Searched</td><td class="success">Region</td><td class="success">Consensus</td><td class="success">&#35; Copies</td><td class="success">Score</td></tr>'
    
            mers = np.array(mers)
            mers = mers.astype(np.uint16)
            mers.sort()
            
            for mer in mers:
                mer_data = data['data'][str(mer)]
                condition_value = False
                
                for region in mer_data:
                    condition = region['modified']
                    if condition == 'false':
                        condition_value = condition_value or False
                    else:
                        condition_value = condition_value or True
                        
                if not condition_value:
                    if length_condition =='true':
                        html_string = html_string + '<tr class="info primary_redundant hide_unhide"><td class="mer-row text-center" colspan="5">' + str(mer) + '-mer Repeats'
                    else:
                        html_string = html_string + '<tr class="info primary_redundant hide_unhide"><td class="mer-row text-center" colspan="5">' + str(mer) + '-mer Repeats (View Window Spectrum: '
                else:
                    if length_condition=='true':
                        html_string = html_string + '<tr class="info"><td class="mer-row text-center primary" colspan="5">' + str(mer) + '-mer Repeats'
                    else:
                        html_string = html_string + '<tr class="info"><td class="mer-row text-center primary" colspan="5">' + str(mer) + '-mer Repeats (View Window Spectrum: '
                
                wspec_files = glob(output_dir + '/*wspec' + str(region_no) + '_' + str(mer) + '.png')
                no_of_files = len(wspec_files)
                fft = 0
                dct = 0
                
                if os.path.exists(output_dir + '/fourier_wspec' +  str(region_no) + '_' + str(mer) + '.png'):
                    html_string = html_string + ' <a href="#" class="view-winspectra" data-spectype="fft" data-mer="' + str(mer) + '" data-region="' + str(region_no) +'">FFT</a> '
                    fft = 1
                    if no_of_files > 1:
                        html_string = html_string + ' | '
                
                if os.path.exists(output_dir + '/dct_wspec' + str(region_no) + '_' + str(mer) + '.png'):
                    html_string = html_string + ' <a href="#" class="view-winspectra" data-spectype="dct" data-mer="' + str(mer) + '" data-region="'+ str(region_no) +'">DCT</a> '
                    if fft == 1 and no_of_files > 2:
                        html_string = html_string + ' | '
                    if fft == 0 and no_of_files > 1:
                        html_string = html_string + ' | '
                    dct = 1
                
                if os.path.exists(output_dir + '/dst_wspec' + str(region_no) + '_' + str(mer) + '.png'):
                    html_string = html_string + ' <a href="#" class="view-winspectra" data-spectype="dst" data-mer="' + str(mer) + '" data-region="'+ str(region_no) +'">DST</a> '
                    if fft == 1 and dct == 1 and no_of_files > 3:
                        html_string = html_string + ' | '
                    if(((fft == 0 and dct == 1) or (fft == 1 and dct == 0)) and no_of_files > 2):
                        html_string = html_string + ' | '
                    if(fft == 0 and dct == 0 and no_of_files > 1):
                        html_string = html_string + ' | '
                
                if os.path.exists(output_dir + '/stransform_wspec' + str(region_no) + '_' + str(mer) + '.png'):
                    html_string = html_string + ' <a href="#" class="view-winspectra" data-spectype="strans" data-mer="' + str(mer) + '" data-region="'+ str(region_no) +'">AST</a> '
                    
                if length_condition=='true':
                    html_string = html_string + '</td></tr>'
                else:
                    html_string = html_string + ')</td></tr>'
                
                mer_data = data['data'][str(mer)]
                
                for region in mer_data:
                    patterns = region['patterns']
                    condition = region['modified']
                    
                    for pattern in patterns:
                        if condition == 'false':
                            html_string = html_string + '<tr class="text-center droid-sans redundant hide_unhide">'
                        else:
                            html_string = html_string + '<tr class="text-center droid-sans">'
                        
                        html_string = html_string + '<td class="dna-sequence">' + pattern['pattern_searched'] + '</td>'
                        tag = pattern['pattern_searched'].lower()
                        tag = tag + str(region_no)
                        html_string = html_string + '<td><a class="in-div-anchor" data-number="' + str(region_no) + '" href="#' + str(tag) + '">' + str(region['from']) + ':' + str(region['upto']) + '</a></td>'
                        html_string = html_string + '<td class="dna-sequence">' + pattern['consensus']
                        pattern_index = mer_data.index(region) + 1
                        html_string = html_string + '&nbsp; &nbsp;<button class="btn logo_button" data-spectype="logo" data-number="' + str(pattern_index) + '" data-region="' + str(region_no) + '"data-mer="' + str(mer) + '"><span class="red_logo">L</span><span class="orange_logo">o</span><span class="green_logo">g</span><span class="blue_logo">o</span></button>'
                        html_string = html_string + '</td>'
                        html_string = html_string + '<td>' + str(pattern['number_of_copies']) + '</td>'
                        html_string = html_string + '<td>' + str(pattern['score']) + '</td>'
                        html_string = html_string + '</tr>'
                        
                        if len(pattern['pattern_file']) > 0:
                            pattern_string = ''
                            if condition == 'false':
                                pattern_string = pattern_string +  '<div class= "redundant_detailed hide_unhide"> '
                            else:
                                pattern_string = pattern_string + '<div> '
                            
                            with open(output_dir+pattern['pattern_file']) as FPO:
                                for line in FPO:
                                    pattern_string = pattern_string + line
                            
                            pattern_string = pattern_string + '</div>'
                            pattern_files.append(pattern_string)
                    
            html_string = html_string + '</tbody></table>'
            html_string = html_string + '<div class="redundant-text hide_unhide">* Rows highlighted in yellow are redundant repeats.</div>'
            html_string = html_string + '<div id="show-details' + str(region_no) + '" class="show-details mer-row primary" data-region="' + str(region_no) + '" data-action="show"><a href="#"><i class="fa fa-angle-double-down"></i>&nbsp;Detailed Results</a></div>'
            html_string = html_string + '<div id="details' + str(region_no) + '" style="display:none;">'
            
            for pattern_file in pattern_files:
                html_string = html_string + pattern_file
            
            html_string = html_string + '</div>'
        
        pattern_repeat_files = glob(output_dir + '/pattern*.html')
        if len(pattern_repeat_files) == 0:
            html_string = '<div class="alert alert-info">'
            html_string = html_string + '<h4 class="alert-heading">No Repeats Found!</h4><p>No repeating patterns were found in the submitted sequence.</p></div>'
        
    final_string = ''
    
    with open('./assets/layout1.html') as FPO:
        for line in FPO:
            final_string = final_string + line
    
    final_string = final_string + html_string
    
    with open('./assets/layout2.html') as FPO:
        for line in FPO:
            final_string = final_string + line
    
    with open(output_dir + '/result_file.html','w') as FPO:
        FPO.write(final_string)
        
        
    	    
                        
                        
                        
                
                    
                    
            
            
        