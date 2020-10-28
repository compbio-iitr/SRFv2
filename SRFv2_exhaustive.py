#!/home/deepak/anaconda3/envs/compbio/bin/python

import os
import sys
import shutil
import re
import math
import numpy as np
import time
from glob import glob
from html_file_exhaustive import html

#import parallel_winspectra
import prettifyoutput_exhaustive
import matplotlib.pyplot as plt
from multiprocessing import Pool as ThreadPool
import multiprocessing as mp
from scipy.fftpack import dct, dst
from stockwell.smt import st

output_dir = ''
window = 100
copy_number = 5
peak_value = 4
PI = math.pi
scale = 2.85535
measure_speed = True
DFT = 1
FFT = 2
algo = FFT
parallel = False


def get_inputs():
    cnt = 0
    inputs = sys.argv
    cnt = len(inputs)
    if cnt < 6:
        print('USAGE: SRFv2_exhaustive.py <fasta sequence file> <min mer> <max mer> <%%match> <output directory> [options]\n\n')
        # TODO add outputs
        kill()
    return inputs

def kill():
    exit('Exiting')

def error_check():
    global output_dir
    output_dir = sys.argv[5]
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    global peak_value
    global copy_number

    for idx, val in enumerate(sys.argv):
        if idx == 1:
            fatsa_file = open(val, 'r')
            line = fatsa_file.readline()
            fatsa_file.close()
            if '>' in line and line.index('>') != 0:
                print('ERROR: Sequence file not in FASTA format\n')
                kill()
        elif idx == 2:
            if re.search(r'[^0-9]', val) or int(val) < 2 or int(val) > int(sys.argv[3]):
                print('ERROR: Range of minimum mer >= 2 and <= maximum mer & an Integer\n\n')
                kill()
        elif idx == 3:
            if re.search(r'[^0-9]', val) or int(val) < int(sys.argv[2]):
                print('ERROR: Range of maximum mer >= minimum mer & an Integer\n\n')
                kill()
        elif idx == 4:
            if int(val) <= 1 or int(val) > 100:
                print('ERROR: Percent Match should be between 1 - 100\n\n')
                kill()
        elif idx == 6:
            cd = val[0:2]
            if cd == '-c':
                copy_number = int(val[3:])
                if int(copy_number) < 2:
                    print('ERROR: Copy number should be integer >= 2\n\n')
                    kill()
            elif cd == '-p':
                peak_value = int(val[3:])
                if int(peak_value) <= 0:
                    print('ERROR: Spectral peak should be a value > 0\n\n')
                    kill()
        elif idx == 7:
            cd = val[0:2]
            if cd == '-c':
                copy_number = int(val[3:])
                if int(copy_number) < 2:
                    print('ERROR: Copy number should be integer >= 2\n\n')
                    kill()
            elif cd == '-p':
                peak_value = int(val[3:])
                if int(peak_value) <= 0:
                    print('ERROR: Spectral peak should be a value > 0\n\n')
                    kill()


def win_fourier_spectra(argument):
    genome = argument[0]
    final = argument[1]
    outfile =argument[2]
    mer = argument[3]
    i14 = argument[4]
    length_of_subseq = len(genome)
    if length_of_subseq <= 600:
        window = 100
    else:
        if mer >=2 and mer <= 150:
            window = 300
        else:
            window = 600
    global PI
    global scale
    # mul_factor = 2.85535
    FWO = open(outfile, 'w')
    f = 1/float(mer)
    rmer = 0
    s = 0

    REGIONS = {}
    for ii in range(length_of_subseq):
        
        iii = ii+1 + i14*19000
        rty = genome[ii:ii+window]
        win_final = final[ii:ii+window]
        if len(rty) < window:
            break
        
        sqA = rty.count('A')/float(window)
        sqC = rty.count('C')/float(window)
        sqG = rty.count('G')/float(window)
        sqT = rty.count('T')/float(window)

        sqA = sqA*sqA
        sqC = sqC*sqC
        sqG = sqG*sqG
        sqT = sqT*sqT

        sq = sqA+sqC+sqG+sqT
        gsq = 1+(1/float(window))
        hsq = gsq - sq
        av_spectra = hsq/float(window)
        
        xA = yA = xC = yC = xG = yG = xT = yT = 0

        for i in range(window):
            count = i+1
            b1 = PI*f
            sine = math.sin(b1*count*2)
            cosine = math.cos(b1*count*2)

            sA = cA = sC = cC = sG = cG = sT = cT = 0

            if rty[i] == 'A':
                sA = 0.1260*sine
                cA = 0.1260*cosine
            else:
                sA = cA = 0
            if rty[i] == 'C':
                sC = 0.1340*sine
                cC = 0.1340*cosine
            else:
                sC = cC = 0
            if rty[i] == 'G':
                sG = 0.0806*sine
                cG = 0.0806*cosine
            else:
                sG = cG = 0
            if rty[i] == 'T':
                sT = 0.1335*sine
                cT = 0.1335*cosine
            else:
                sT = cT = 0

            xA += sA
            yA += cA
            xC += sC
            yC += cC
            xG += sG
            yG += cG
            xT += sT
            yT += cT

        xA**=2
        yA**=2
        xC**=2
        yC**=2
        xG**=2
        yG**=2
        xT**=2
        yT**=2

        pA = (xA+yA)/(window**2)
        pC = (xC+yC)/(window**2)
        pG = (xG+yG)/(window**2)
        pT = (xT+yT)/(window**2)

        P = pA + pC + pG + pT

        Pk = P*scale*100

        Peak = Pk/av_spectra
        FWO.write('%d %f\n'%(iii, Peak))

        if Peak >= peak_value and s == 0:
            if REGIONS.get(rmer) == None:
                REGIONS[rmer] = {}
            REGIONS[rmer][0] = iii
            s = 1
        if Peak < peak_value and s == 1:
            REGIONS[rmer][1] = iii - 1
            rmer += 1
            s = 0
    if s == 1:
        REGIONS[rmer][1] = length_of_subseq + i14*19000
        rmer += 1

    FWO.write('# Mers %d were found in the regions\n' % (mer))

    REG = {}
    gomer = 0
    for i in range(rmer):
        if not (REGIONS[i][1] - REGIONS[i][0]) < mer*2:
            REG[gomer] = {}
            REG[gomer][0] = REGIONS[i][0]
            REG[gomer][1] = REGIONS[i][1]
            gomer += 1

    GK = {}
    half_window = window

    for i in range(gomer):
        GK[i] = {}
        GK[i][0] = max(1+i14*19000, REG[i][0])
        GK[i][1] = min(length_of_subseq+i14*19000, REG[i][1]+half_window)

    bner = 0
    start, stop = {}, {}
    for i in range(gomer):
        if i > 0:
            if GK[i][0] <= stop[bner-1]:
                stop[bner-1] = GK[i][1]
            else:
                start[bner] = GK[i][0]
                stop[bner] = GK[i][1]
                bner += 1
        else:
            start[bner] = GK[i][0]
            stop[bner] = GK[i][1]
            bner += 1

    for i in range(bner):
        if mer <= 50:
            if (stop[i]+1) - start[i] < copy_number*mer:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))
        else:
            if (stop[i]+1)-start[i] < 100:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))

    FWO.close()

def win_dct_spectra(argument):
    genome = argument[0]
    final = argument[1]
    outfile =argument[2]
    mer = argument[3]
    i14 = argument[4]
    
    
    length_of_subseq = len(genome)
    if length_of_subseq <= 600:
        window = 100
    else:
        if mer >=2 and mer <= 150:
            window = 300
        else:
            window = 600
    global PI
    global scale
    FWO = open(outfile, 'w')
    f = 1/float(mer)
    rmer = 0
    s = 0

    REGIONS = {}
    peak_values = []
    for ii in range(length_of_subseq):
        
        iii = ii+1 + i14*19000
        rty = genome[ii:ii+window]
        win_final = final[ii:ii+window]
        if len(rty) < window:
            break
        
        power_spectrum = dct(win_final)
        freq = np.fft.fftfreq(power_spectrum.shape[0])
    
        power_spectrum = ((1.0/window**2)*abs(power_spectrum)**2)[0:1+int(len(freq)/2)]

        av_spectra = np.mean(power_spectrum[1:])

        yA = yC = yG = yT = 0

        for i in range(window):
            count = (2*i+1)/2
            b1 = PI*f
            cosine = math.cos(b1*count)

            cA = cC = cG = cT = 0

            if rty[i] == 'A':
                cA = 0.1260*cosine
            else:
                cA = 0
            if rty[i] == 'C':
                cC = 0.1340*cosine
            else:
                cC = 0
            if rty[i] == 'G':
                cG = 0.0806*cosine
            else:
                cG = 0
            if rty[i] == 'T':
                cT = 0.1335*cosine
            else:
                cT = 0

            yA += cA
            yC += cC
            yG += cG
            yT += cT

        yA**=2
        yC**=2
        yG**=2
        yT**=2

        pA = yA/(window**2)
        pC = yC/(window**2)
        pG = yG/(window**2)
        pT = yT/(window**2)

        P = pA + pC + pG + pT
        

        Pk = P*scale

        Peak = Pk/av_spectra
        peak_values.append(Peak)

    y_coords = np.array(peak_values)

    for i in range(y_coords.shape[0]):
        if i <= 6:
            y_coords[i] = (y_coords[i] + y_coords[i+1] + y_coords[i+2] + y_coords[i+3] + y_coords[i+4])/5
        elif i >= y_coords.shape[0]-7:
            y_coords[i] = (y_coords[i-4] + y_coords[i-3] + y_coords[i-2] + y_coords[i-1] + y_coords[i])/5
        else:
            y_coords[i] = (y_coords[i-3] + y_coords[i-2] + y_coords[i-1] + y_coords[i] + y_coords[i+1] + y_coords[i+2] + y_coords[i+3])/7


    for ii in range(length_of_subseq):
        iii = ii + 1 + i14*19000
        rty = genome[ii:ii+window]
        win_final = final[ii:ii+window]
        if len(rty) < window:
            break
        FWO.write('%d %f\n'%(iii, y_coords[ii]))
                  
        if y_coords[ii] >= peak_value and s == 0:
            if REGIONS.get(rmer) == None:
                REGIONS[rmer] = {}
            REGIONS[rmer][0] = iii
            s = 1
        if y_coords[ii] < peak_value and s == 1:
            REGIONS[rmer][1] = iii - 1
            rmer += 1
            s = 0
    if s == 1:
        REGIONS[rmer][1] = length_of_subseq + i14*19000
        rmer += 1

    FWO.write('# Mers %d were found in the regions\n' % (mer))

    REG = {}
    gomer = 0
    for i in range(rmer):
        if not (REGIONS[i][1] - REGIONS[i][0]) < mer*2:
            REG[gomer] = {}
            REG[gomer][0] = REGIONS[i][0]
            REG[gomer][1] = REGIONS[i][1]
            gomer += 1

    GK = {}
    half_window = window

    for i in range(gomer):
        GK[i] = {}
        GK[i][0] = max(1 + i14*19000, REG[i][0])
        GK[i][1] = min(length_of_subseq + i14*19000, REG[i][1]+half_window)

    bner = 0
    start, stop = {}, {}
    for i in range(gomer):
        if i > 0:
            if GK[i][0] <= stop[bner-1]:
                stop[bner-1] = GK[i][1]
            else:
                start[bner] = GK[i][0]
                stop[bner] = GK[i][1]
                bner += 1
        else:
            start[bner] = GK[i][0]
            stop[bner] = GK[i][1]
            bner += 1

    for i in range(bner):
        if mer <= 50:
            if (stop[i]+1) - start[i] < copy_number*mer:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))
        else:
            if (stop[i]+1)-start[i] < 100:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))

    FWO.close()

def win_dst_spectra(argument):
    genome = argument[0]
    final = argument[1]
    outfile =argument[2]
    mer = argument[3]
    i14 = argument[4]
    length_of_subseq = len(genome)
    if length_of_subseq <= 600:
        window = 100
    else:
        if mer >=2 and mer <= 150:
            window = 300
        else:
            window = 600
    global PI
    global scale
    FWO = open(outfile, 'w')
    f = 1/float(mer)
    rmer = 0
    s = 0

    REGIONS = {}
    peak_values = []
    for ii in range(length_of_subseq):
        
        iii = ii+1 + i14*19000
        rty = genome[ii:ii+window]
        win_final = final[ii:ii+window]
        if len(rty) < window:
            break
        
        power_spectrum = dst(win_final)
        freq = np.fft.fftfreq(power_spectrum.shape[0])

        power_spectrum = ((1.0/window**2)*abs(power_spectrum)**2)[0:1+int(len(freq)/2)]

        av_spectra = np.mean(power_spectrum[1:])

        yA = yC = yG = yT = 0

        for i in range(window):
            count = (2*i+1)/2
            b1 = PI*f
            sine = math.sin(b1*count)

            cA = cC = cG = cT = 0

            if rty[i] == 'A':
                cA = 0.1260*sine
            else:
                cA = 0
            if rty[i] == 'C':
                cC = 0.1340*sine
            else:
                cC = 0
            if rty[i] == 'G':
                cG = 0.0806*sine
            else:
                cG = 0
            if rty[i] == 'T':
                cT = 0.1335*sine
            else:
                cT = 0

            yA += cA
            yC += cC
            yG += cG
            yT += cT

        yA**=2
        yC**=2
        yG**=2
        yT**=2

        pA = yA/(window**2)
        pC = yC/(window**2)
        pG = yG/(window**2)
        pT = yT/(window**2)

        P = pA + pC + pG + pT
        
        

        Pk = P*scale

        Peak = Pk/av_spectra
        
        peak_values.append(Peak)

    y_coords = np.array(peak_values)

    for i in range(y_coords.shape[0]):
        if i <= 8:
            y_coords[i] = (y_coords[i] + y_coords[i+1] + y_coords[i+2] + y_coords[i+3] + y_coords[i+4])/5
        elif i >= y_coords.shape[0]-9:
            y_coords[i] = (y_coords[i-4] + y_coords[i-3] + y_coords[i-2] + y_coords[i-1] + y_coords[i])/5
        else:
            y_coords[i] = (y_coords[i-4] + y_coords[i-3] + y_coords[i-2] + y_coords[i-1] + y_coords[i] + y_coords[i+1] + y_coords[i+2] + y_coords[i+3] + y_coords[i+4])/9


    for ii in range(length_of_subseq):
        iii = ii + 1 + i14*19000
        rty = genome[ii:ii+window]
        win_final = final[ii:ii+window]
        if len(rty) < window:
            break
        FWO.write('%d %f\n'%(iii, y_coords[ii]))
        if y_coords[ii] >= peak_value and s == 0:
            if REGIONS.get(rmer) == None:
                REGIONS[rmer] = {}
            REGIONS[rmer][0] = iii
            s = 1
        if y_coords[ii] < peak_value and s == 1:
            REGIONS[rmer][1] = iii - 1
            rmer += 1
            s = 0
    if s == 1:
        REGIONS[rmer][1] = length_of_subseq + i14*19000
        rmer += 1

    FWO.write('# Mers %d were found in the regions\n' % (mer))

    REG = {}
    gomer = 0
    for i in range(rmer):
        if not (REGIONS[i][1] - REGIONS[i][0]) < mer*2:
            REG[gomer] = {}
            REG[gomer][0] = REGIONS[i][0]
            REG[gomer][1] = REGIONS[i][1]
            gomer += 1

    GK = {}
    half_window = window

    for i in range(gomer):
        GK[i] = {}
        GK[i][0] = max(1 + i14*19000, REG[i][0])
        GK[i][1] = min(length_of_subseq + i14*19000, REG[i][1]+half_window)

    bner = 0
    start, stop = {}, {}
    for i in range(gomer):
        if i > 0:
            if GK[i][0] <= stop[bner-1]:
                stop[bner-1] = GK[i][1]
            else:
                start[bner] = GK[i][0]
                stop[bner] = GK[i][1]
                bner += 1
        else:
            start[bner] = GK[i][0]
            stop[bner] = GK[i][1]
            bner += 1

    for i in range(bner):
        if mer <= 50:
            if (stop[i]+1) - start[i] < copy_number*mer:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))
        else:
            if (stop[i]+1)-start[i] < 100:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))

    FWO.close()
    
def win_stransform_spectra(argument):
    genome = argument[0]
    final = argument[1]
    outfile =argument[2]
    mer = argument[3]
    i14 = argument[4]
    if mer == 2:
        return
    
    length_of_subseq = len(genome)
    if length_of_subseq <= 600:
        window = 100
    else:
        if mer >=2 and mer <= 150:
            window = 300
        else:
            window = 600
    global PI
    global scale
    FWO = open(outfile, 'w')
    f = 1/float(mer)
    
    index = f*100
    floor = math.floor(index)
    ceil = math.ceil(index)
    rmer = 0
    m = 0
    REGIONS = {}
    Peak = 0

    for ii in range(length_of_subseq):
        tickoff = 0
        peak_values = []
        plot_values = []
        iii = ii+1 + i14*19000
        rty = genome[ii:ii+window]
        win_final = final[ii:ii+window]
        if len(rty) < window:
            break
        
        
        stock = st(win_final)
        stock = abs(stock)
    
        stock = stock[1:]
        binary = np.zeros(stock.shape)
        median = np.mean(stock,axis=-1)

        for i in range(win_final.shape[0]):
            for j in range(stock.shape[0]):
                if stock[j,i]>median[j]:
                    binary[j,i] = 1
            
            
        s = np.sum(binary,axis=-1)
        
        freq = np.fft.fftfreq(window)
        freq = np.around(freq,4)
        median_value = np.median(s)
        ind_beg = np.argwhere(freq == floor/100).reshape(-1)[0]
        ind_end = np.argwhere(freq == ceil/100).reshape(-1)[0]
        for k in range(ind_beg,ind_end+1):
            plot_values.append(s[k])
            if s[k] > median_value:
                tickoff = 1
                peak_values.append(s[k])
        
        if len(peak_values) > 0:
            
            peak_values = np.array(peak_values)
                
            Peak = np.max(peak_values)
            p = np.min(peak_values)
        else:
            try:
                p
            except NameError:
                Peak = 0
            else:
                Peak = p
        
        peak_plot = np.array(plot_values)
        peak_plot = np.max(peak_plot)
        Peak = (4*peak_plot)/median_value
        FWO.write('%d %f\n'%(iii, Peak))

        if tickoff==1 and m == 0:
            if REGIONS.get(rmer) == None:
                REGIONS[rmer] = {}
            REGIONS[rmer][0] = iii
            m = 1
        if tickoff==0 and m == 1:
            REGIONS[rmer][1] = iii - 1
            rmer += 1
            m = 0
    if m == 1:
        REGIONS[rmer][1] = length_of_subseq + i14*19000
        rmer += 1

    FWO.write('# Mers %d were found in the regions\n' % (mer))

    REG = {}
    gomer = 0
    for i in range(rmer):
        if not (REGIONS[i][1] - REGIONS[i][0]) < mer*2:
            REG[gomer] = {}
            REG[gomer][0] = REGIONS[i][0]
            REG[gomer][1] = REGIONS[i][1]
            gomer += 1

    GK = {}
    half_window = window

    for i in range(gomer):
        GK[i] = {}
        GK[i][0] = max(1 + i14*19000, REG[i][0])
        GK[i][1] = min(length_of_subseq + i14*19000, REG[i][1]+half_window)

    bner = 0
    start, stop = {}, {}
    for i in range(gomer):
        if i > 0:
            if GK[i][0] <= stop[bner-1]:
                stop[bner-1] = GK[i][1]
            else:
                start[bner] = GK[i][0]
                stop[bner] = GK[i][1]
                bner += 1
        else:
            start[bner] = GK[i][0]
            stop[bner] = GK[i][1]
            bner += 1

    for i in range(bner):
        if mer <= 50:
            if (stop[i]+1) - start[i] < copy_number*mer:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))
        else:
            if (stop[i]+1)-start[i] < 100:
                continue
            FWO.write('# %d %d\n'%(start[i], stop[i]))

    FWO.close()
def pattern_finder(argument):
    genome = argument[0]
    outfile = argument[1]
    mer = argument[2]
    rstart = argument[3]
    substart = argument[4]
    
    

    length_of_subseq = len(genome)
    thres = float(sys.argv[4])/100.
    patt = ['xxx']


    for i6 in range(length_of_subseq):
        string = genome[i6:i6+mer]
        if len(string) < mer:
            break
        s = 0
        i7 = 0
        while i7 < len(patt):
            if string == patt[i7]:
                s = 1
                break
            i7 += 1
        if s == 0:
            patt.append(string)
    patt = patt[1:]

    max_ = 0
    mxval = 0

    for pat in patt:
        Ppat = ['xx']
        Ppos = [0]
        Pval = [0]
        for i8 in range(length_of_subseq):
            string = genome[i8:i8+mer]
            if len(string) < mer: break
            value = 0
            for i9 in range(mer):
                if pat[i9] == string[i9]:
                    value += 1
            pos = rstart-1 + substart-1 + i8 + 1;
            value = float(value)/mer
            if value == 1:
                if (pos - Ppos[-1]) < mer and Ppos[-1] != 0:
                    if value > Pval[-1]:
                        Pval[-1] = value
                        Ppos[-1] = pos
                        Ppat[-1] = string
                if (pos - Ppos[-1]) < mer and Ppos[-1] == 0:
                    Pval.append(value)
                    Ppos.append(pos)
                    Ppat.append(string)
                if (pos - Ppos[-1]) >= mer:
                    Pval.append(value)
                    Ppos.append(pos)
                    Ppat.append(string)
                i8 = i8 + mer - 1
            elif value >= thres:
                if (pos - Ppos[-1]) < mer and Ppos[-1] != 0:
                    if value > Pval[-1]:
                        Pval[-1] = value
                        Ppos[-1] = pos
                        Ppat[-1] = string
                if (pos - Ppos[-1]) >= mer or Ppos[-1] == 0:
                    Pval.append(value)
                    Ppos.append(pos)
                    Ppat.append(string)
        Ppos = Ppos[1:]
        Pval = Pval[1:]
        Ppat = Ppat[1:]

        meval = 0

        for gal in Pval:
            meval += gal

        copy_num = len(Ppos)
        meval = meval/float(copy_num)

        if copy_num >= copy_number:
            if copy_num > max_:
                consensus = cons(Ppat, mer) 
                with open(outfile, 'w') as FPO:
                    FPO.write('<:%s:%d:%s\n'%(pat, copy_num, consensus))
                    FPO.write('>:'+' '.join(Ppat)+'\n')
                    FPO.write('>:'+' '.join([str(_) for _ in Pval])+'\n')
                    FPO.write('>:'+' '.join([str(_) for _ in Ppos])+'\n')
            elif copy_num == max_:
                if meval == mxval:
                    consensus = cons(Ppat, mer)
                    with open(outfile, 'a') as FPO:
                        FPO.write('<:%s:%d:%s\n'%(pat, copy_num, consensus))
                        FPO.write('>:'+' '.join(Ppat)+'\n')
                        FPO.write('>:'+' '.join([str(_) for _ in Pval])+'\n')
                        FPO.write('>:'+' '.join([str(_) for _ in Ppos])+'\n')
                elif meval > mxval:
                    consensus = cons(Ppat, mer)
                    with open(outfile, 'w') as FPO:
                        FPO.write('<:%s:%d:%s\n'%(pat, copy_num, consensus))
                        FPO.write('>:'+' '.join(Ppat)+'\n')
                        FPO.write('>:'+' '.join([str(_) for _ in Pval])+'\n')
                        FPO.write('>:'+' '.join([str(_) for _ in Ppos])+'\n')
        if max_ < copy_num:
            max_ = copy_num
            mxval = meval
        if max_ == copy_num and mxval < meval:
            max_ = copy_num
            mxval = meval
      
    

def cons(patterns, mer):
    consensus = ''

    for ik1 in range(mer):
        a = c = t = g = 0
        for ik2 in range(len(patterns)):
            fd = patterns[ik2][ik1:ik1+1]
            if fd.upper() == 'A':
                a += 1
            elif fd.upper() == 'C': 
                c += 1
            elif fd.upper() == 'G':
                g += 1
            elif fd.upper() == 'T':
                t += 1
        if a>t and a>g and a>c:
            consensus += 'A'
        if c>t and c>g and c>a:
            consensus += 'C'
        if g>t and g>a and g>c:
            consensus += 'G'
        if t>a and t>g and t>c:
            consensus += 'T'
        if a==c and a>t and a>g:
            consensus += 'M'
        if t==c and t>a and t>g:
            consensus += 'Y'
        if g==c and g>a and g>t:
            consensus += 'S'
        if g==t and g>a and g>c:
            consensus += 'K'
        if g==a and g>c and g>t:
            consensus += 'R'
        if t==a and t>c and t>g:
            consensus += 'W'
        if t==a and t==c and t>g:
            consensus += 'H'
        if t==a and t==g and t>c:
            consensus += 'D'
        if c==a and c==g and c>t:
            consensus += 'V'
        if c==t and c==g and c>a:
            consensus += 'B'
        if c==t and c==g and c==a:
            consensus += 'N'
    return consensus


def fast_fourier_spectra(genome,final, outfile):
    x_coord = []
    y_coord = []

    global peak_value
    minmer, maxmer = int(sys.argv[2]), int(sys.argv[3])
    length_of_subseq = final.shape[0]

    global PI
    global scale
    len_by_3 = length_of_subseq/3

    

    power = np.fft.fft(final)
    
    freq = np.fft.fftfreq(power.shape[-1])

    power_spectrum = ((1.0/length_of_subseq**2)*abs(power)**2)[0:1+int(len(freq)/2)]
    
    av_spectra = np.mean(power_spectrum[1:])
    
    

    FO = open(outfile, 'w')

    MERS = {}
    PEAKS = {}
    vmer = 0

    for idx in range(1, 1+int(len(freq)/2)):
        f = abs(freq[idx])
        Pk = scale*power_spectrum[idx]
        FO.write('%f %f\n'%(f, Pk))
        x_coord.append(f)
        y_coord.append(Pk)
        Peak = Pk/av_spectra
        if Peak > peak_value:
            cntp = 1/f
            cng = int(cntp)

            if cng <= len_by_3:
                if cntp-cng >= 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng+1:
                            s = 1
                            break
                    if s==0 and cng+1 >= minmer and cng+1 <= maxmer:
                        MERS[vmer] = cng+1
                        PEAKS[vmer] = Peak
                        vmer += 1
                elif cntp-cng < 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng:
                            s = 1
                            break
                    if s == 0 and cng >= minmer and cng <= maxmer: 
                        flag = 0
                        if cng > 20 and cng < 60 and MERS.get(vmer-1) is not None and cng < (MERS[vmer-1]+4):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        elif cng >= 60 and MERS.get(vmer-1) is not None and cng<(MERS[vmer-1]+10):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        if flag == 0:
                            MERS[vmer] = cng
                            PEAKS[vmer] = Peak
                            vmer+=1
                   
    pera = '%6.3f'%(genome.count('A')*100/float(length_of_subseq))
    perc = '%6.3f'%(genome.count('C')*100/float(length_of_subseq))
    pert = '%6.3f'%(genome.count('T')*100/float(length_of_subseq))
    perg = '%6.3f'%(genome.count('G')*100/float(length_of_subseq))

    FO.write('# Average Spectra :%f:%s:%s:%s:%s\n'%(av_spectra, pera, perg, perc, pert))
    FO.write('# Mers Found in the DNA sequence\n')
    bmer7 = 0

    for i5 in range(vmer-1,-1,-1):
        if MERS[i5] == bmer7:
            continue
        bmer7 = MERS[i5]
        FO.write('# %d %f\n'%(MERS[i5], PEAKS[i5]))
    FO.close()
    x_coord = np.array(x_coord)
    y_coord = np.array(y_coord)
    if copy_number >= 5:
        y_coord[0] = 0
        y_coord[1] = 0
        y_coord[2] = 0
        y_coord[3] = 0
    else:
        y_coord[0] = 0
        y_coord[1] = 0
    np.save(outfile[:-5]+'_x.npy',x_coord)
    np.save(outfile[:-5]+'_y.npy',y_coord)

   
def dct_spectra(genome,final, outfile):
    x_coord = []
    y_coord = []

    global peak_value
    minmer, maxmer = int(sys.argv[2]), int(sys.argv[3])
    length_of_subseq = final.shape[0]

    global PI
    global scale
    len_by_3 = length_of_subseq/3 

    
    power = dct(final)
    
    freq = np.fft.fftfreq(power.shape[-1])

    power_spectrum = ((1.0/length_of_subseq**2)*abs(power)**2)[0:1+int(len(freq)/2)]
    
    av_spectra = np.mean(power_spectrum[1:])

    FO = open(outfile, 'w')

    MERS = {}
    PEAKS = {}
    vmer = 0

    for idx in range(1, 1+int(len(freq)/2)):
        f = abs(freq[idx])
        Pk = scale*power_spectrum[idx] 
        FO.write('%f %f\n'%(f, Pk))
        x_coord.append(f)
        y_coord.append(Pk)
        Peak = Pk/av_spectra
        if Peak > peak_value:
            cntp = 1/f
            cng = int(cntp)

            if cng <= len_by_3:
                if cntp-cng >= 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng+1:
                            s = 1
                            break
                    if s==0 and cng+1 >= minmer and cng+1 <= maxmer:
                        MERS[vmer] = cng+1
                        PEAKS[vmer] = Peak
                        vmer += 1
                elif cntp-cng < 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng:
                            s = 1
                            break
                    if s == 0 and cng >= minmer and cng <= maxmer: 
                        flag = 0
                        if cng > 20 and cng < 60 and MERS.get(vmer-1) is not None and cng < (MERS[vmer-1]+4):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        elif cng >= 60 and MERS.get(vmer-1) is not None and cng<(MERS[vmer-1]+10):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        if flag == 0:
                            MERS[vmer] = cng
                            PEAKS[vmer] = Peak
                            vmer+=1
    pera = '%6.3f'%(genome.count('A')*100/float(length_of_subseq))
    perc = '%6.3f'%(genome.count('C')*100/float(length_of_subseq))
    pert = '%6.3f'%(genome.count('T')*100/float(length_of_subseq))
    perg = '%6.3f'%(genome.count('G')*100/float(length_of_subseq))

    FO.write('# Average Spectra :%f:%s:%s:%s:%s\n'%(av_spectra, pera, perg, perc, pert))
    FO.write('# Mers Found in the DNA sequence\n')
    bmer7 = 0

    for i5 in range(vmer-1,-1,-1):
        if MERS[i5] == bmer7:
            continue
        bmer7 = MERS[i5]
        FO.write('# %d %f\n'%(MERS[i5], PEAKS[i5]))
    FO.close()
    x_coord = np.array(x_coord)
    y_coord = np.array(y_coord)
    if copy_number >= 5:
        y_coord[0] = 0
        y_coord[1] = 0
        y_coord[2] = 0
        y_coord[3] = 0
    else:
        y_coord[0] = 0
        y_coord[1] = 0
    np.save(outfile[:-5]+'_x.npy',x_coord)
    np.save(outfile[:-5]+'_y.npy',y_coord)

        
def dst_spectra(genome,final, outfile):
    x_coord = []
    y_coord = []

    global peak_value
    global copy_number
    minmer, maxmer = int(sys.argv[2]), int(sys.argv[3])
    length_of_subseq = final.shape[0]

    global PI
    global scale
    len_by_3 = length_of_subseq/3 
    power = dst(final)
    
    freq = np.fft.fftfreq(power.shape[-1])

    power_spectrum = ((1.0/length_of_subseq**2)*abs(power)**2)[0:1+int(len(freq)/2)]
    
    av_spectra = np.mean(power_spectrum[1:])

    FO = open(outfile, 'w')

    MERS = {}
    PEAKS = {}
    vmer = 0

    for idx in range(1, 1+int(len(freq)/2)):
        f = abs(freq[idx])
        Pk = scale*power_spectrum[idx]
        FO.write('%f %f\n'%(f, Pk))
        x_coord.append(f)
        y_coord.append(Pk)
        Peak = Pk/av_spectra
        if Peak > peak_value:
            cntp = 1/f
            cng = int(cntp)

            if cng <= len_by_3:
                if cntp-cng >= 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng+1:
                            s = 1
                            break
                    if s==0 and cng+1 >= minmer and cng+1 <= maxmer:
                        MERS[vmer] = cng+1
                        PEAKS[vmer] = Peak
                        vmer += 1
                elif cntp-cng < 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng:
                            s = 1
                            break
                    if s == 0 and cng >= minmer and cng <= maxmer: 
                        flag = 0
                        if cng > 20 and cng < 60 and MERS.get(vmer-1) is not None and cng < (MERS[vmer-1]+4):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        elif cng >= 60 and MERS.get(vmer-1) is not None and cng<(MERS[vmer-1]+10):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        if flag == 0:
                            MERS[vmer] = cng
                            PEAKS[vmer] = Peak
                            vmer+=1
    pera = '%6.3f'%(genome.count('A')*100/float(length_of_subseq))
    perc = '%6.3f'%(genome.count('C')*100/float(length_of_subseq))
    pert = '%6.3f'%(genome.count('T')*100/float(length_of_subseq))
    perg = '%6.3f'%(genome.count('G')*100/float(length_of_subseq))

    FO.write('# Average Spectra :%f:%s:%s:%s:%s\n'%(av_spectra, pera, perg, perc, pert))
    FO.write('# Mers Found in the DNA sequence\n')
    bmer7 = 0

    for i5 in range(vmer-1,-1,-1):
        if MERS[i5] == bmer7:
            continue
        bmer7 = MERS[i5]
        FO.write('# %d %f\n'%(MERS[i5], PEAKS[i5]))
    FO.close()
    x_coord = np.array(x_coord)
    y_coord = np.array(y_coord)
    if copy_number >= 5 and length_of_subseq >= 10:
        y_coord[0] = 0
        y_coord[1] = 0
        y_coord[2] = 0
        y_coord[3] = 0
    elif length_of_sequence >= 6:
        y_coord[0] = 0
        y_coord[1] = 0
    np.save(outfile[:-5]+'_x.npy',x_coord)
    np.save(outfile[:-5]+'_y.npy',y_coord)

  
def stransform_spectra(genome,final, outfile):
    x_coord = []
    y_coord = []

    global peak_value
    minmer, maxmer = int(sys.argv[2]), int(sys.argv[3])
    length_of_subseq = final.shape[0]


    global PI
    global scale
    len_by_3 = length_of_subseq/3 

    stock = st(final)
    stock = abs(stock)
    
    binary = np.zeros(stock.shape)
    median = np.mean(stock,axis=-1)

    for i in range(final.shape[0]):
        for j in range(stock.shape[0]):
            if stock[j,i]>median[j]:
                binary[j,i] = 1
            
            
    s = np.sum(binary,axis=-1)
    
    freq = np.fft.fftfreq(length_of_subseq)
    
    power_spectrum = s
    
    av_spectra = np.median(power_spectrum[1:])
    
    FO = open(outfile, 'w')

    MERS = {}
    PEAKS = {}
    vmer = 0

    for idx in range(1, 1+int(len(freq)/2)):
        f = abs(freq[idx])
        Pk = power_spectrum[idx] 
        FO.write('%f %f\n'%(f, Pk))
        x_coord.append(f)
        y_coord.append(Pk)

        Peak = Pk/av_spectra
        if Peak >= 1:
            cntp = 1/f
            cng = int(cntp)

            if cng <= len_by_3:
                if cntp-cng >= 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng+1:
                            s = 1
                            break
                    if s==0 and cng+1 >= minmer and cng+1 <= maxmer:
                        MERS[vmer] = cng+1
                        PEAKS[vmer] = Peak
                        vmer += 1
                elif cntp-cng < 0.5:
                    s = 0
                    for i4 in range(vmer):
                        if MERS[i4] == cng:
                            s = 1
                            break
                    if s == 0 and cng >= minmer and cng <= maxmer: 
                        flag = 0
                        if cng > 20 and cng < 60 and MERS.get(vmer-1) is not None and cng < (MERS[vmer-1]+4):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        elif cng >= 60 and MERS.get(vmer-1) is not None and cng<(MERS[vmer-1]+10):
                            if Peak > PEAKS[vmer-1]:
                                MERS[vmer-1] = cng
                                PEAKS[vmer-1] = Peak
                            flag = 1
                        if flag == 0:
                            MERS[vmer] = cng
                            PEAKS[vmer] = Peak
                            vmer+=1
    pera = '%6.3f'%(genome.count('A')*100/float(length_of_subseq))
    perc = '%6.3f'%(genome.count('C')*100/float(length_of_subseq))
    pert = '%6.3f'%(genome.count('T')*100/float(length_of_subseq))
    perg = '%6.3f'%(genome.count('G')*100/float(length_of_subseq))

    FO.write('# Average Spectra :%f:%s:%s:%s:%s\n'%(peak_value, pera, perg, perc, pert))
    FO.write('# Mers Found in the DNA sequence\n')
    bmer7 = 0
    for i5 in range(vmer-1,-1,-1):
        if MERS[i5] == bmer7:
            continue
        bmer7 = MERS[i5]
        FO.write('# %d %f\n'%(MERS[i5], PEAKS[i5]))
    FO.close()
    x_coord = np.array(x_coord)
    y_coord = np.array(y_coord)
    np.save(outfile[:-5]+'_x.npy',x_coord)
    np.save(outfile[:-5]+'_y.npy',y_coord)

  
def generate_full_spectra_plot(transform, filenum, length_of_subseq, length):
    spectra_file = (output_dir + '/'+transform+'_spectra%d.file'%filenum)
    png_name = spectra_file[:spectra_file.rfind('.')]
    x_coords = np.load(output_dir + '/'+transform+'_spectra%d_x.npy'%filenum)
    y_coords = np.load(output_dir + '/'+transform+'_spectra%d_y.npy'%filenum)
    max_y = np.max(y_coords)
    
    with open(spectra_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('#'):
                x, y = map(float, line.split(' '))
              
            elif not (line.startswith('# Mers') or line.startswith('# Ave')):
                line = line.replace('#', '').strip()
                mer = int(line.split(' ')[0])
                wspec_file = spectra_file[:spectra_file.rfind('/')] + '/'+transform+'_wspec%d_%d.file'%(filenum, mer) 
                generate_win_spectra_plot(wspec_file, mer, length_of_subseq,transform)
    plt.figure()
    count = 0
    m = max_y
    while m < 1:
        m = m * 10
        count = count + 1

    y_coords = y_coords*math.pow(10,count)
    
    plt.plot(x_coords, y_coords)
    if transform == 'fourier':
        plt.title("Full Sequence FFT Spectrum")
    elif transform == 'stransform':
        plt.title("Full Sequence AST Spectrum")
    else:
        plt.title("Full Sequence " + transform.upper() + " Spectrum")
    
    
    plt.xlabel("Frequency")
    if transform != "stransform":
        plt.ylabel("Power Spectrum (x10^-" + str(count)+")")
    else:
        plt.ylabel("Power Spectrum")
        
    plt.xlim((0, 0.5))
    png_name = '%s.png'%png_name
    
    plt.savefig(png_name)
    plt.close()

def generate_win_spectra_plot(win_spectra_file, mer, length,transform):
    if mer == 2 and transform == 'stransform':
        return
    png_name = win_spectra_file[:win_spectra_file.rfind('.')]
    x_coords = []
    y_coords = []
    max_y = -1
    count = 0
    with open(win_spectra_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('#'):
                x, y = map(float, line.split(' '))
                x_coords.append(x)
                y_coords.append(y)
                max_y = max(max_y, y)
            if line.startswith('#'):
                count = count + 1
    plt.figure()
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)

   

    if max_y < 4:
        return
    if count == 1:
        return
    
        
    plt.plot(x_coords, y_coords)
    if transform == 'fourier':
        plt.title("Window FFT Spectrum for %d-mer" % mer)
    elif transform == 'stransform':
        plt.title("Window AST Spectrum for %d-mer" % mer)
    else:
        plt.title("Window " + transform.upper() + " Spectrum for %d-mer" % mer)
    plt.xlabel("Nucleotide Position")
    plt.ylabel("Signal-to-noise ratio")
    plt.xlim((np.min(x_coords), np.max(x_coords)))
    png_name = '%s.png'%png_name
    plt.savefig(png_name)
    plt.close()

def generate_full_spectra_plot_300(transform, filenum, length_of_subseq):
    spectra_file = (output_dir + '/'+transform+'_spectra%d.file'%filenum)
    png_name = spectra_file[:spectra_file.rfind('.')]
    x_coords = np.load(output_dir + '/'+transform+'_spectra%d_x.npy'%filenum)
    y_coords = np.load(output_dir + '/'+transform+'_spectra%d_y.npy'%filenum)
    max_y = np.max(y_coords)
    
    with open(spectra_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('#'):
                x, y = map(float, line.split(' '))
            elif not (line.startswith('# Mers') or line.startswith('# Ave')):
                line = line.replace('#', '').strip()
                
    plt.figure()
    count = 0
    m = max_y
    while m < 1 and m!=0:
        m = m * 10
        count = count + 1

    y_coords = y_coords*math.pow(10,count)
    
    plt.plot(x_coords, y_coords)
    if transform == 'fourier':
        plt.title("Full Sequence FFT Spectrum")
    elif transform == 'stransform':
        plt.title("Full Sequence AST Spectrum")
    else:
        plt.title("Full Sequence " + transform.upper() + " Spectrum")
    
    
    plt.xlabel("Frequency")
    if transform != "stransform":
        plt.ylabel("Power Spectrum (x10^-" + str(count)+")")
    else:
        plt.ylabel("Power Spectrum")
        
    plt.xlim((0, 0.5))
    png_name = '%s.png'%png_name
    print(png_name)
    plt.savefig(png_name)
    plt.close()

progpath = os.getcwd()

def parallel_pattern(args):
    global window
    global peak_value
    global copy_number

    pool = ThreadPool(mp.cpu_count())
    
    pool.map(pattern_finder, args)
    pool.close()
    pool.join()
    
def parallel_fourier(args):
    global window
    global peak_value
    global copy_number

    pool = ThreadPool(mp.cpu_count())
    
    pool.map(win_fourier_spectra, args)
    pool.close()
    pool.join()

def parallel_dct(args):
    global window
    global peak_value
    global copy_number

    pool = ThreadPool(mp.cpu_count())
    
    pool.map(win_dct_spectra, args)
    pool.close()
    pool.join()

def parallel_dst(args):
    global window
    global peak_value
    global copy_number

    pool = ThreadPool(mp.cpu_count())
    
    pool.map(win_dst_spectra, args)
    pool.close()
    pool.join()

def parallel_stransform(args):
    global window
    global peak_value
    global copy_number

    pool = ThreadPool(mp.cpu_count())
    
    pool.map(win_stransform_spectra, args)
    pool.close()
    pool.join()


if __name__ == '__main__':
    t1 = time.time()
    inputs = get_inputs()
    print('[i] Arguments: %s'%(', '.join(inputs)))
    error_check()
    print('[*] No Input Errors Detected\n[*] Starting Program')
    genome = ''
    fatsa_file = open(sys.argv[1], 'r')
    fatsa_file.readline()
    for line in fatsa_file:
        genome += line.strip()
    fatsa_file.close()
    genome = genome.upper()
    genome = re.sub(r'[^ACGT]', '', genome)
    
    length_of_sequence = len(genome)
    arr = []
    for i in genome:
        arr.append(i)
    arr = np.array(arr)
    final = np.zeros((length_of_sequence,))
    final[arr=='A'] = 0.1260
    final[arr=='T'] = 0.1335
    final[arr=='C'] = 0.1340
    final[arr=='G'] = 0.0806

    step = 1
    print('[i] Length of Sequence: {}'.format(length_of_sequence))

    if length_of_sequence > 25000:
        filenum = 0
        WED = []
        fir = 0
        for i14 in range(0, length_of_sequence, 19000):
            string = genome[i14:i14+20000]
            balance = len(string)

            if balance <= 5000: #why this?
                string += genome[i14+20000:]

            genome1 = string

            lenstr = len(genome1)
            filenum += 1
            
            

            patfile = output_dir + '/pat%d.file' % filenum

            if os.path.exists(patfile):
                os.remove(patfile)

            initreg = i14 + 1
            finreg = initreg + lenstr - 1
            fir += 1

            spectra_file = output_dir + '/spectra%d.file' % filenum
            
            arr = []
            for i in genome1:
                arr.append(i)
            arr = np.array(arr)
            final = np.zeros(arr.shape)
            final[arr=='A'] = 0.1260
            final[arr=='T'] = 0.1335
            final[arr=='C'] = 0.1340
            final[arr=='G'] = 0.0806
            fast_fourier_spectra(genome1, final, output_dir + '/fourier_spectra%d.file'%filenum)
            dct_spectra(genome1,final,output_dir+'/dct_spectra%d.file'%filenum)
            dst_spectra(genome1,final,output_dir+'/dst_spectra%d.file'%filenum)
            stransform_spectra(genome1,final,output_dir+'/stransform_spectra%d.file'%filenum)
            
            if os.path.exists(output_dir + '/pat%d.file' % filenum):
                os.remove(output_dir + '/pat%d.file' % filenum)
            
            pattern_args = []
            fft_window = []
            with open(output_dir + '/fourier_spectra%d.file'%filenum,'r') as FII:
                for line2 in FII:
                    line2 = line2.strip()
                    col1 = re.split(r'\s+', line2)
                        
                    if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                        mer = col1[1]
                        fft_window.append([genome1, final, output_dir + '/fourier_wspec%d_%s.file'%(filenum, mer), int(mer),filenum-1])
            parallel_fourier(fft_window)
        
            dct_window = []
            with open(output_dir + '/dct_spectra%d.file'%filenum,'r') as FII:
                for line2 in FII:
                    line2 = line2.strip()
                    col1 = re.split(r'\s+', line2)
                        
                    if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                        mer = col1[1]
                        dct_window.append([genome1, final, output_dir + '/dct_wspec%d_%s.file'%(filenum, mer), int(mer),filenum-1])
            
            parallel_dct(dct_window)
            
            dst_window = []
            with open(output_dir + '/dst_spectra%d.file'%filenum,'r') as FII:
                for line2 in FII:
                    line2 = line2.strip()
                    col1 = re.split(r'\s+', line2)
                        
                    if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                        mer = col1[1]
                        dst_window.append([genome1, final, output_dir + '/dst_wspec%d_%s.file'%(filenum, mer), int(mer),filenum-1])
            parallel_dst(dst_window)
            
            stransform_window = []


            with open(output_dir + '/stransform_spectra%d.file'%filenum,'r') as FII:
                for line2 in FII:
                    line2 = line2.strip()
                    col1 = re.split(r'\s+', line2)
                        
                    if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                        mer = col1[1]
                        stransform_window.append([genome1, final, output_dir + '/stransform_wspec%d_%s.file'%(filenum, mer), int(mer),filenum-1])
            parallel_stransform(stransform_window)
            
            fourier_windows = glob(output_dir+'/fourier_wspec'+str(filenum)+'*')
            dct_windows = glob(output_dir+'/dct_wspec'+str(filenum)+'*')
            dst_windows = glob(output_dir+'/dst_wspec'+str(filenum)+'*')
            stransform_windows = glob(output_dir+'/stransform_wspec'+str(filenum)+'*')
            fourier_windows.sort()
            dct_windows.sort()
            dst_windows.sort()
            stransform_windows.sort()
            mers_length = []
            print(fourier_windows)
            for f in fourier_windows:
                mers_length.append(int(f.split('_')[2].split('.')[0]))
            for f in dct_windows:
                mers_length.append(int(f.split('_')[2].split('.')[0]))
            for f in dst_windows:
                mers_length.append(int(f.split('_')[2].split('.')[0]))
            for f in stransform_windows:
                mers_length.append(int(f.split('_')[2].split('.')[0]))
            mers_length = np.array(mers_length)
            mers_length = np.unique(mers_length)
            mers_length.sort()
            final_regions=[]
            for mer in mers_length:
                regions = []
                if os.path.exists(output_dir + '/fourier_wspec%d_%s.file'%(filenum, mer)):
                    with open(output_dir + '/fourier_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                        for lsg in FG:
                            lsg = lsg.strip()
                            codw = re.split(r'\s+', lsg)
                            if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                                regions.append([int(codw[1]),int(codw[2])])
                if os.path.exists(output_dir + '/dct_wspec%d_%s.file'%(filenum, mer)):
                    with open(output_dir + '/dct_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                        for lsg in FG:
                            lsg = lsg.strip()
                            codw = re.split(r'\s+', lsg)
                            if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                                regions.append([int(codw[1]),int(codw[2])])
                if os.path.exists(output_dir + '/dst_wspec%d_%s.file'%(filenum, mer)):
                    with open(output_dir + '/dst_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                        for lsg in FG:
                            lsg = lsg.strip()
                            codw = re.split(r'\s+', lsg)
                            if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                                regions.append([int(codw[1]),int(codw[2])])
                if os.path.exists(output_dir + '/stransform_wspec%d_%s.file'%(filenum, mer)):
                    with open(output_dir + '/stransform_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                        for lsg in FG:
                            lsg = lsg.strip()
                            codw = re.split(r'\s+', lsg)
                            if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                                regions.append([int(codw[1]),int(codw[2])])
                                
                if len(regions) > 0:
                
                    regions = np.array(regions)
                    regions = np.sort(regions,axis=0)
                    merged_regions = []

                    merge = np.copy(regions)
                    for i in range(1,merge.shape[0]):
                        for j in range(0,i):
                    
                            if merge[i,0]>=merge[j,0] and merge[i,0]<=merge[j,1]+100:
                                merge[i,0] = merge[j,0] 
                                if merge[i,1] > merge[j,1]:
                                    merge[j,1] = merge[i,1]
            
                    uniq_initial = np.unique(merge[:,0])
                    uniq_final = np.unique(merge[:,1])
            
                    for i in range(uniq_initial.shape[0]):
                        merged_regions.append([uniq_initial[i],uniq_final[i]])
                
                    merged_regions = np.array(merged_regions)
                    for i in range(merged_regions.shape[0]):
                
                        dif = int(merged_regions[i,1])-(int(merged_regions[i,0])-1)
                        string2 = genome[(int(merged_regions[i,0])-1):(int(merged_regions[i,0])-1)+dif]
                        merfile = output_dir + '/subpat%d-%s.res' % (filenum, mer)
                        parallel_arg = [string2,merfile,int(mer),initreg-i14,int(merged_regions[i,0])]
                        final_regions.append([int(merged_regions[i,0]),int(merged_regions[i,1]),int(mer)])
                        pattern_args.append(parallel_arg)
            final_regions = np.array(final_regions)
        
            parallel_pattern(pattern_args)
            list_mers = glob(output_dir+'/subpat%d-'%(filenum)+'*.res')
            final_mers = []
            for l in list_mers:
                final_mers.append(int(l.split('-')[1][:-4]))
                
            final_mers.sort()
                
    
            with open(output_dir + '/pat%d.file'%filenum,'a') as FPO:
                for j in range(len(final_mers)):
                    with open(output_dir+'/subpat%d-%s.res'%(filenum,final_mers[j]),'r') as sub:
                        pattern = ''
                        for line in sub:
                            line = line.strip()
                            pattern += line
                        
                        pattern_split = pattern.split('<:')
                        for i in range(1,len(pattern_split)):
                            pattern_full = pattern_split[i]
                            splits = pattern_full.split('>:')
                            regions = splits[-1]
                            split_regions = regions.split(' ')
                            beg = int(split_regions[0])
                            end = int(split_regions[-1]) + final_mers[j] - 1
                            copy_number = int(splits[0].split(':')[1])
                        
                            
                            first = '<:' + splits[0]
                            second = '>:' + splits[1]
                            third = '>:' + splits[2]
                            fourth = '>:' + splits[3]
                            FPO.write('#:REGION:%s:%s:%s\n'%(beg, end, final_mers[j]))
                            FPO.write(first+'\n')
                            FPO.write(second+'\n')
                            FPO.write(third+'\n')
                            FPO.write(fourth+'\n')
                            weblogo = splits[1]
                            list_patterns = weblogo.split(' ')
                            with open(output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i),'a') as logo:    
                                for k in range(len(list_patterns)):
                                    logo.write('>'+str(k+1)+'\n')
                                    logo.write(list_patterns[k]+'\n')
                            os.system("weblogo --format PNG --color" + " '#32CD32' A 'Purine' --color" + " '#0000ff' C 'Pyrimidine' --color red T 'Pyrimidine' --color orange G 'Purine' < " + output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i) + ' > ' + output_dir+'/weblogo%d_%s_%s.png'%(filenum,final_mers[j],i))
                    
                            
                            
                    
                        del pattern, pattern_full, pattern_split
                        
            
            prettifyoutput_exhaustive.pat2json(output_dir + '/pat%d.file'%filenum, output_dir + '/data%d.json'%filenum,output_dir + '/modified%d.json'%filenum,output_dir+'/SRFv2_results%d.file'%(filenum),filenum, 1, len(genome1), initreg, initreg + len(genome1)-1)
            generate_full_spectra_plot('fourier', filenum, len(genome1),'long') 
            generate_full_spectra_plot('dct', filenum, len(genome1),'long') 
            generate_full_spectra_plot('dst', filenum, len(genome1),'long') 
            generate_full_spectra_plot('stransform', filenum, len(genome1),'long')
            html(output_dir)

            if balance <= 5000:
                i14 = length_of_sequence
    elif length_of_sequence > 300 and length_of_sequence <= 25000:
        filenum = 1
       
        fast_fourier_spectra(genome, final, output_dir + '/fourier_spectra%d.file'%filenum)
        dct_spectra(genome,final,output_dir+'/dct_spectra%d.file'%filenum)
        dst_spectra(genome,final,output_dir+'/dst_spectra%d.file'%filenum)
        stransform_spectra(genome,final,output_dir+'/stransform_spectra%d.file'%filenum)
        
        if os.path.exists(output_dir + '/pat%d.file' % filenum):
            os.remove(output_dir + '/pat%d.file' % filenum)

        pattern_args = []
        fft_window = []
        with open(output_dir + '/fourier_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    fft_window.append([genome, final, output_dir + '/fourier_wspec%d_%s.file'%(filenum, mer), int(mer),0])
        parallel_fourier(fft_window)
        
        dct_window = []
        with open(output_dir + '/dct_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    dct_window.append([genome, final, output_dir + '/dct_wspec%d_%s.file'%(filenum, mer), int(mer),0])
        
        parallel_dct(dct_window)
        
        dst_window = []
        with open(output_dir + '/dst_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    dst_window.append([genome, final, output_dir + '/dst_wspec%d_%s.file'%(filenum, mer), int(mer),0])
        parallel_dst(dst_window)
        
        stransform_window = []


        with open(output_dir + '/stransform_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    stransform_window.append([genome, final, output_dir + '/stransform_wspec%d_%s.file'%(filenum, mer), int(mer),0])
        parallel_stransform(stransform_window)
        
        fourier_windows = glob(output_dir+'/fourier_wspec*')
        dct_windows = glob(output_dir+'/dct_wspec*')
        dst_windows = glob(output_dir+'/dst_wspec*')
        stransform_windows = glob(output_dir+'/stransform_wspec*')
        fourier_windows.sort()
        dct_windows.sort()
        dst_windows.sort()
        stransform_windows.sort()
        mers_length = []
        for f in fourier_windows:
            mers_length.append(int(f.split('_')[2].split('.')[0]))
        for f in dct_windows:
            mers_length.append(int(f.split('_')[2].split('.')[0]))
        for f in dst_windows:
            mers_length.append(int(f.split('_')[2].split('.')[0]))
        for f in stransform_windows:
            mers_length.append(int(f.split('_')[2].split('.')[0]))
        mers_length = np.array(mers_length)
        mers_length = np.unique(mers_length)
        mers_length.sort()
        final_regions=[]
        for mer in mers_length:
            regions = []
            if os.path.exists(output_dir + '/fourier_wspec%d_%s.file'%(filenum, mer)):
                with open(output_dir + '/fourier_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                    for lsg in FG:
                        lsg = lsg.strip()
                        codw = re.split(r'\s+', lsg)
                        if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                            regions.append([int(codw[1]),int(codw[2])])
            if os.path.exists(output_dir + '/dct_wspec%d_%s.file'%(filenum, mer)):
                with open(output_dir + '/dct_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                    for lsg in FG:
                        lsg = lsg.strip()
                        codw = re.split(r'\s+', lsg)
                        if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                            regions.append([int(codw[1]),int(codw[2])])
            if os.path.exists(output_dir + '/dst_wspec%d_%s.file'%(filenum, mer)):
                with open(output_dir + '/dst_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                    for lsg in FG:
                        lsg = lsg.strip()
                        codw = re.split(r'\s+', lsg)
                        if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                            regions.append([int(codw[1]),int(codw[2])])
            if os.path.exists(output_dir + '/stransform_wspec%d_%s.file'%(filenum, mer)):
                with open(output_dir + '/stransform_wspec%d_%s.file'%(filenum, mer),'r') as FG:
                    for lsg in FG:
                        lsg = lsg.strip()
                        codw = re.split(r'\s+', lsg)
                        if '#' in lsg and lsg.index('#') == 0 and codw[1] != 'Mers':
                            regions.append([int(codw[1]),int(codw[2])])
            
            if len(regions) > 0:
                
                regions = np.array(regions)
                regions = np.sort(regions,axis=0)
                merged_regions = []

                merge = np.copy(regions)
                for i in range(1,merge.shape[0]):
                    for j in range(0,i):
                    
                        if merge[i,0]>=merge[j,0] and merge[i,0]<=merge[j,1]+100:
                            merge[i,0] = merge[j,0] 
                            if merge[i,1] > merge[j,1]:
                                merge[j,1] = merge[i,1]
              
            
                uniq_initial = np.unique(merge[:,0])
                uniq_final = np.unique(merge[:,1])
            
                for i in range(uniq_initial.shape[0]):
                    merged_regions.append([uniq_initial[i],uniq_final[i]])
                
                merged_regions = np.array(merged_regions)
                for i in range(merged_regions.shape[0]):
                
                    dif = int(merged_regions[i,1])-(int(merged_regions[i,0])-1)
                    string2 = genome[(int(merged_regions[i,0])-1):(int(merged_regions[i,0])-1)+dif]
                    merfile = output_dir + '/subpat%d-%s.res' % (filenum, mer)
                    parallel_arg = [string2,merfile,int(mer),1,int(merged_regions[i,0])]
                    final_regions.append([int(merged_regions[i,0]),int(merged_regions[i,1]),int(mer)])
                    pattern_args.append(parallel_arg)
        final_regions = np.array(final_regions)
        
        parallel_pattern(pattern_args)
         
        list_mers = glob(output_dir+'/subpat%d-'%(filenum)+'*.res')
        final_mers = []
        for l in list_mers:
            final_mers.append(int(l.split('-')[1][:-4]))
            
        final_mers.sort()
            

        with open(output_dir + '/pat%d.file'%filenum,'a') as FPO:
            for j in range(len(final_mers)):
                with open(output_dir+'/subpat%d-%s.res'%(filenum,final_mers[j]),'r') as sub:
                    pattern = ''
                    for line in sub:
                        line = line.strip()
                        pattern += line
                    
                    pattern_split = pattern.split('<:')
                    for i in range(1,len(pattern_split)):
                        pattern_full = pattern_split[i]
                        splits = pattern_full.split('>:')
                        regions = splits[-1]
                        split_regions = regions.split(' ')
                        beg = int(split_regions[0])
                        end = int(split_regions[-1]) + final_mers[j] - 1
                        copy_number = int(splits[0].split(':')[1])
                    
                        
                        first = '<:' + splits[0]
                        second = '>:' + splits[1]
                        third = '>:' + splits[2]
                        fourth = '>:' + splits[3]
                        FPO.write('#:REGION:%s:%s:%s\n'%(beg, end, final_mers[j]))
                        FPO.write(first+'\n')
                        FPO.write(second+'\n')
                        FPO.write(third+'\n')
                        FPO.write(fourth+'\n')
                        weblogo = splits[1]
                        list_patterns = weblogo.split(' ')
                        with open(output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i),'a') as logo:    
                            for k in range(len(list_patterns)):
                                logo.write('>'+str(k+1)+'\n')
                                logo.write(list_patterns[k]+'\n')
                        os.system("weblogo --format PNG --color" + " '#32CD32' A 'Purine' --color" + " '#0000ff' C 'Pyrimidine' --color red T 'Pyrimidine' --color orange G 'Purine' < " + output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i) + ' > ' + output_dir+'/weblogo%d_%s_%s.png'%(filenum,final_mers[j],i))
                    
                        
                    
                    del pattern, pattern_full, pattern_split
        prettifyoutput_exhaustive.pat2json(output_dir + '/pat%d.file'%filenum, output_dir + '/data%d.json'%filenum,output_dir + '/modified%d.json'%filenum,output_dir+'/SRFv2_results%d.file'%filenum,filenum, 1, len(genome),1, len(genome))
        generate_full_spectra_plot('fourier', filenum, len(genome),'short') 
        generate_full_spectra_plot('dct', filenum, len(genome),'short') 
        generate_full_spectra_plot('dst', filenum, len(genome),'short') 
        generate_full_spectra_plot('stransform', filenum, len(genome),'short') 
        html(output_dir)
#
    elif length_of_sequence <= 300:
        filenum = 1
    

        fast_fourier_spectra(genome, final, output_dir + '/fourier_spectra%d.file'%filenum)
        dct_spectra(genome,final,output_dir+'/dct_spectra%d.file'%filenum)
        dst_spectra(genome,final,output_dir+'/dst_spectra%d.file'%filenum)
        stransform_spectra(genome,final,output_dir+'/stransform_spectra%d.file'%filenum)
        
        if os.path.exists(output_dir + '/pat%d.file' % filenum):
            os.remove(output_dir + '/pat%d.file' % filenum)

        pattern_args = []
        mers = []
        with open(output_dir + '/fourier_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    mers.append(int(mer))
                    
        with open(output_dir + '/dct_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    mers.append(int(mer))
                    
        with open(output_dir + '/dst_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    mers.append(int(mer))
                    
        with open(output_dir + '/stransform_spectra%d.file'%filenum,'r') as FII:
            for line2 in FII:
                line2 = line2.strip()
                col1 = re.split(r'\s+', line2)
                    
                if '#' in line2 and line2.index('#') == 0 and col1[1] != 'Average' and col1[1] != 'Mers':
                    mer = col1[1]
                    mers.append(int(mer))
        
        mers = np.array(mers)
        mers = np.unique(mers)
        
        for m in mers:
            merfile = output_dir + '/subpat%d-%s.res' % (filenum, m)
            parallel_arg = [genome,merfile,int(m),1,1]
            pattern_args.append(parallel_arg)
            
        parallel_pattern(pattern_args)
        mers.sort()
        list_mers = glob(output_dir+'/subpat%d-'%(filenum)+'*.res')
        final_mers = []
        for l in list_mers:
            final_mers.append(int(l.split('-')[1][:-4]))
            
        final_mers.sort()
            

        with open(output_dir + '/pat%d.file'%filenum,'a') as FPO:
            for j in range(len(final_mers)):
                with open(output_dir+'/subpat%d-%s.res'%(filenum,final_mers[j]),'r') as sub:
                    pattern = ''
                    for line in sub:
                        line = line.strip()
                        pattern += line
                    
                    pattern_split = pattern.split('<:')
                    for i in range(1,len(pattern_split)):
                        pattern_full = pattern_split[i]
                        # print(pattern_full)
                        splits = pattern_full.split('>:')
                        # print(splits)
                        regions = splits[-1]
                        split_regions = regions.split(' ')
                        beg = int(split_regions[0])
                        end = int(split_regions[-1]) + final_mers[j] - 1
                        copy_number = int(splits[0].split(':')[1])
                    
                        
                        first = '<:' + splits[0]
                        second = '>:' + splits[1]
                        third = '>:' + splits[2]
                        fourth = '>:' + splits[3]
                        FPO.write('#:REGION:%s:%s:%s\n'%(beg, end, final_mers[j]))
                        FPO.write(first+'\n')
                        FPO.write(second+'\n')
                        FPO.write(third+'\n')
                        FPO.write(fourth+'\n')
                        weblogo = splits[1]
                        list_patterns = weblogo.split(' ')
                        with open(output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i),'a') as logo:    
                            for k in range(len(list_patterns)):
                                logo.write('>'+str(k+1)+'\n')
                                logo.write(list_patterns[k]+'\n')
                        os.system('/home/deepak/anaconda3/envs/compbio/bin/weblogo --format PNG < ' + output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i) + ' > ' + output_dir+'/weblogo%d_%s_%s.png'%(filenum,final_mers[j],i))
                        os.system('/home/deepak/anaconda3/envs/compbio/bin/weblogo --format PNG < ' + output_dir+'/weblogo%d_%s_%s.fa'%(filenum,final_mers[j],i) + ' > ' + output_dir+'/download/Web_Logos/weblogo%d_%s_%s.png'%(filenum,final_mers[j],i))
                    
                        
        prettifyoutput_exhaustive.pat2json(output_dir + '/pat%d.file'%filenum, output_dir + '/data%d.json'%filenum,output_dir + '/modified%d.json'%filenum,output_dir+'/SRFv2_results%d.file'%filenum,filenum, 1, len(genome),1,len(genome))
        generate_full_spectra_plot_300('fourier', filenum, len(genome)) 
        generate_full_spectra_plot_300('dct', filenum, len(genome)) 
        generate_full_spectra_plot_300('dst', filenum, len(genome)) 
        generate_full_spectra_plot_300('stransform', filenum, len(genome))
        html(output_dir)
                    
    t2 = time.time()
    if measure_speed:
        print('[i] Total Run Time: %.3fs'%(t2-t1))
    with open(output_dir + '/completed', 'w') as compl:
        compl.write('1')
    with open(output_dir + '/exhaustive.txt', 'w') as compl:
        compl.write('1')


