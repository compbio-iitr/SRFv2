import re
import json
from reading_json import modification_json
from downloadable_file_exhaustive import downloadable

def _get_nucleotide_contributions(line):
    patterns = line.split(' ')
    contribution = [{'A': 0, 'C': 0, 'T': 0, 'G': 0} for i in range(len(patterns[0]))]

    for pattern in patterns:
        # print pattern
        for idx, nt in enumerate(pattern):
            contribution[idx][nt] += 1
    return contribution

def pattern2html(seq_file, locations,region_no, mer, no, ps, cons, noc, score):
    fatsa_file = open(seq_file, 'r')
    fatsa_file.readline()
    genome = ''
    for line in fatsa_file:
        genome += line.strip()
    fatsa_file.close()
    genome = genome.upper()
    genome = re.sub(r'[^ACGT]', '', genome)
    if len(locations) >= 2:
        p1, p2 = locations[0]-1, locations[-1]+mer-1
        fmt = '%%0%dd'%(len(str(locations[-1])))
        html_string = '<div class="pattern" id="%s%d">' % (ps.lower(), region_no)
        html_string += '<div class="pattern-title primary">'
        html_string += ps
        html_string += '</div><div class="pattern-info">'
        html_string += '<span class="pattern-detail">Pattern Searched: <span class="match dna-sequence">%s</span></span><br />' % ps
        html_string += '<span class="pattern-detail">Consensus Pattern: <span class="dna-sequence">%s</span></span><br />' % cons
        html_string += '<span class="pattern-detail">Number of Copies: %d</span><br />' % noc
        html_string += '<span class="pattern-detail">Score: %s</span><br />' % score
        html_string += '</div><div class="dna-sequence">'
        idx = p1
        loc_idx = 0
        count = 0
        while idx < p2:
            if loc_idx > 0 and idx == locations[loc_idx-1] - 1 + mer:
                html_string += '</span>'
            if loc_idx < len(locations) and idx == locations[loc_idx] - 1:
                html_string += '<span class="match">'
                loc_idx += 1
            if count%10 == 0 and count > 0:
                html_string += "&nbsp;"
            if count%60 == 0:
                if not count == 0:
                    html_string += '<span class="droid-sans">' + fmt%(count+locations[0]-1) + '</span>' + "&nbsp;<br />" + '<span class="droid-sans">' +  fmt%(count+1+locations[0]-1) + '</span>' + "&nbsp;"
                else:
                    html_string += '<span class="droid-sans">' + fmt%(count+locations[0]) + '</span>' + "&nbsp;"
            html_string += genome[idx]
            idx += 1
            count += 1
        pattern_file = seq_file.replace('.fasta', '/') + 'pattern%d%d%d.html'%(region_no, mer, no)
        html_string += '<br /><br /></div></div>'

        with open(pattern_file, 'w') as pf:
            pf.write(html_string)
        return pattern_file[pattern_file.rfind('/'):]
    return None

def split_and_append(text, split_by, append_text):
    html_string = ''
    length = len(text)
    div = length/split_by
    for k in xrange(1, div+1):
        html_string = html_string + text[(k-1)*split_by:k*split_by] + append_text
    html_string = html_string + text[div*split_by:length]
    return html_string

#def pat2json(filename, jsonfilename):
#    json_string = {'data' : {}, 'mers': []}
#    mers = []
#    with open(filename, 'r') as pat:
#        idx = 0
#        lines = pat.readlines()
#        while idx < (len(lines)):
#            line = lines[idx].strip()
#            if line.startswith('#'):
#                result = re.search(r'^#:REGION:(\d+):(\d+):(\d+)$', line)
#                frm, upto, mer = result.group(1), result.group(2), result.group(3)
#                
#                
#                region_dict = {}
#                region_dict['from'] = int(frm)
#                region_dict['upto'] = int(upto)
#                region_dict['patterns'] = []
#                
#                idx += 1
#                if idx < len(lines):
#
#                    line = lines[idx].strip()
#                hello = 0
#                while idx < len(lines) and not line.startswith('#'):
#                    hello = 1
#                    pattern_dict = {}
#
#                    result = re.search(r'^<:([ACTG]+):(\d+):([ACTGactg/]+)$', line)
#
#                    pattern_dict['pattern_searched'] = result.group(1)
#                    pattern_dict['consensus'] = result.group(3)
#                    pattern_dict['number_of_copies'] = int(result.group(2))
#
#                    idx += 1
#
#                    line = lines[idx].strip()
#
#                    pattern_dict['contribution'] = _get_nucleotide_contributions(line[2:])
#
#                    idx += 1
#
#                    scores = [float(score) for score in lines[idx][2:].split(' ')]
#                    score = sum(scores)/len(scores)
#
#                    pattern_dict['score'] = '%.4f' % score
#
#                    region_dict['patterns'].append(pattern_dict)
#
#                    idx += 2
#
#                    if idx < len(lines):
#                        line = lines[idx].strip()
#                idx -= 1
#                if hello == 1:
#                  mers.append(mer)
#                  mer_data = json_string['data'].get(mer)
#                  if mer_data is None:
#                    json_string['data'][mer] = []
#                  json_string['data'][mer].append(region_dict)
#                  
#
#            idx += 1
#        mers = sorted(list(set(mers)))
#        json_string['mers'] = mers
#        with open(jsonfilename, 'w') as jsonfile:
#            jsonfile.write(json.dumps(json_string))#, sort_keys=True, indent=4, separators=(',', ': ')))

def pat2json(filename, jsonfilename, modified_name, download_filename, region_no, region_begin, region_end,start,end):
    seq_file = filename[:filename.rfind('/')]
    seq_file = seq_file +'.fasta'
    json_string = {'data' : {}, 'mers': [], 'region_info': {'region_no':region_no, 'region_begin':start, 'region_end':end}}
    mers = []
    pattern_no = 1
    with open(filename, 'r') as pat:
        idx = 0
        lines = pat.readlines()
        while idx < (len(lines)):
            line = lines[idx].strip()
            if line.startswith('#'):
                if idx == len(lines) - 1:
                    idx += 1
                    continue

                line_next = lines[idx+1].strip()
                if line_next.startswith('#'):
                    idx += 1
                    continue

                result = re.search(r'^#:REGION:(\d+):(\d+):(\d+)$', line)
                frm, upto, mer = result.group(1), result.group(2), result.group(3)
                
                region_dict = {}
                region_dict['from'] = int(frm) + region_begin - 1
                region_dict['upto'] = int(upto) + region_begin - 1
                region_dict['patterns'] = []
                

                idx += 1

                if idx < len(lines):
                    line = lines[idx].strip()
                hello = 0
                while idx < len(lines) and not line.startswith('#'):
                    hello = 1
                    pattern_dict = {}

                    result = re.search(r'^<:([ACTG]+):(\d+):([ACTGWSRYKMBDHVN]+)$', line)

                    pattern_dict['pattern_searched'] = result.group(1)
                    pattern_dict['consensus'] = result.group(3)
                    pattern_dict['number_of_copies'] = int(result.group(2))

                    idx += 1

                    line = lines[idx].strip()

                    pattern_dict['contribution'] = _get_nucleotide_contributions(line[2:])

                    idx += 1

                    scores = [float(score) for score in lines[idx][2:].split(' ')]
                    score = sum(scores)/len(scores)

                    pattern_dict['score'] = '%.4f' % score

                    idx += 1

                    line = lines[idx].strip()[2:]
                    locations = [int(loc) + region_begin - 1 for loc in line.split(' ')]
                    pattern_file = pattern2html(seq_file, locations, region_no, int(mer), pattern_no, pattern_dict['pattern_searched'], pattern_dict['consensus'], pattern_dict['number_of_copies'], pattern_dict['score'])
                    if pattern_file:
                        pattern_dict['pattern_file'] = pattern_file
                    else:
                        pattern_dict['pattern_file'] = ''
                    pattern_no += 1

                    region_dict['patterns'].append(pattern_dict)

                    idx += 1

                    if idx < len(lines):
                        line = lines[idx].strip()
                idx -= 1
                if hello == 1:
                    mers.append(mer)
                    mer_data = json_string['data'].get(mer)
                    if mer_data is None:
                        json_string['data'][mer] = []
                    json_string['data'][mer].append(region_dict)
            idx += 1
        mers = sorted(list(set(mers)))
        json_string['mers'] = mers
        with open(jsonfilename, 'w') as jsonfile:
            jsonfile.write(json.dumps(json_string))#, sort_keys=True, indent=4, separators=(',', ': ')))
        modification_json(jsonfilename,modified_name,region_end)
        downloadable(modified_name,filename,download_filename)
        




#if __name__ == '__main__':
#    pat2json('pat1.file', 'data.json')
#filename = './greater/pat2.file'
#jsonfilename = './pretty/data2.json'
#modified_name = './pretty/modified1.json'
#region_no = 1
#region_begin = 1
#region_end = 10000
#start = 9500
#end = 19500
#pat2json(filename, jsonfilename, modified_name, region_no, region_begin, region_end,start,end)