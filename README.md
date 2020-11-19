# SRFv2
SRFv2 is a program to find repeats using various digital signal processing (DSP) techniques on a given DNA sequence


## Requirements
- Operating system: Linux/Mac
- Python
- Numpy
- Scipy
- Matplotlib
- weblogo-3.7.1
- stockwell 0.1.0

### For Linux:
- stockwell 0.1.0 (pip install stockwell)

### For MacOS:
- Download stockwell-0.1.0.tar.gz. Replace linux2 with darwin in 'setup.py' and then install it
- python3.8 ./setup.py build
- sudo python3.8 ./setup.py install

### General Instructions
- File should be in FASTA format (file name should have .fasta extension)
- Output directory should not have any underscore symbol

### Run the following command on command line

#### For fast algorithm
SRFv2_fast.py fasta_sequence_file min_mer max_mer percent_match output_directory [options]

#### For exhaustive algorithm
SRFv2_exhaustive.py fasta_sequence_file min_mer max_mer percent_match output_directory [options]

#### Options are
- -c: Copy Number
- -p: Peak value threshold

## Graphical Abstract
![Graphical Abstract](https://drive.google.com/uc?export=view&id=1ZEf3ojPCj27ng97LUnrFvjfL38j83jhX)
