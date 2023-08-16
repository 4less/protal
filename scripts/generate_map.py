import math
import sys
import glob
import os
from os.path import isfile, join
import textdistance
import argparse



def OmitExtension(file: str, extensions=None):
    if extensions == None:
        extensions = {".fq", ".fastq", ".fq.gz", ".fastq.gz"}
        
    for ext in extensions:
        if file.endswith(ext):
            return file.replace(ext, "")
    
    return file
    
    
def get_all_files_recursive(directory):
    all_files = []

    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            all_files.append(file_path)

    return all_files


def IsPair(read1: str, read2: str):
    r1 = os.path.split(read1)
    r2 = os.path.split(read2)
    
    if ' '.join(r1[:-1]) != ' '.join(r2[:-1]):
        return False
    
    if textdistance.hamming(read1, read2) != 1 or len(read1) != len(read2):
        return False

    for i in range(1, len(read1)):
        c1 = read1[i]
        c2 = read2[i]
        if c1 != c2 and not read1[i-1].isdigit() and not read2[i-1].isdigit()\
           and (i == (len(read1) - 1) or not read1[i+1].isdigit() and not read2[i+1].isdigit()):
            if c1 == "1" and c2 == "2" or c1 == "2" and c2 == "1":
                return True
            
    return False

def GetRelevantIndices(path_list):
    path_list_split = [path.split('/') for path in path_list]
    
    min_len = min([len(path) for path in path_list_split])
    max_len = max([len(path) for path in path_list_split])
    
    
    relevant_indices = [*range(min_len, max_len)]
    
    for i in range(0, min_len):
        if not all([path[i] == path_list_split[0][i] for path in path_list_split]):
            relevant_indices.append(i)
    
    return relevant_indices


def FirstDiffIndex(str1, str2):
    if len(str1) != len(str2):
        return -1
    
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            return i
    return -1

def RemovePair(read1, read2):
    diff_idx = FirstDiffIndex(read1, read2)
    sample = read1
    remove_set = {'_', '/', '.'}
    if diff_idx == len(read1) - 1:
        i = diff_idx - 1
        while i > 0 and read1[i] in remove_set:
            i -= 1
        sample = read1[:i+1]
        
    return sample


def GenerateSampleName(read1: str, read2: str, indices):
    r1 = read1.split('/')[:-1]
    r2 = read2.split('/')[:-1]
    
    r1_base = OmitExtension(os.path.basename(read1))
    r2_base = OmitExtension(os.path.basename(read2))
    
    sample_name = '_'.join([r1[i] for i in indices if i < len(r1)])
    sample_name += '_' + RemovePair(r1_base, r2_base)
    
    return sample_name
    
def GenerateSamplesFromReads(read_files):
    read_files = sorted(read_files)
    
    indices = GetRelevantIndices(read_files)
    
    samples = []
    
    for i in range(0, len(read_files) - 1):
        first_path = read_files[i]
        potential_second_path = read_files[i+1]
        
        first = OmitExtension(os.path.basename(first_path))
        potential_second = OmitExtension(os.path.basename(potential_second_path))
        
        if IsPair(first, potential_second):
            sample_name = GenerateSampleName(first_path, potential_second_path, indices)
            
            samples.append({
                "#SAMPLEID": sample_name,
                "FIRST": first_path,
                "SECOND": potential_second_path,
                "SAM": "{}.sam".format(sample_name),
                "PROFILE": "{}.profile".format(sample_name),
                "PREFIX": sample_name,
            })
            i += 1
        # else:
        #     #print(first)
        
    return samples

class ProtalMap:
    MISC_OUTPUT_DIR = "#OUTPUT_DIR"
    STRAIN_OUTPUT_DIR = "#STRAIN_OUTPUT_DIR"
    SAM_OUTPUT_DIR = "#SAM_OUTPUT_DIR"
    PROFILE_OUTPUT_DIR = "#PROFILE_OUTPUT_DIR"
    INPUT_DIR = "#INPUT_DIR"


    DEFAULT_MISC_OUTPUT_DIR = "misc"
    DEFAULT_STRAIN_OUTPUT_DIR = "strains"
    DEFAULT_PROFILE_OUTPUT_DIR = "profiles"
    DEFAULT_SAM_OUTPUT_DIR = "alignments"
    
    SAMPLEID_COL = "#SAMPLEID"
    FIRST_COL = "FIRST"
    SECOND_COL = "SECOND"
    SAM_COL = "SAM"
    PREFIX_COL = "PREFIX"
    PROFILE_COL = "PROFILE"
    
    HEADER_VARS = [SAMPLEID_COL, FIRST_COL, SECOND_COL, SAM_COL, PROFILE_COL, PREFIX_COL]
    
    def __init__(self, input_dir):
        self.misc_output_dir = ""
        self.strain_output_dir = ""
        self.profile_output_dir = ""
        self.sam_output_dir = ""
        self.input_dir = input_dir
        
        self.column_header = [self.SAMPLEID_COL, self.FIRST_COL, self.SECOND_COL, self.SAM_COL, self.PROFILE_COL, self.PREFIX_COL]
        self.samples = []
        
    def SetOutputsToDefault(self, output_base):
        self.strain_output_dir = os.path.join(output_base, self.DEFAULT_STRAIN_OUTPUT_DIR)
        self.profile_output_dir = os.path.join(output_base, self.DEFAULT_PROFILE_OUTPUT_DIR)
        self.sam_output_dir = os.path.join(output_base, self.DEFAULT_SAM_OUTPUT_DIR)
        self.misc_output_dir = os.path.join(output_base, self.DEFAULT_MISC_OUTPUT_DIR)

    def ValidSample(self, sample):
        keys = set(sample.keys())
        if keys.issubset(set(self.column_header)) and set(self.column_header).issubset(keys):
            return True
        return False

    def AddSample(self, sample):
        if not self.ValidSample(sample):
            return False
        
        self.samples.append(sample)


    def GenerateHeader(self):
        header = ""
        header += "{}\t{}\n".format(self.MISC_OUTPUT_DIR, self.misc_output_dir)
        header += "{}\t{}\n".format(self.STRAIN_OUTPUT_DIR, self.strain_output_dir)
        header += "{}\t{}\n".format(self.PROFILE_OUTPUT_DIR, self.profile_output_dir)
        header += "{}\t{}\n".format(self.SAM_OUTPUT_DIR, self.sam_output_dir)
        return header
        
    def GenerateBody(self):
        body = '\t'.join(self.column_header) + '\n'
        
        for sample in self.samples:
            body += '\t'.join([sample[key] for key in self.column_header]) + '\n'
            
        return body
        
    def GenerateMap(self):
        return self.GenerateHeader() + self.GenerateBody()
        
        
def parse_arguments():
    parser = argparse.ArgumentParser(description="Argument Parser for Genomes, Output, and Coverages")

    # Adding arguments
    parser.add_argument('-i', '--input', type=str, help='Folder that contains all the reads')
    parser.add_argument('-o', '--output', type=str, help='Output base for protal')

    args = parser.parse_args()

    return args


args = parse_arguments()

        
read_base = args.input
output_base = args.output

all_files = get_all_files_recursive(read_base)
all_files_fastq = [f for f in all_files if f.endswith(".fq")]
all_files_fastq = [os.path.normpath(path).replace(os.path.normpath(read_base) + '/', "") for path in all_files_fastq if os.path.normpath(path).startswith(read_base)]

map = ProtalMap(read_base)
map.SetOutputsToDefault(output_base)


samples = GenerateSamplesFromReads(all_files_fastq)

for sample in samples:
    map.AddSample(sample)

print(map.GenerateMap())
