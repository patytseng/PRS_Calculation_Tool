import os
import gzip
import concurrent.futures
import threading

#generate list of all compressed VCF files
def get_vcf_file_path_list(vcf_path):
    vcf_file_path_list = []

    #grabbing all compressed .vcf.gz files and their directory
    for filename in os.listdir(vcf_path):
        if filename.startswith('ALL.chr') and filename.endswith('.vcf.gz'):
            vcf_file_path_list.append(vcf_path + filename)

    #sorting all obtained paths
    vcf_file_path_list.sort()

    #removing MT X and Y data
    return vcf_file_path_list[:-3]

def get_PRS_targets(PRS_path):
    #path to PRS file

    #generating a 22 item long list of blank lists
    PRS_targets = []
    for numerical_counter in range(0,22):
        PRS_targets.append([])

    #opening file
    with open(PRS_path,'r') as PRS_file:

        #running through file
        for line in PRS_file:

            #filtering for relevant lines
            if line[0:2] == 'rs':

                #converting text string of each line to a list
                line = list(line.split('\t'))[1:]

                #taking the index line out cuz string can't indicate index
                if line[0] == 'chr_name':
                    pass

                #checking if the chromosome number is a number
                elif line[0].isnumeric:
                    PRS_targets[int(line[0])-1].append(line[1:6])

                #for debugging just in case something goes wrong
                else:
                    print('FLAG: PRS no chromosome index')
    return PRS_targets

#obtaining subject ID list
def get_ID_library(vcf_file_path_list):
    #pulling relevant line from vcf file header
    with gzip.GzipFile(vcf_file_path_list[0],'r') as decompressed_file:
        linecount = 0
        for line in decompressed_file:
            linecount += 1
            line = str(line)[2:-3]
            if line[0:6] == '#CHROM':
                data_index = list(line.split('\\t'))
                IDs = data_index[9:]
                break

    #converting to dictionary for easier processing
    ID_library = {}
    for ID in IDs:
        ID_library[ID] = 0

    return ID_library

def get_partial_PRS(vcf_file_path,chromosome_number,PRS_targets,ID_library):
    #opening vcf file
    with gzip.GzipFile(vcf_file_path,'r') as decompressed_file:
        linecount = 0
        PRS_target_index = 0
        for line in decompressed_file:
            linecount += 1
            line = str(line)[2:-3]
            if line[0:6] == '#CHROM':
                data_index = list(line.split('\\t'))
                IDs = data_index[9:]
            elif line[0].isnumeric and line[0] != '#':
                line = list(line.split('\\t'))
                if line[1] == PRS_targets[PRS_target_index][0]:
                    if line[3] == PRS_targets[PRS_target_index][2] and line[4] == PRS_targets[PRS_target_index][1] and line[8] == 'GT':
                        line = line[9:]

                        GT_index = 0
                        for GT in line:
                            line[GT_index] = sum(list(map(int, line[GT_index].split('|'))))
                            GT_index += 1

                        GT_index = 0
                        for ID in IDs:
                            ID_library[ID] = ID_library[ID] + (float(PRS_targets[PRS_target_index][3])*line[GT_index])
                            GT_index += 1

                    PRS_target_index += 1
                    print(f'CHR{chromosome_number} @ {round(PRS_target_index/len(PRS_targets)*100,3)}%')
                    
                    #break condition
                    if PRS_target_index == len(PRS_targets):
                        break
        return ID_library
    
def task(vcf_file_path,chromosome_number,PRS_targets,ID_library):
    print(f'starting CHR{chromosome_number}')
    output = get_partial_PRS(vcf_file_path,chromosome_number,PRS_targets,ID_library)
    print(f'finishing CHR{chromosome_number}')
    return output

def main():
    #collecting basic variables
    vcf_path = '/Users/patytseng/bin/1000Genomes/'
    PRS_path = '/Users/patytseng/Desktop/Labs/Torkamani_Lab/20210413_Project_1_BMIPRSFor1000GenomesProject/PGS000299.txt'
    vcf_file_path_list = get_vcf_file_path_list(vcf_path)
    PRS_targets = get_PRS_targets(PRS_path)
    ID_library = get_ID_library(vcf_file_path_list)
    #multithreading over all the chromosomes
    chromosome_count = 22
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        future_list = []
        for i in range(0,chromosome_count):
            future = executor.submit(task,vcf_file_path_list[i],i+1,PRS_targets[i],ID_library)
            future_list.append(future)
        PRS_ID_library_list = []
        for i in range(0,chromosome_count):
            PRS_ID_library_list.append(future_list[i].result())
    for key in ID_library:
        for i in range(0,chromosome_count):
            ID_library[key] = round(ID_library[key] + PRS_ID_library_list[i][key],4)
    #opening file to write to
    try:
        output_file = open('PRS_output.txt','x')
    except:
        output_file = open('PRS_output.txt','w')
    #writing to file: header then data
    output_file.write('ID\tScore\n')
    for key in ID_library:
        output_file.write(f'{key}\t{ID_library[key]}\n')

if __name__ == '__main__':
    main()
