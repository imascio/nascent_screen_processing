import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import pickle
import re
from Bio.Seq import Seq

    
def UMI_attach_read2_barcode_list(sample, input_folder, output_folder, ligation_barcode_list, RT_barcode_list, randomN_barcodes):
    #open the read1, read2, and output file
    Read1 = input_folder + "/" + sample + ".R1.fastq.gz"
    Read2 = input_folder + "/" + sample + ".R3.fastq.gz"
    Read3 = input_folder + "/" + sample + ".R2.fastq.gz"
    output_file1 = output_folder + "/" + sample + ".R2.fastq.gz"
    output_file2 = output_folder + "/" + sample + ".R1.fastq.gz"
    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    f3 = gzip.open(output_file1, 'wb')
    f4 = gzip.open(Read3)
    f5 = gzip.open(output_file2, 'wb')
    
    line1 = f1.readline()
    line2 = f2.readline()
    line3 = f4.readline()
    total_line = 0
    filtered_line = 0
    
    while (line1):
        total_line += 1
        # read in sequence line from R1 fastq (UMI + RT barcode)
        line1 = f1.readline()
        # read in sequence line from i5 (R2 fastq, ligation barcode)
        line3 = f4.readline()
        # first check if the ligation barcode match with the barcode
        tmp_lig = line3[0:10]
        # convert to a string
        tmp_lig_str = tmp_lig.decode()

        if tmp_lig_str in ligation_barcode_list:
            #print('Yay!')
            ligation_bc_match_str = ligation_barcode_list[tmp_lig_str]
            ligation_bc_match = bytes(ligation_bc_match_str, 'utf-8')
            # check RT barcode
            target_RT = line1[8:18]
            #convert to a string
            target_RT_str = target_RT.decode()
            
            # check if the RT barcode is a match and correct within 1 editing distance
            if target_RT_str in RT_barcode_list:
                barcode_str = RT_barcode_list[target_RT_str]
                #convert back to bytes
                barcode = bytes(barcode_str, 'utf-8')
                UMI = line1[:8]
                #convert UMI to a string
                UMI_str = UMI.decode()

                #check if the RT barcode is a random barcode
                if barcode in randomN_barcodes:
                    
                    fq_header_split = line2.split(b" ")
                    #print(line2)
                    #print(fq_header_split)
                    #print(fq_header_split[0])
                    first_line_3 = fq_header_split[0] + b'_' + ligation_bc_match + barcode + b'_' + UMI + b' ' + fq_header_split[1]
                    #print(first_line_3)
                    first_line_5 = fq_header_split[0] + b'_' + ligation_bc_match + barcode + b'_' + UMI + b' ' + fq_header_split[1]
                    #first_line_3 = '@' + ligation_bc_match_str + barcode + ',' + UMI_str + ',' + re.sub("3:N:","2:N:",line2[1:])
                    #first_line_5 = '@' + ligation_bc_match + barcode + ',' + UMI + ',' + re.sub("3:N:","1:N:",line2[1:])

                    second_line_5 = line1[18:]                
                    result_adaptor = re.search(b"CTGTCTCTTATACACAT", second_line_5)
                    if result_adaptor == None:
                        second_line_5 = second_line_5
                    else:
                        second_line_5 = second_line_5[:result_adaptor.start()] + b"\n"

                    seq = Seq(target_RT)
                    seqRC = str(seq.reverse_complement())
                    second_line_3 = f2.readline()
                    result_barcode = re.search(bytes(seqRC, 'utf-8'), second_line_3)
                    if result_barcode == None:
                        second_line_3 = second_line_3
                    else:
                        second_line_3 = second_line_3[:result_barcode.start()] + b"\n"

                    third_line_3 = f2.readline()
                    third_line_5 = f1.readline()

                    four_line_5 = f1.readline()
                    four_line_5 = four_line_5[18:]
                    if result_adaptor == None:
                        four_line_5 = four_line_5
                    else:
                        four_line_5 = four_line_5[:result_adaptor.start()] + b"\n"

                    four_line_3 = f2.readline()
                    if result_barcode == None:
                        four_line_3 = four_line_3
                    else:
       	       	       	four_line_3 = four_line_3[:result_barcode.start()] + b"\n"

                    line1 = f1.readline()
                    line2 = f2.readline()

                    if len(second_line_3) > 20 and len(second_line_5) > 0 :
                        filtered_line += 1
                        f3.write(first_line_3)
                        f5.write(first_line_5)
       	       	       	f3.write(second_line_3)
       	       	       	f5.write(second_line_5)
       	       	       	f3.write(third_line_3)
       	       	       	f5.write(third_line_5)
       	       	       	f3.write(four_line_3)
       	       	       	f5.write(four_line_5) 

                else:            
                    # add the ligation BC, RT BC, and UMI to the fastq header in a UMI tools commpatible format
                    fq_header_split = line2.split(b" ")
                    #print(line2)
                    #print(fq_header_split)
                    #print(fq_header_split[0])
                    first_line_3 = fq_header_split[0] + b'_' + ligation_bc_match + barcode + b'_' + UMI + b' ' + fq_header_split[1]
                    #print(first_line_3)
                    first_line_5 = fq_header_split[0] + b'_' + ligation_bc_match + barcode + b'_' + UMI + b' ' + fq_header_split[1]
                    #first_line_3 = '@' + ligation_bc_match_str + barcode + ',' + UMI_str + ',' + re.sub("3:N:","2:N:",line2[1:])
                    #first_line_5 = '@' + ligation_bc_match + barcode + ',' + UMI + ',' + re.sub("3:N:","1:N:",line2[1:])

                    
                    # remove nextera read 2 adaptor from read 1(?)
                    second_line_5 = line1[5:]
                    result_adaptor = re.search(b"CTGTCTCTTATACACAT", second_line_5)
                    if result_adaptor == None:
                        #print('No adaptor')
                        second_line_5 = second_line_5
                    else:
                        #print('Found adaptor')
                        second_line_5 = second_line_5[:result_adaptor.start()] + "\n"

                    
                    # check for reverse complement of RT barcode in read 2 (actual read 2 of gene body)
                    seq = Seq(target_RT)
                    seqRC = str(seq.reverse_complement())
                    #seqRC = seq.reverse_complement()
                    #print(seq)
                    #print(seqRC)
                    second_line_3 = f2.readline()
                    result_barcode = re.search(bytes(seqRC, 'utf-8'), second_line_3)
                    if result_barcode == None:
                        #print('No RT barcode in read 2')
                        second_line_3 = second_line_3
                    else:
                        #print('Found RT barcode in read 2')
                        second_line_3 = second_line_3[:(result_barcode.start() - 15)] + b"\n"

                    # read in third fastq line of read 2 and read 1
                    third_line_3 = f2.readline()
                    third_line_5 = f1.readline()

                    # removed adaptor from the fourth fastq line (quality score)
                    four_line_5 = f1.readline()
                    four_line_5 = four_line_5[5:]
                    if result_adaptor == None:
                        four_line_5 = four_line_5
                    else:
                        four_line_5 = four_line_5[:result_adaptor.start()] + b"\n"

                    four_line_3 = f2.readline()
                    if result_barcode == None:
                        four_line_3 = four_line_3
                    else:
                        four_line_3 = four_line_3[:(result_barcode.start() - 15)] + b"\n"

                    line1 = f1.readline()
                    line2 = f2.readline()
                    # write fastq file
                    if len(second_line_3) > 20 and len(second_line_5) > 0 :
                        filtered_line += 1
                        f3.write(first_line_3)
                        f5.write(first_line_5)
                        f3.write(second_line_3)
                        f5.write(second_line_5)
                        f3.write(third_line_3)
                        f5.write(third_line_5)
                        f3.write(four_line_3)
                        f5.write(four_line_5)


            else:
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()

                
        else:
            #print('Uh Oh!')
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
       	    line1 = f1.readline()
            line1 = f1.readline()
            line1 = f1.readline()


        line3 = f4.readline() 
        line3 = f4.readline()
        line3 = f4.readline()

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f" 
          %(sample, total_line, filtered_line, float(filtered_line) / float(total_line)))

# this function accept an input folder and a output folder and then generate the output file with the index
def attach_UMI_files(input_folder, sample, output_folder, ligation_barcode_file, RT_barcode_file, core, randomN_barcode_file):
    print(f'''
    --------------------------start attaching UMI-----------------------------
    Sample ID: {sample}
    Input folder: {input_folder}
    Output folder: {output_folder}
    Ligation barcode file: {ligation_barcode_file}
    RT barcode file: {RT_barcode_file}
    ___________________________________________________________________________
    ''')

    print("Load ligation barcode dictionary...")
    with open(ligation_barcode_file, 'rb') as f:
        ligation_barcode_list = pickle.load(f)

    print("Load RT barcode dictionary...")
    with open(RT_barcode_file, 'rb') as f:
        RT_barcode_list = pickle.load(f)

    print("Load randomN barcode list...")
    randomN_barcodes = []
    with open(randomN_barcode_file, 'rb') as f:
        for line in f:
            randomN_barcodes.append(line.strip())

    # Process this single sample
    UMI_attach_read2_barcode_list(sample, input_folder, output_folder, ligation_barcode_list, RT_barcode_list, randomN_barcodes)

    print("~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~")
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sample = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode_file = sys.argv[4]
    RT_barcode_file = sys.argv[5]
    core = sys.argv[6]
    randomN_barcode_file = sys.argv[7]
    attach_UMI_files(input_folder, sample, output_folder, ligation_barcode_file, RT_barcode_file, core, randomN_barcode_file)
