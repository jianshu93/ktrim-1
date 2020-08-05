import sys
import os
import subprocess
import tempfile
import argparse
import multiprocessing
import re

#This contains code which generates a complete list of illumina adapters from scratch
def generate_adapters_temporary_file():

	print("Preparing adapter file for you.")
	adapters_dict = {}
	
	'''
	I identify the adapter families here with comments. Any adapter recognized in one of these during preprocessing will include 
	all of the members of its family in final, e.g. seeing Illumina_Single_End_Apapter_1 will include the following:
	Illumina_Single_End_Apapter_1, Illumina_Single_End_Apapter_2, Illumina_Single_End_PCR_Primer_1, Illumina_Single_End_PCR_Primer_2, and Illumina_Single_End_Sequencing_Primer
	in the final filtering fasta
	'''
	
	#Single end family
	adapters_dict[">Illumina_Single_End_Apapter_1"] = "ACACTCTTTCCCTACACGACGCTGTTCCATCT"
	adapters_dict[">Illumina_Single_End_Apapter_2"] = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
	adapters_dict[">Illumina_Single_End_PCR_Primer_1"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Single_End_PCR_Primer_2"] = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
	adapters_dict[">Illumina_Single_End_Sequencing_Primer"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	
	#Paired end family
	adapters_dict[">Illumina_Paired_End_Adapter_1"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Paired_End_Adapter_2"] = "CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	adapters_dict[">Illumina_Paried_End_PCR_Primer_1"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Paired_End_PCR_Primer_2"] = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	adapters_dict[">Illumina_Paried_End_Sequencing_Primer_1"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Paired_End_Sequencing_Primer_2"] = "CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT"
	
	#DpnII family
	adapters_dict[">Illumina_DpnII_expression_Adapter_1"] = "ACAGGTTCAGAGTTCTACAGTCCGAC"
	adapters_dict[">Illumina_DpnII_expression_Adapter_2"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_DpnII_expression_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_DpnII_expression_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict[">Illumina_DpnII_expression_Sequencing_Primer"] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict[">Illumina_NlaIII_expression_Adapter_1"] = "ACAGGTTCAGAGTTCTACAGTCCGACATG"
	adapters_dict[">Illumina_NlaIII_expression_Adapter_2"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_NlaIII_expression_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_NlaIII_expression_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict[">Illumina_NlaIII_expression_Sequencing_Primer"] = "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"
	
	#Small RNA family
	adapters_dict[">Illumina_Small_RNA_Adapter_1"] = "GTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict[">Illumina_Small_RNA_Adapter_2"] = "TCGTATGCCGTCTTCTGCTTGT"
	adapters_dict[">Illumina_Small_RNA_RT_Primer"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_Small_RNA_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_Small_RNA_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict[">Illumina_Small_RNA_Sequencing_Primer"] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	
	
	#Multiplexing Family
	adapters_dict[">Illumina_Multiplexing_Adapter_1"] = "GATCGGAAGAGCACACGTCT"
	adapters_dict[">Illumina_Multiplexing_Adapter_2"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Multiplexing_PCR_Primer_1.01"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Multiplexing_PCR_Primer_2.01"] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
	adapters_dict[">Illumina_Multiplexing_Read1_Sequencing_Primer"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">Illumina_Multiplexing_Index_Sequencing_Primer"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	adapters_dict[">Illumina_Multiplexing_Read2_Sequencing_Primer"] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
	
	
	#PCR primer family
	adapters_dict[">Illumina_PCR_Primer_Index_1"] = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_2"] = "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_3"] = "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_4"] = "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_5"] = "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_6"] = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_7"] = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_8"] = "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_9"] = "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_10"] = "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_11"] = "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"
	adapters_dict[">Illumina_PCR_Primer_Index_12"] = "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"
	
	
	#DpnII Gex family
	adapters_dict[">Illumina_DpnII_Gex_Adapter_1"] = "GATCGTCGGACTGTAGAACTCTGAAC"
	adapters_dict[">Illumina_DpnII_Gex_Adapter_1.01"] = "ACAGGTTCAGAGTTCTACAGTCCGAC"
	adapters_dict[">Illumina_DpnII_Gex_Adapter_2"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_DpnII_Gex_Adapter_2.01"] = "TCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">Illumina_DpnII_Gex_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_DpnII_Gex_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict[">Illumina_DpnII_Gex_Sequencing_Primer"] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict[">Illumina_NlaIII_Gex_Adapter_1.01"] = "TCGGACTGTAGAACTCTGAAC"
	adapters_dict[">Illumina_NlaIII_Gex_Adapter_1.02"] = "ACAGGTTCAGAGTTCTACAGTCCGACATG"
	adapters_dict[">Illumina_NlaIII_Gex_Adapter_2.01"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_NlaIII_Gex_Adapter_2.02"] = "TCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">Illumina_NlaIII_Gex_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict[">Illumina_NlaIII_Gex_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict[">Illumina_NlaIII_Gex_Sequencing_Primer"] = "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"
	
	#Other RNA family
	adapters_dict[">Illumina_5p_RNA_Adapter"] = "GTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict[">Illumina_RNA_Adapter1"] = "TCGTATGCCGTCTTCTGCTTGT"
	adapters_dict[">Illumina_Small_RNA_3p_Adapter_1"] = "ATCTCGTATGCCGTCTTCTGCTTG"
	
	#TrueSeq family
	adapters_dict[">TruSeq_Universal_Adapter"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict[">TruSeq_Adapter_Index_1"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_2"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_3"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_4"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_5"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_6"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_7"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_8"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_9"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_10"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_11"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq_Adapter_Index_12"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"
	
	#RNA PCR family
	adapters_dict[">Illumina_RNA_RT_Primer"] = "GCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">Illumina_RNA_PCR_Primer"] = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
	adapters_dict[">RNA_PCR_Primer_Index_1"] = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_2"] = "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_3"] = "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_4"] = "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_5"] = "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_6"] = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_7"] = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_8"] = "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_9"] = "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_10"] = "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_11"] = "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_12"] = "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_13"] = "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_14"] = "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_15"] = "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_16"] = "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_17"] = "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_18"] = "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_19"] = "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_20"] = "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_21"] = "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_22"] = "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_23"] = "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_24"] = "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_25"] = "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_26"] = "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_27"] = "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_28"] = "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_29"] = "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_30"] = "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_31"] = "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_32"] = "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_33"] = "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_34"] = "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_35"] = "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_36"] = "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_37"] = "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_38"] = "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_39"] = "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_40"] = "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_41"] = "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_42"] = "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_43"] = "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_44"] = "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_45"] = "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_46"] = "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_47"] = "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict[">RNA_PCR_Primer_Index_48"] = "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	
	#ABI family
	adapters_dict[">ABI_Dynabead_EcoP_Oligo"] = "CTGATCTAGAGGTACCGGATCCCAGCAGT"
	adapters_dict[">ABI_Solid3_Adapter_A"] = "CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"
	adapters_dict[">ABI_Solid3_Adapter_B"] = "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"
	adapters_dict[">ABI_Solid3_5_AMP_Primer"] = "CCACTACGCCTCCGCTTTCCTCTCTATG"
	adapters_dict[">ABI_Solid3_3_AMP_Primer"] = "CTGCCCCGGGTTCCTCATTCT"
	adapters_dict[">ABI_Solid3_EF1_alpha_Sense_Primer"] = "CATGTGTGTTGAGAGCTTC"
	adapters_dict[">ABI_Solid3_EF1_alpha_Antisense_Primer"] = "GAAAACCAAAGTGGTCCAC"
	adapters_dict[">ABI_Solid3_GAPDH_Forward_Primer"] = "TTAGCACCCCTGGCCAAGG"
	adapters_dict[">ABI_Solid3_GAPDH_Reverse_Primer"] = "CTTACTCCTTGGAGGCCATG"
	
	#TrueSeq2 family
	adapters_dict[">TruSeq2_SE"] = "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict[">TruSeq2_PE_f"] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	adapters_dict[">TruSeq2_PE_r"] = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"
	adapters_dict[">TruSeq3_IndexedAdapter"] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	adapters_dict[">TruSeq3_UniversalAdapter"] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
	
	#Nextera Family
	adapters_dict[">Nextera_PE_PrefixNX/1"] = "AGATGTGTATAAGAGACAG"
	adapters_dict[">Nextera_PE_PrefixNX/2"] = "AGATGTGTATAAGAGACAG"
	adapters_dict[">Nextera_PE_Trans1"] = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
	adapters_dict[">Nextera_PE_Trans1_rc"] = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
	adapters_dict[">Nextera_PE_Trans2"] = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
	adapters_dict[">Nextera_PE_Trans2_rc"] = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
	
	all_adapters = tempfile.NamedTemporaryFile(mode = "w", delete = False)
		
	for adapt in adapters_dict:
		print(adapt, file = all_adapters)
		print(adapters_dict[adapt], file = all_adapters)
	
	name = all_adapters.name
	all_adapters.close()
	
	return adapters_dict, name

#identifies the same adapters as in the full file with a family of origin, so that all adapters in a family can be selected.
def create_adapter_families():
	adapters_fam_dict = {}
	
	adapters_fam_dict["Illumina_Single_End_Apapter_1"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_Apapter_2"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_PCR_Primer_1"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_PCR_Primer_2"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_Sequencing_Primer"] = "singleend"
	adapters_fam_dict["Illumina_Paired_End_Adapter_1"] = "pairedend"
	adapters_fam_dict["Illumina_Paired_End_Adapter_2"] = "pairedend"
	adapters_fam_dict["Illumina_Paried_End_PCR_Primer_1"] = "pairedend"
	adapters_fam_dict["Illumina_Paired_End_PCR_Primer_2"] = "pairedend"
	adapters_fam_dict["Illumina_Paried_End_Sequencing_Primer_1"] = "pairedend"
	adapters_fam_dict["Illumina_Paired_End_Sequencing_Primer_2"] = "pairedend"
	adapters_fam_dict["Illumina_DpnII_expression_Adapter_1"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_Adapter_2"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_PCR_Primer_1"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_PCR_Primer_2"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_Sequencing_Primer"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_Adapter_1"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_Adapter_2"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_PCR_Primer_1"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_PCR_Primer_2"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_Sequencing_Primer"] = "dpnII"
	adapters_fam_dict["Illumina_Small_RNA_Adapter_1"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_Adapter_2"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_RT_Primer"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_PCR_Primer_1"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_PCR_Primer_2"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_Sequencing_Primer"] = "smallrna"
	adapters_fam_dict["Illumina_Multiplexing_Adapter_1"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Adapter_2"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_PCR_Primer_1.01"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_PCR_Primer_2.01"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Read1_Sequencing_Primer"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Index_Sequencing_Primer"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Read2_Sequencing_Primer"] = "multiplex"
	adapters_fam_dict["Illumina_PCR_Primer_Index_1"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_2"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_3"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_4"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_5"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_6"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_7"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_8"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_9"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_10"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_11"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_12"] = "pcr"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_1"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_1.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_2"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_2.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_PCR_Primer_1"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_PCR_Primer_2"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Sequencing_Primer"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_1.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_1.02"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_2.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_2.02"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_PCR_Primer_1"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_PCR_Primer_2"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Sequencing_Primer"] = "dpnIIgex"
	adapters_fam_dict["Illumina_5p_RNA_Adapter"] = "otherrna"
	adapters_fam_dict["Illumina_RNA_Adapter1"] = "otherrna"
	adapters_fam_dict["Illumina_Small_RNA_3p_Adapter_1"] = "otherrna"
	adapters_fam_dict["TruSeq_Universal_Adapter"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_1"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_2"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_3"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_4"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_5"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_6"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_7"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_8"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_9"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_10"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_11"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_12"] = "trueseq"
	adapters_fam_dict["Illumina_RNA_RT_Primer"] = "rnapcr"
	adapters_fam_dict["Illumina_RNA_PCR_Primer"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_1"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_2"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_3"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_4"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_5"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_6"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_7"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_8"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_9"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_10"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_11"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_12"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_13"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_14"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_15"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_16"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_17"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_18"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_19"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_20"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_21"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_22"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_23"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_24"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_25"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_26"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_27"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_28"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_29"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_30"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_31"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_32"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_33"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_34"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_35"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_36"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_37"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_38"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_39"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_40"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_41"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_42"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_43"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_44"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_45"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_46"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_47"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_48"] = "rnapcr"
	adapters_fam_dict["ABI_Dynabead_EcoP_Oligo"] = "abi"
	adapters_fam_dict["ABI_Solid3_Adapter_A"] = "abi"
	adapters_fam_dict["ABI_Solid3_Adapter_B"] = "abi"
	adapters_fam_dict["ABI_Solid3_5_AMP_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_3_AMP_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_EF1_alpha_Sense_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_EF1_alpha_Antisense_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_GAPDH_Forward_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_GAPDH_Reverse_Primer"] = "abi"
	adapters_fam_dict["TruSeq2_SE"] = "trueseq2"
	adapters_fam_dict["TruSeq2_PE_f"] = "trueseq2"
	adapters_fam_dict["TruSeq2_PE_r"] = "trueseq2"
	adapters_fam_dict["TruSeq3_IndexedAdapter"] = "trueseq2"
	adapters_fam_dict["TruSeq3_UniversalAdapter"] = "trueseq2"
	adapters_fam_dict["Nextera_PE_PrefixNX/1"] = "nextera"
	adapters_fam_dict["Nextera_PE_PrefixNX/2"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans1"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans1_rc"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans2"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans2_rc"] = "nextera"
	
	return(adapters_fam_dict)
	
#Checks the command. SeqTK needs piped, so it has to have some special handling from subprocess.
def seqtk_manager(command):
	if command[0] == "fastqc":
		print("Running FastQC for", command[4])
		subprocess.call(command)
		return "fastqc"
	else:
		print("Subsampling", command[4], "with SeqTK")
		temp = tempfile.NamedTemporaryFile("w", delete=False)
		#Handle the wierd issues with piping.
		ps = subprocess.run(command, stdout=subprocess.PIPE, universal_newlines = True)
		temp.write(ps.stdout)
		
		name = temp.name
		temp.close()
				
		return name

'''
Complex function. This does FastQC pre-trim reports, subsamples reads using SeqTK, uses FaQCs on the subsamples, 
identifies the detected adapters, based on the FaQCs report and returns the set of detected adapters for use in
creating a subsetted fasta for use in the full trim.
'''
def identify_adapters(artificial_artifacts, forward = "", reverse = "", unpaired = "", complete = "", threads = 1, output = ".", minimum_presence = 0.75):
	print("Performing pre-trim QC with FastQC and checking for specific adapter presence with FaQCs...")
		
	f_check = False
	r_check = False
	u_check = False
	
	count = 0
	
	if forward != "":
		f_check = True
		count += 1
	if reverse != "":
		r_check = True
		count += 1
	if unpaired != "":
		u_check = True
		count += 1
	
	fastqc_command_f = ["fastqc", "--quiet", "-o", output, forward]
	fastqc_command_r = ["fastqc", "--quiet", "-o", output, reverse]
	fastqc_command_u = ["fastqc", "--quiet", "-o", output, unpaired]
	
	seqtk_command_f = ["seqtk", "sample", "-s", "100", forward, "100000"]
	seqtk_command_r = ["seqtk", "sample", "-s", "100", reverse, "100000"]
	seqtk_command_u = ["seqtk", "sample", "-s", "100", unpaired, "100000"]
		
	overall_command = []
		
	if f_check:
		overall_command.append(fastqc_command_f)
		overall_command.append(seqtk_command_f)
	if r_check:
		overall_command.append(fastqc_command_r)
		overall_command.append(seqtk_command_r)
	if u_check:
		overall_command.append(fastqc_command_u)
		overall_command.append(seqtk_command_u)
	
	p = multiprocessing.Pool(min(threads, count*2))
	set = p.map(seqtk_manager, overall_command)
	
	p.close()
	
	clean_set = []
	
	for item in set:
		if item != "fastqc":
			clean_set.append(item)
	
	#Has to be a string because subprocess won't take it otherwise
	faqcs_subset_command = ["FaQCs", "-t", str(threads), "--qc_only", "-d", output, "--artifactFile", artificial_artifacts, "--prefix", "Subsample_Adapter_Detection"]
	
	if f_check:
		faqcs_subset_command.append("-1")
		faqcs_subset_command.append(clean_set[0])
	if r_check:
		faqcs_subset_command.append("-2")
		faqcs_subset_command.append(clean_set[1])
	if u_check:
		faqcs_subset_command.append("-u")
		#If this is the only one supplied, it's 0, if it's part of a set of 3, it's 2
		faqcs_subset_command.append(clean_set[len(clean_set)-1])
		

	#Runs the faqcs on however many threads it has access to.
	print("Detecting adapters now...")
	ps = subprocess.Popen(faqcs_subset_command)
	ps.wait()
	
	detection_report = open(output + "/" + "Subsample_Adapter_Detection.stats.txt")
	
	detected_adapters = {}
	begin_assessment = False
	for line in detection_report:
		if not begin_assessment:
			if line.strip() == "Reads with Adapters/Primers: 23183 (11.59 %)":
				begin_assessment = True
		else:
			segment = line.strip().split()
			detected_adapters[segment[0]] = float(re.findall("\d+\.\d+", segment[3])[0])
	
	detection_report.close()
	
	clean_detection = []
	
	for adapter in detected_adapters:
		if detected_adapters[adapter] >= minimum_presence:
			clean_detection.append(adapter)
	
	#Cleans up after itself.
	os.remove(artificial_artifacts)
	for item in clean_set:
		os.remove(item)
	
	print("Detection done!")
	
	return clean_detection
	
#Figures out which adapter families to include in the final trim and generates a fasta on that basis
def parse_adapters(full_list, detected_adapters, output):
	print("Creating specific adapters file for you.")
	
	#This should be a dict of groups of adapters whose set presence is identified through the presence any of the members
	recognized_adapter_collections = create_adapter_families()
	
	#This is a file I don't want to be temporary. It both helps identify the adapters present in a dataset and provides a fasta for a user to reuse
	subset = open(output+"/"+"detected_adapters.fasta", "w")
	
	families_detected = []
	
	#Needs to actually check the things
	for adapter in detected_adapters:
		if adapter in recognized_adapter_collections:
			if recognized_adapter_collections[adapter] not in families_detected:
				families_detected.append(recognized_adapter_collections[adapter])
	
	
	finalized_adapters = []
	
	for adapter in recognized_adapter_collections:
		if recognized_adapter_collections[adapter] in families_detected:
			finalized_adapters.append(">"+adapter)
	
	for adapter in full_list:
		if adapter in finalized_adapters:
			print(adapter, file = subset)
			print(full_list[adapter], file = subset)
			
	subset.close()
	
	return(output+"/"+"detected_adapters.fasta")

#This should handle it when the artifacts file is empty
def full_trim(subsetted_artifacts, forward = "", reverse = "", unpaired = "", complete = "", threads = 1, output = ".", prefix = "fully_trimmed_reads"):
	print("Performing the full trim of your reads. This will (probably) take a while.")
	
	f_check = False
	r_check = False
	u_check = False
	
	count = 0
	
	if forward != "":
		f_check = True
		count += 1
	if reverse != "":
		r_check = True
		count += 1
	if unpaired != "":
		u_check = True
		count += 1
		
	faqcs_command = ["FaQCs", "--mode", "HARD", "-t", str(threads), "-d", output, "--artifactFile", subsetted_artifacts]
	
	if prefix != "":
		faqcs_command.append("--prefix")
		faqcs_command.append(prefix)
	
	if f_check:
		faqcs_command.append("-1")	
		faqcs_command.append(forward)
	if r_check:
		faqcs_command.append("-2")
		faqcs_command.append(reverse)
	if u_check:
		faqcs_command.append("-u")
		#If this is the only one supplied, it's 0, if it's part of a set of 3, it's 2
		faqcs_command.append(unpaired)
		

	#Runs the faqcs on however many threads it has access to.
	ps = subprocess.Popen(faqcs_command)
	ps.wait()
	
	print("Reads Trimmed!")
	
	return None
	
#One last FastQC on the trimmed reads
def final_QC(forward = "", reverse = "", unpaired = "", threads = 1, output = ".", prefix = "fully_trimmed_reads"):

	print("Performing FastQC on trimmed reads...")
	f_check = False
	r_check = False
	u_check = False
	
	count = 0
		
	if forward != "":
		f_check = True
		forward = output+"/"+prefix+".1.trimmed.fastq"
		count += 1
	if reverse != "":
		r_check = True
		reverse = output+"/"+prefix+".2.trimmed.fastq"
		count += 1
	if unpaired != "":
		u_check = True
		unpaired = output+"/"+prefix+".unpaired.trimmed.fastq"
		count += 1
	
	fastqc_command_f = ["fastqc", "--quiet", "-o", output, forward]
	fastqc_command_r = ["fastqc", "--quiet", "-o", output, reverse]
	fastqc_command_u = ["fastqc", "--quiet", "-o", output, unpaired]
	
	overall_command = []
		
	if f_check:
		overall_command.append(fastqc_command_f)
	if r_check:
		overall_command.append(fastqc_command_r)
	if u_check:
		overall_command.append(fastqc_command_u)
	
	p = multiprocessing.Pool(min(threads, count))
	set = p.map(seqtk_manager, overall_command)
	
	p.close()
	
	return None
	
#Stolen from a SO thread on how to issue usage information on an error.
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def gather_opts():
	parser = MyParser(description=''' This program is designed to facilitate effective trimming of your reads.
	It will help to identify the presence of adapters in your reads, trim those adapters and the reads efficiently,
	and produce several bfore and after quality reports in addition to the trimmed reads. This is a pipeline incorporating
	FaQCs, FastQC, and seqtk commands, in addition to several python operations which exist to facilitate adapter finding and
	subsetting.''')
	#Use all available cores.
	parser.add_argument("--UNLIMITED_POWER", dest = "Sheev", action = 'store_true', help = "Have you ever heard the tragedy of Darth Parallelegius?")
	
	parser.add_argument("--threads", "-t", dest = "threads", default = 1, help = "How many threads do you want to use? Use --UNLIMITED_POWER flag to use all of them.")

	parser.add_argument("--forward", "-1", dest = "f", default = "", help = "Forward Strand Reads (use -u for unpaired reads)")
	parser.add_argument("--reverse", "-2", dest = "r", default = "", help = "Reverse Strand Reads (use -u for unpaired reads)")
	parser.add_argument("--unpaired", "-u", dest = "u", default = "", help = "Unpaired Reads")
	
	parser.add_argument("--output", "-o", dest = "outdir", default = ".", help = "Directory to send final outputs. Must exist.")
	
	parser.add_argument("--min_adapt_pres", "-m", dest = "minpres", default = 0.75, help = "Minimum presence of an adapter for it to be considered present in a set of reads.")
	
	
	return(parser, parser.parse_args())

#Program Control
def main():
	#Keep the parser on hand so I can prent usage as needed.
	help_message, options = gather_opts()
	
	if len(sys.argv)==1:
		help_message.print_help(sys.stderr)
		sys.exit(1)
		
	final_output = options.outdir
	
	#Get the reads
	f = options.f
	r = options.r
	u = options.u
	
	
	#Check to make sure it actually has data, or exits
	if f == "" and r == "" and u == "":
		print("I need to be given reads! Exiting program.")
		help_message.print_help(sys.stderr)
		quit()
	
	if f == "" and r != "" or f != "" and r == "":
		print("If you have paired reads, I need both the forward and reverse files. If you just want to process one, use -u to specify it. Exiting program.")
		help_message.print_help(sys.stderr)
		quit()
	
	threads = options.threads
	
	minpres = float(options.minpres)
	
	#It's a joke, get off my back!
	#Uses all the threads a system has available. This is written specifically to permit function in a child process on a unix server in which
	#there are fewer threads/processors available to this program than exist on the system.
	if options.Sheev:
		print("\nI AM THE SENATE!\n")
		threads = len(os.sched_getaffinity(0))
		
	#Tell the user how we're doing and what we're working on.
	print("Working with", threads, "threads.")
	print("Placing results in:", final_output)
	print("Primary Strand Reads:", f, "\nReverse Strand Reads:", r, "\nUnpaired Reads:", u)
	
	print("Adapters considered detected if present in "+ str(minpres) + " % of reads.")

	
	#necessary preparation - creates all adapters file
	adapter_set, complete_adapter_file_name = generate_adapters_temporary_file()

	#Perform pre QC, detect present adapters
	correct_adapters = identify_adapters(forward = f, reverse = r, unpaired = u, threads = threads, artificial_artifacts = complete_adapter_file_name, output = final_output, minimum_presence = minpres)

	#Generate fasta with all of the detected adapters
	cleaned_adapters = parse_adapters(adapter_set, correct_adapters, final_output)
	
	#Trim according to the adapter set, parameters
	full_trim(subsetted_artifacts = cleaned_adapters, forward = f, reverse = r, unpaired = u, threads = threads, output = final_output)
	
	#one more FastQC on the final outputs
	final_QC(forward = f, reverse = r, unpaired = u, threads = threads, output = final_output)
	
	print("IT'S OVER, [Insert User Name], I HAVE THE HIGH GROUND!")
	
	
	
	
if __name__ == "__main__":
	main()

	

	
#Leftover creation functions.	
	
#Regenerate adatapers file output from a fasta. This is a utility function I do not expect to see used in the final product.	
def fasta_to_permanent_python(original_adapters_fasta):
	fasta = open(original_adapters_fasta, "r")
	
	fasta_seq_dict = {}
	
	current_line = fasta.readline().strip()
	
	current_id = current_line
	current_seq = ""
	
	current_line = fasta.readline().strip()

	while current_line:
		if current_line.startswith(">"):
			fasta_seq_dict[current_id] = current_seq
			current_id = current_line
			current_seq = ""
		else:
			current_seq += current_line
		
		current_line = fasta.readline().strip()
		
	#Finally, python needs this logic.
	fasta_seq_dict[current_id] = current_seq
	
	fasta.close()
	
	for contig in fasta_seq_dict:
		print("adapters_dict[\""+contig+"\"] = \""+ fasta_seq_dict[contig]+"\"")
		
#These spit out spoofed python code for the in-built creation of an adapters file to supply tools, without the external need for this file.
#I just copy-paste the results to tbe generate_adapters_temporary_file function's body to get the results.
#OG = "/mnt/c/Users/Kenji/Desktop/NovaSeq_invest/pre_QC/all_adapters.txt"
#fasta_to_permanent_python(OG)


#As above, this prints out a python-correct set of commands for me to copy-paste.
#This one produces a set of "families" for adapters, where each is the set of adapters in a kit.
def fasta_to_families(original_adapters_fasta):
	fasta = open(original_adapters_fasta, "r")

	whichfam = [5, 6, 10, 6, 7, 12, 14, 3, 13, 50, 9, 5, 6]
	
	families = []
	families.append("singleend")
	families.append("pairedend")
	families.append("dpnII")
	families.append("smallrna")
	families.append("multiplex")
	families.append("pcr")
	families.append("dpnIIgex")
	families.append("otherrna")
	families.append("trueseq")
	families.append("rnapcr")
	families.append("abi")
	families.append("trueseq2")
	families.append("nextera")
	
	current_fam = 0
	
	famlist = []
	
	for family_size in whichfam:
		for i in range(0, family_size):
			famlist.append(families[current_fam])
		current_fam += 1
	
	fasta_fam_dict = {}
	
	current_fam = 0
	
	for line in fasta:
		if line.strip().startswith(">"):
			fasta_fam_dict[line.strip()[1:]] = famlist[current_fam]
			current_fam += 1
			
	fasta.close()
	
	for contig in fasta_fam_dict:
		print("adapters_fam_dict[\""+contig+"\"] = \""+ fasta_fam_dict[contig]+"\"")


#OG = "/mnt/c/Users/Kenji/Desktop/NovaSeq_invest/pre_QC/all_adapters.txt"		
#fasta_to_families(OG)
		

