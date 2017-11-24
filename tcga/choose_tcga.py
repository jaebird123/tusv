#     file: choose_tcga.py
#   author: Jesse Eaton
#  created: 11/5/2017
# modified: 11/5/2017
#  purpose: to select which TCGA samples to run tusv.py on


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import vcf      # for reading and writing VCF files
import numpy as np
import pandas as pd

sys.path.insert(0, '../help/')
import file_manager as fm   # sanitizes file and directory arguments

# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

CLINICAL_FEATURES_FNAME = 'clinical_features.tsv'


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	# print args['tcga_patients_directory']
	# exit()

	mut_dic = get_mut_dic(args['tcga_patients_directory'])
	met_dic = get_ismet_dic(args['clinical_features_file'])

	mut_dic, met_dic = keep_mutal_keys(mut_dic, met_dic)

	no_met_mut_dic, met_mut_dic = split_mut_by_met(mut_dic, met_dic)

	print_mut_dic(no_met_mut_dic, args['tcga_patients_directory'], args['nmt_dir'])
	print '\n - - - - - - METASTATIC - - - - -\n'
	print_mut_dic(met_mut_dic, args['tcga_patients_directory'], args['met_dir'])

def print_mut_dic(mut_dic, from_master_dir = None, to_dir = None):
	mut_dic = inv_dic(mut_dic)
	i = 0
	for (num_bps, num_cnvs) in sorted(mut_dic.keys(), key = lambda x: x[1]): # sort by num_cnvs
		i += 1
		print_sample(mut_dic[(num_bps, num_cnvs)], num_bps, num_cnvs)
		if from_master_dir != None and to_dir != None:
			from_sample_dir = from_master_dir + mut_dic[(num_bps, num_cnvs)] + '/'
			to_sample_dir = to_dir + str(i).zfill(2) + '_' + mut_dic[(num_bps, num_cnvs)] + '/'
			fm.make_dir(to_sample_dir)
			# fm.cp_dir_contents(from_sample_dir, to_sample_dir)

def inv_dic(dic):
	out = {}
	for key, val in dic.iteritems():
		out[val] = key
	return out

def split_mut_by_met(mut_dic, met_dic):
	no_mets, mets = {}, {}
	for barcode, is_met in met_dic.iteritems():
		if int(is_met) == 0:
			no_mets[barcode] = mut_dic[barcode]
		else:
			mets[barcode] = mut_dic[barcode]
	return no_mets, mets

def keep_mutal_keys(dic1, dic2):
	keys = set(dic1.keys()) & set(dic2.keys())
	out1, out2 = {}, {}
	for key in list(keys):
		if key in dic1 and key in dic2:
			out1[key] = dic1[key]
			out2[key] = dic2[key]
	return out1, out2

def print_sample(barcode, num_bps, num_cnvs):
	print barcode + ':\tnum_bps: ' + str(num_bps) + ',\tnum_cnvs: ' + str(num_cnvs) + '\ttotal: ' + str(num_bps + num_cnvs)


#
#   INPUT
#

def get_mut_dic(d):
	file_dic = get_vcf_file_dic(d)
	dic = {}
	for subdir, reader in file_dic.iteritems():
		barcode = subdir[:-1] # remove '/' from end of subdir string
		num_cnvs, num_bps = 0, 0
		for rec in reader:
			if is_bp(rec):
				num_bps += 1
			if is_cnv(rec):
				num_cnvs += 1
		dic[barcode] = (num_bps, num_cnvs)
	return dic

def is_bp(rec):
	return type(rec.ALT[0]) == vcf.model._Breakend

def is_cnv(rec):
	return type(rec.ALT[0]) == vcf.model._SV

def get_vcf_file_dic(d):
	subdirs = fm.get_subdir_names(d)
	dic = {}
	for subdir in subdirs:
		dic[subdir] = get_reader(d + subdir + fm.get_fnames_in_dir(d + subdir)[0])
	return dic

def get_reader(fname):
	return vcf.Reader(open(fname, 'r'))

def get_file(d, ext):
	fnames = [ d + fname for fname in fnames ]
	for fname in fnames:
		if fname.endswith(ext):
			return open(fname, 'r')

def get_ismet_dic(file):
	df = removeSparseCol(pd.read_csv(file, sep = '\t', index_col = 0))
	barcodes = list(df.index.values)
	is_metas = list(df.as_matrix(['isMetaRecurAny']).flatten())
	dic = {}
	for i, barcode in enumerate(barcodes):
		dic[barcode] = is_metas[i]
	return dic

## also converts all values to float
def removeSparseCol(myDF):
    DF = myDF.copy()
    DF = DF.astype(float)
    sparseCol = []
    for colName in DF:
        if DF.describe().loc['count'][colName] < DF.shape[0]/2:
            sparseCol.append(colName)
    for spCol in sparseCol:
        del DF[spCol]
    return DF

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'choose_tcga.py', description = "to select which TCGA samples to run tusv.py on")
	parser.add_argument('-t', '--tcga_patients_directory', required = True, type = lambda x: fm.valid_master_dir_with_files_and_ext(parser, x, [], '.vcf'), help = 'master directory containing subdirectories. each subdirectory should contain a single .vcf file')
	parser.add_argument('-c', '--clinical_features_file', default = CLINICAL_FEATURES_FNAME, type = lambda x: is_valid_file(parser, x), help = 'file containing clinical information on each of the samples from the tcga patients directory')
	parser.add_argument('-m', '--met_dir', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'file where all TCGA metastatic TCGA directories will be copied to')
	parser.add_argument('-n', '--nmt_dir', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'file where all TCGA non metastatic TCGA directories will be copied to')
	return vars(parser.parse_args(argv))

	# valid_master_dir_with_files_and_ext(parser, arg, fnames, ext)

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return open(arg, 'r')


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
