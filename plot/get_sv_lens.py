#     file: template.py
#   author: authorname
#  created: date
# modified: date
#  purpose: purpose


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import vcf
import numpy as np

sys.path.insert(0, '../help/')
import file_manager as fm   # sanitizes file and directory arguments


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

# here


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	
	sv_lens = []
	for dataset_dir in fm.get_subdir_names(args['directory']):
		dataset_dir = args['directory'] + dataset_dir
		for patient_dir in fm.get_subdir_names(dataset_dir):
			patient_dir = dataset_dir + patient_dir
			for fname in fm.get_fnames_in_dir_with_ext(patient_dir, '.vcf'):
				fname = patient_dir + fname
				id_blacklist = set()
				for rec in get_reader(fname):
					if is_sv_rec(rec):
						my_id = rec.ID
						mt_id = rec.INFO['MATEID'][0]
						if my_id not in id_blacklist and mt_id not in id_blacklist:
							id_blacklist.add(my_id)
							id_blacklist.add(mt_id)
							if on_same_chrm(rec):
								sv_lens.append(get_sv_len(rec))
	print sv_lens

def get_reader(fname):
	return vcf.Reader(open(fname, 'r'))

def is_sv_rec(rec):
	return rec.ID.startswith('sv')

def on_same_chrm(rec):
	my_chrm = rec.CHROM
	mt_chrm = str(rec.ALT[0]).replace('[', '').replace(']', '').split(':')[0]
	return my_chrm == mt_chrm

def get_sv_len(rec):
	my_pos = rec.POS
	mt_pos = int(str(rec.ALT[0]).replace('[', '').replace(']', '').split(':')[1])
	return np.abs(my_pos - mt_pos) + 1


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'template.py', description = "purpose")
	parser.add_argument('-d', '--directory', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'entire simulated data set with results directory removed')
	return vars(parser.parse_args(argv))

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
