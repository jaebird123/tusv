#     file: tusv.py
#   author: Jesse Eaton
#  created: 10/13/2017
# modified: 10/14/2017
#  purpose: Unmixes mixed copy numbers for breakpoints and segments and infers phylogeny
#             with various phylogenetic constraints


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import random
import numpy as np
import multiprocessing as mp

# custom modules
sys.path.insert(0, 'model/')
import solver as sv


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

# here


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	
	F, Q, G, A, H, n, c_max, lamb, alpha = get_vars()

	num_restarts = 4
	max_iters = 2

	p = mp.Pool(processes = 2)

	arg_set = (F, Q, G, A, H, n, c_max, lamb, alpha, max_iters)
	arg_sets_to_use = [ arg_set for _ in xrange(0, num_restarts) ]

	Us, Cs, Es, obj_vals = [], [], [], []
	num_complete = 0
	while arg_sets_to_use:
		arg_sets = arg_sets_to_use
		arg_sets_to_use = []

		for res in p.imap_unordered(setup_get_UCE, arg_sets):
			U, C, E, obj_val, err_msg = res
			if err_msg != None:
				arg_sets_to_use.append(arg_set)
				printnow('failure... putting task back on the queue\n')
			else:
				num_complete += 1
				printnow(str(num_complete) + ' of ' + str(num_restarts) + ' complete\n')
				Us.append(U)
				Cs.append(C)
				Es.append(E)
				obj_vals.append(obj_val)

	print ''
	for i, obj_val in enumerate(obj_vals):
		print obj_val
		print Us[i]
		print Cs[i]
		print Es[i]
		print ''

def setup_get_UCE(args):
	return sv.get_UCE(*args)

def printnow(s):
	sys.stdout.write(s)
	sys.stdout.flush()

# temporary functions. REMOVE LATER!!!

def get_vars():
	m = 1       # samples
	n = 3       # leaves
	l = 6       # breakpoints
	r = 10      # segments
	c_max = 7  # maximum copy number
	f_scale = 5

	F = f_scale * np.random.rand(m, l + r)
	Q = np.array([ np.arange(0, r) == random.randint(0, r-1) for bp in xrange(l) ], dtype = int)
	G = gen_G(l)
	A = np.random.binomial(100, 0.25, [m, l])
	H = 100 * np.ones([m, l])
	lamb = 1.0
	alpha = 1.5

	return F, Q, G, A, H, n, c_max, lamb, alpha, 

def gen_G(l):
	G = np.zeros((l, l))
	I = [ x for x in xrange(0, l) ]   # list of all indicies
	random.shuffle(I)                 # randomly permut to make random pairs
	I = np.array(I).reshape((l/2, 2)) # make a l/2 by 2 numpy array of mated pairs
	for i, j in I:
		G[i, j] = 1
		G[j, i] = 1
		G[i, i] = 1
		G[j, j] = 1
	return G

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'template.py', description = "purpose")
	# parser.add_argument('input_file', help = 'input .txt file', type = lambda x: is_valid_file(parser, x))
	return vars(parser.parse_args(argv))

# def is_valid_file(parser, arg):
# 	if not os.path.exists(arg):
# 		parser.error('The file \"' + str(arg) + '\" could not be found.')
# 	else:
# 		return open(arg, 'r')


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
