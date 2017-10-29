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
from graphviz import Digraph

# custom modules
sys.path.insert(0, 'model/')
sys.path.insert(0, 'help/')
import solver as sv
import file_manager as fm      # sanitizes file and directory arguments
import generate_matrices as gm # gets F, Q, G, A, H from .vcf files
import printer as pt


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

MAX_NUM_LEAVES = 10
MAX_COPY_NUM = 20
MAX_CORD_DESC_ITERS = 1000
MAX_RESTART_ITERS = 1000
NUM_CORES = mp.cpu_count()


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	
	F, Q, G, A, H = gm.get_mats(args['input_directory'])
	l, r = Q.shape
	m = len(F)

	if args['lambda1'] == None:
		# args['lambda1'] = float(r) / float((l + r) * m) # r/(l+r) * 1/m
		args['lambda1'] = 0.05
	if args['lambda2'] == None:
		args['lambda2'] = 0.05
		# args['lambda2'] = float(l) / float(l + r)

	n, c_max, lamb1, lamb2 = args['num_leaves'], args['c_max'], args['lambda1'], args['lambda2']
	num_restarts, num_cd_iters, num_processors = args['restart_iters'], args['cord_desc_iters'], args['processors']

	check_valid_input(Q, G, A, H)

	p = mp.Pool(processes = num_processors)

	arg_set = (F, Q, G, A, H, n, c_max, lamb1, lamb2, num_cd_iters)
	arg_sets_to_use = [ arg_set for _ in xrange(0, num_restarts) ]

	Us, Cs, Es, obj_vals = [], [], [], []
	Rs, Ws = [], []
	num_complete = 0
	while arg_sets_to_use:
		arg_sets = arg_sets_to_use
		arg_sets_to_use = []

		for res in p.imap_unordered(setup_get_UCE, arg_sets):
			U, C, E, R, W, obj_val, err_msg = res
			if err_msg != None:
				# arg_sets_to_use.append(arg_set)
				# printnow('failure... putting task back on the queue\n')
				printnow('failure... task terminated\n')
			else:
				num_complete += 1
				printnow(str(num_complete) + ' of ' + str(num_restarts) + ' complete\n')
				Us.append(U)
				Cs.append(C)
				Es.append(E)
				Rs.append(R)
				Ws.append(W)
				obj_vals.append(obj_val)

	best_i = 0
	best_obj_val = obj_vals[best_i]
	for i, obj_val in enumerate(obj_vals):
		if obj_val < best_obj_val:
			best_obj_val = obj_val
			best_i = i

	write_to_files(args['output_directory'], Us[best_i], Cs[best_i], Es[best_i], Rs[best_i], Ws[best_i], F, obj_vals[best_i])

def setup_get_UCE(args):
	return sv.get_UCE(*args)

def printnow(s):
	sys.stdout.write(s)
	sys.stdout.flush()

# d (str) is local directory path. all others are np.array
def write_to_files(d, U, C, E, R, W, F, obj_val):
	fnames = [ d + fname for fname in ['U.tsv', 'C.tsv', 'T.dot', 'F.tsv', 'obj_val.txt'] ]
	for fname in fnames:
		fm.touch(fname)
	np.savetxt(fnames[0], U, delimiter = '\t', fmt = '%.8f')
	np.savetxt(fnames[1], C, delimiter = '\t', fmt = '%i')
	np.savetxt(fnames[3], F, delimiter = '\t', fmt = '%.8f')
	np.savetxt(fnames[4], np.array([obj_val]), delimiter = '\t', fmt = '%.8f')
	dot = to_dot(E, R, W)
	open(fnames[2], 'w').write(dot.source) # write tree T in dot format
	dot.render(d + 'T')                    # display tree T in .png


#  input: E (np.array of int) [2n-1, 2n-1] 0 if no edge, 1 if edge between nodes i and j
#         R (np.array of int) [2n-1, 2n-1] cost of each edge in the tree
#         W (np.array of int) [2n-1, 2n-1] number of breakpoints appearing along each edge in tree
# output: dot (graphviz.dot.Digraph) directed tree representation of E
def to_dot(E, R, W):
	N = len(E)
	dot = Digraph(format = 'png')
	dot.node(str(N-1))
	for i in xrange(N-1, -1, -1):
		for j in xrange(N-1, -1, -1):
			if int(E[i, j]) == 1:
				edge_label = ' ' + str(int(R[i, j])) + '/' + str(int(W[i, j]) / 2)
				dot.node(str(j))
				dot.edge(str(i), str(j), label = edge_label)
	return dot

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

# input: Q (np.array of 0 or 1) [l, r] q_b,s == 1 if breakpoint b is in segment s. 0 otherwise
#        G (np.array of 0 or 1) [l, l] g_s,t == 1 if breakpoints s and t are mates. 0 otherwise
#        A (np.array of int) [m, l] a_p,b is number of mated reads for breakpoint b in sample p
#        H (np.array of int) [m, l] h_p,b is number of total reads for breakpoint b in sample p
#  does: exits with error message if any of the input is not valid
def check_valid_input(Q, G, A, H):
	l, r = np.shape(Q)
	m = np.shape(A)[0]
	Q_msg = 'There is an issue with input binary matrix Q (indicates which segment each breakpoint belongs to). Each breakpoint must belong to exactly one segment.'
	G_msg = 'There is an issue with input binary matrix G (indicates which breakpoints are mates). Each breakpoint must be mated into pairs.'
	A_msg = 'There is an issue with input integer matricies A and H (indicating the number of reads mapped to each mated breakpoint and the number of total reads mapping to a breakpoint). The number of mated reads must be less or equal to the total reads and both should be non negative.'

	raiseif(not np.all(np.sum(Q, 1) == 1), Q_msg)

	raiseif(not np.all(np.sum(G, 0) == 2) or not np.all(np.sum(G, 0) == 2), G_msg)
	for i in xrange(0, l):
		for j in xrange(0, l):
			raiseif(G[i, j] != G[j, i], G_msg)
			raiseif(i == j and G[i, i] != 1, G_msg)

	for p in xrange(0, m):
		for b in xrange(0, l):
			raiseif(A[p, b] < 0 or A[p, b] > H[p, b], A_msg)

# raises exception if boolean is true
def raiseif(should_raise, msg):
	if should_raise:
		raise Exception(msg)

	# condition for G and A and H

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'template.py', description = "unmixes mixed copy numbers for breakpoints and segments and infers phylogeny with various phylogenetic constraints")
	parser.add_argument('-i', '--input_directory', required = True, type = lambda x: fm.valid_dir_ext(parser, x, '.vcf'), help = 'directory containing a .vcf for each sample from a single patient')
	parser.add_argument('-o', '--output_directory', required = True, type = lambda x: fm.valid_dir(parser, x), help = 'empty directory for output U.tsv, C.tsv, and T.dot files to go')
	parser.add_argument('-n', '--num_leaves', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 2, MAX_NUM_LEAVES), help = 'number of leaves for inferred binary tree. total number of nodes will be 2*n-1')
	parser.add_argument('-c', '--c_max', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, MAX_COPY_NUM), help = 'maximum allowed copy number at any node in the tree')
	parser.add_argument('-l', '--lambda1', type = lambda x: fm.valid_float_above(parser, x, 0.0), help = 'regularization term to weight total tree cost against unmixing error in objective function. setting as 0.0 will put no tree cost constraint. setting as 1.0 will equally consider tree cost and unmixing error.')
	parser.add_argument('-a', '--lambda2', type = lambda x: fm.valid_float_above(parser, x, 0.0), help = 'regularization term to weight error in inferred ratio between copy number of a breakpoint and the copy number of the segment originally containing the position of breakpoint')
	parser.add_argument('-t', '--cord_desc_iters', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, MAX_CORD_DESC_ITERS), help = 'maximum number of cordinate descent iterations for each initialization of U')
	parser.add_argument('-r', '--restart_iters', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, MAX_RESTART_ITERS), help = 'number of random initializations for picking usage matrix U')
	parser.add_argument('-p', '--processors', required = True, type = lambda x: fm.valid_int_in_range(parser, x, 1, NUM_CORES), help = 'number of processors to use')
	return vars(parser.parse_args(argv))

# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
