#     file: solver.py
#   author: Jesse Eaton
#  created: 9/30/2017
# modified: 9/30/2017
#  purpose: Linear program solver of single instance of mixed copy number F, solving for either
#              copy number C or mixture U where the other is assumed constant


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import math     # it's math. we're gonna need it
import cvxpy as cvx
import numpy as np


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

U_MIN = 1*10**(-5)
MAX_PRB_ITERS = 10000


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

#  input: F (np.array of float) [m, l+r] mixed copy number f_p,s of mutation s in sample p
#         C (np.array of int) [2n-1, l+r] int copy number c_k,s of mutation s in clone k
#         n (int) number of leaves in phylogeny. 2n-1 is total number of nodes
# output: U (np.array of float) [m, 2n-1] 0 <= u_p,k <= 1. percent of sample p made by clone k
def get_U(F, C, n):
	m = len(F)
	U = cvx.Variable(m, 2*n-1)
	
	cst = [0 <= U, U <= 1]
	for i in xrange(m):
		U_row = U[i, :]
		cst.append(cvx.sum_entries(U_row) == 1)

	obj = cvx.Minimize(cvx.sum_entries(cvx.abs(F - U * C)))
	prb = cvx.Problem(obj, cst)
	prb.solve()

	# remove any numbers extremely close to zero
	U = U.value
	for i in xrange(m):
		for j in xrange(2*n-1):
			if U[i, j] <= U_MIN:
				U[i, j] = 0.0

	# renormalize U so all rows sum to 1
	rowsums = np.sum(U, 1)
	for i in xrange(m):
		U[i, :] = U[i, :] / rowsums[i]

	return U

#  input: F (np.array of float) [m, l+r] mixed copy number f_p,s of mutation s in sample p
#         U (np.array of float) [m, 2n-1] 0 <= u_p,k <= 1. percent of sample p made by clone k
#         Q (np.array of 0 or 1) [l, r] q_b,s == 1 if breakpoint b is in segment s. 0 otherwise
#         G (np.array of 0 or 1) [l, l] g_s,t == 1 if breakpoints s and t are mates. 0 otherwise
#         A (np.array of int) [m, l] a_p,b is number of mated reads for breakpoint b in sample p
#         H (np.array of int) [m, l] h_p,b is number of total reads for breakpoint b in sample p
#         n (int) number of leaves in phylogeny. 2n-1 is total number of nodes
#         c_max (int) maximum allowed copy number for any element in output C
#         lamb (float) regularization term to weight total tree cost against unmixing error
#         alpha (float) number of standard deviations allowed for estimator of breakpoint frequency
# output: C (np.array of int) [2n-1, l+r] int copy number c_k,s of mutation s in clone k
#  notes: l (int) is number of structural variants. r (int) is number of copy number regions
def get_C(F, U, Q, G, A, H, n, c_max, lamb, alpha):
	l, r = Q.shape
	C = cvx.Int(2 * n - 1, l + r)

	E = cvx.Int(2 * n - 1, 2 * n - 1) # binary edge indicators
	R = cvx.Int(2 * n - 1, 2 * n - 1) # rho. cost across each edge

	# add constraints
	cst = []
	cst += _get_copy_num_constraints(C, c_max, l, r, n)
	cst += _get_tree_constraints(E, n)
	cst += _get_cost_constraints(R, C, E, n, l, r, c_max)

	# REMOVEE COMMENT LATER!!!!
	Y = cvx.Int(2*n-1, l+r)
	cst += _bin(C, Y, 2*n-1, l+r, c_max)

	obj = cvx.Minimize(cvx.sum_entries(cvx.abs(F - U * C)) + lamb * cvx.sum_entries(R))
	# obj = cvx.Minimize(cvx.sum_entries(cvx.abs(F - U * C)))
	prb = cvx.Problem(obj, cst)

	prb.solve(max_iters = MAX_PRB_ITERS)
	print prb.value

	print E.value.round()
	return C.value.round()

def _get_copy_num_constraints(C, c_max, l, r, n):
	cst = [0 <= C, C <= c_max]

	# copy number of breakpoints at root are all zero
	for b in xrange(0, l):
		cst.append(C[2*n-2, b] == 0)

	# copy number of segments at root are all 2
	for s in xrange(l, l+r):
		cst.append(C[2*n-2, s] == 2)

	return cst

def _get_tree_constraints(E, n):
	cst = [0 <= E, E <= 1]

	# no outgoing edges for leaves
	cst.append(E[:n, :] == 0)

	# now only build constraints for internal nodes

	# no edges from descendents to root
	cst.append(E[n : 2*n-1, 2*n-2] == 0)

	# no edges to self (only needed for non leaf, non root nodes b/c no constraints yet)
	for i in xrange(n, 2*n-2):
		cst.append(E[i, i] == 0)

	# internal nodes have 2 outgoing edges
	for i in xrange(n, 2*n-1):
		cst.append(cvx.sum_entries(E[i, :]) == 2)

	# only one edge from ancestors allowed (for every node but root)
	for j in xrange(0, 2*n-2):
		cst.append(cvx.sum_entries(E[n:, j]) == 1)

	# no 2 node cycles
	for i in xrange(n, 2*n-1):
		for j in xrange(n, 2*n-1):
			cst.append(E[i, j] + E[j, i] < 2)

	return cst

def _get_cost_constraints(R, C, E, n, l, r, c_max):
	N = 2*n-1
	X = {}
	for s in xrange(0, r): # x_i,j,s is absolute difference for seg s from node i to j if edge (i,j) exists
		X[s] = cvx.Int(N, N)
	
	cst = []
	for s, _ in X.iteritems():
		cst.append(0 <= X[s])                 # all x_ijs must be >= 0
		for i in xrange(0, N):
			for j in xrange(0, N):
				cst.append(X[s][i, j] <= c_max * E[i, j]) # x_ijs for all segments set to zero if edge (i, j) doesnt exist
				cst.append(X[s][i, j] >= C[i, s+l] - C[j, s+l] - (c_max+1) * (1-E[i, j])) # abs val if edge exists
				cst.append(X[s][i, j] >= C[j, s+l] - C[i, s+l] - (c_max+1) * (1-E[i, j]))

	# define R as cost of transforming copy number profile (cnp) from node i to cnp fr j
	for i in xrange(0, N):
		for j in xrange(0, N):
			cst.append(R[i, j] >= sum([ X[s][i, j] for s in xrange(0, r) ]))

	return cst

#  input: X (cvx.Int) [m, n] matrix to binarize
#         Y (cvx.Int) [m, n] matrix to add constraints to. this will be the binary matrix version of X
#         m, n (int) number of rows and cols respectively for both inputs X and Y
#         x_max (int) maximum allowed value in input X. this is used to create bit representation of X
# output: cst (list of cvx.Constraint) constraints making Y vals = 1 iff X > 0 and Y vals == 0 otherwise
def _bin(X, Y, m, n, x_max):
	num_bits = int(math.floor(math.log(x_max, 2))) + 2
	cst = [0 <= Y, Y <= 1]

	Z = {}
	for b in xrange(0, num_bits): # create binary variables Z
		Z[b] = cvx.Int(m, n)
		cst.append(0 <= Z[b])
		cst.append(Z[b] <= 1)
	for i in xrange(0, m):
		for j in xrange(0, n):

			# set Z as binary representation of X
			cst.append( sum([ Z[b][i, j] * 2**b for b in xrange(0, num_bits) ]) == X[i, j] )

			# constrain Y to be 1 if any bits are 1
			# for b in xrange(0, num_bits):
			# 	cst.append( Z[b][i, j] <= Y[i, j] )

			# # constrain Y to be 0 if all bits are 0
			# cst.append( Y[i, j] <= sum([ Z[b][i, j] for b in xrange(0, num_bits) ]) )

	return cst
