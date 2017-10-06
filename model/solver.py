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
import cvxpy as cvx
import numpy as np


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

U_MIN = 1*10**(-5)


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

	# E = cvx.Int(2 * n - 1, 2 * n - 1)

	# add constraints
	cst = []
	cst += _get_copy_num_constraints(C, c_max, l, r, n)
	# cst += _get_tree_constraints(E, n)

	obj = cvx.Minimize(cvx.sum_entries(cvx.abs(F - U * C)))
	prb = cvx.Problem(obj, cst)

	prb.solve()

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
	
	# root
	for i in xrange(0, 2*n-1):
		cst.append(E[i, 2*n-2] == 0)
	
	# leaves
	for i in xrange(0, n):
		for j in xrange(0, 2*n-1):
			cst.append(E[i, j] == 0)

	# non-root incoming edges
	E_ir = E[:2*n-2, :] # E incoming rootless
	for j in xrange(0, 2*n-1):
		cst.append(cvx.sum_entries(E_ir[:, j]) == 1)

	# internal node outgoing edges
	for i in xrange(n, 2*n-1):
		cst.append(cvx.sum_entries(E[i, :]) == 2)
	
	return cst
