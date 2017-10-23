#     file: solver.py
#   author: Jesse Eaton
#  created: 9/30/2017
# modified: 10/14/2017
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
import gurobipy


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

U_MIN = 1*10**(-5)
MAX_SOLVER_ITERS = 5000


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

#  input: F (np.array of float) [m, l+r] mixed copy number f_p,s of mutation s in sample p
#         Q (np.array of 0 or 1) [l, r] q_b,s == 1 if breakpoint b is in segment s. 0 otherwise
#         G (np.array of 0 or 1) [l, l] g_s,t == 1 if breakpoints s and t are mates. 0 otherwise
#         A (np.array of int) [m, l] a_p,b is number of mated reads for breakpoint b in sample p
#         H (np.array of int) [m, l] h_p,b is number of total reads for breakpoint b in sample p
#         n (int) number of leaves in phylogeny. 2n-1 is total number of nodes
#         c_max (int) maximum allowed copy number for any element in output C
#         lamb1 (float) regularization term to weight total tree cost against unmixing error
#         lamb2 (float) regularization term to weight breakpoint frequency error
#         max_iters (int) maximum number of iterations to predict U then C if convergence not reached
# output: U (np.array of float) [m, 2n-1] 0 <= u_p,k <= 1. percent of sample p made by clone k
#         C (np.array of int) [2n-1, l+r] int copy number c_k,s of mutation s in clone k
#         E (np.array of int) [2n-1, 2n-1] e_i,j == 1 iff edge (i,j) is in tree. 0 otherwise
#         obj_val (float) objective value of final solution
#         err_msg (None or str) None if no error occurs. str with error message if one does
#  notes: l (int) is number of structural variants. r (int) is number of copy number regions
def get_UCE(F, Q, G, A, H, n, c_max, lamb1, lamb2, max_iters):
	np.random.seed() # sets seed for running on multiple processors
	m = len(F)

	for i in xrange(0, max_iters):

		if i == 0:
			U = gen_U(m, n)
		else:
			U = get_U(F, C, n)

		obj_val, C, E, R, err_msg = get_C(F, U, Q, G, A, H, n, c_max, lamb1, lamb2)

		# handle errors
		if err_msg != None:
			return None, None, None, None, err_msg

		# C_diff = abs((C - C_prv)).sum() # could possibly return

	return U, C, E, obj_val, None


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
	prb.solve(solver = cvx.GUROBI)

	# remove any numbers extremely close to zero
	U = np.array(U.value)
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
#         lamb1 (float) regularization term to weight total tree cost against unmixing error
#         lamb2 (float) regularization term to weight breakpoint frequency error
# output: obj_val (float) objective value of solution
#         C (np.array of int) [2n-1, l+r] int copy number c_k,s of mutation s in clone k
#         E (np.array of int) [2n-1, 2n-1] e_i,j == 1 iff edge (i,j) is in tree. 0 otherwise
#         R (np.array of int) [2n-1, 2n-1] cost of each edge in the tree
#         err_msg (None or str) None if no error occurs. str with error message if one does
#  notes: l (int) is number of structural variants. r (int) is number of copy number regions
def get_C(F, U, Q, G, A, H, n, c_max, lamb1, lamb2):
	l, r = Q.shape
	m, _ = U.shape
	C = cvx.Int(2 * n - 1, l + r)

	E = cvx.Int(2 * n - 1, 2 * n - 1) # binary edge indicators
	R = cvx.Int(2 * n - 1, 2 * n - 1) # rho. cost across each edge
	W = _get_bp_appearance_var(n, l)  # w_bin (dict). binary breakpoint appearance indicator
	C_bin = cvx.Int(2*n-1, l+r)       # binary representation of C. C_bin = (C == 1)
	Gam = cvx.Int(2*n-1, l)           # gamma. copy number of segment containing breakpoint

	F_seg = F[:, l:].dot(np.transpose(Q)) # [m, l] mixed copy number of segment containing breakpoint
	Pi = np.divide(F[:, :l], F_seg)       # [m, l] expected bpf (ratio of bp copy num to segment copy num)

	# add constraints
	cst = []
	cst += _get_copy_num_constraints(C, c_max, l, r, n)
	cst += _get_tree_constraints(E, n)
	cst += _get_cost_constraints(R, C, E, n, l, r, c_max)
	cst += _bin(C, C_bin, 2*n-1, l+r, c_max)
	cst += _get_bp_appearance_constraints(C_bin, W, E, G, n, l)
	cst += _get_sum_condition_constraints(C_bin, W, U, m, n, l)
	cst += _get_segment_copy_num_constraints(Gam, C, Q, W, m, n, l, r)
	# cst += _get_bp_frequency_constraints(C, U, Q, A, H, m, n, l, r, alpha)

	coef_1 = 1.0 / float(m * (l + r)) # normalizing for worst case F - U * C
	coef_2 = 1.0 / float(r * (2*n-1)) # normalizing for worst case cost of tree (R)
	coef_3 = 1.0 / float(m * l)       # normalizing for worst case difference in segment estimate using bpf for bp copy num
	S = cvx.abs(cvx.mul_elemwise(Pi, cvx.sum_entries(U * Gam)) - cvx.sum_entries(U * C[:, :l]))
	obj = cvx.Minimize(coef_1 * cvx.sum_entries(cvx.abs(F - U * C)) + lamb1 * coef_2 * cvx.sum_entries(R) + lamb2 * coef_3 * cvx.sum_entries(S))
	prb = cvx.Problem(obj, cst)

	try:
		prb.solve(solver = cvx.GUROBI, max_iters = MAX_SOLVER_ITERS)
	except:
		return None, None, None, None, "ERROR: solving with " + str(cvx.GUROBI) + " did not work"

	if prb.status == cvx.OPTIMAL:
		C = np.array(C.value.round())
		R = np.array(R.value.round())
		E = np.array(E.value.round())
		return prb.value, C, E, R, None

	return prb.value, None, None, None, "error when solving: " + str(prb.status)

# # # # # # # # # # # # # # # # # # # # # # # #
#   C O N S T R A I N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # #

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
			cst.append(E[i, j] + E[j, i] <= 1)

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

#  input: C_bin (cvx.Int) [2n-1, l+r] binary copy number c_k,s of mutation s in clone k. is 0 if cp# is 0. 1 otherwise
#         W (dict of cvx.Int) key is bp_index (int). val is w_i,j (cvx.Int) bp appearance indicator variables
#                             this should have no constraints yet
#         E (cvx.Int) [2n-1, 2n-1] binary edge indicator. 1 if (i,j) is an edge. 0 otherwise
#         G (np.array of int) g_s,t == 1 iff breakpoints s and t are mates. 0 otherwise
#         n (int) number of leaves in tree. total number of nodes is 2n-1
#         l (int) number of breakpoints
# output: cst (list of cvx.Constraint) W_i,j,b == 1 iff bp b appears on edge (i,j). for each b, W_i,j,b == 1 once
def _get_bp_appearance_constraints(C_bin, W, E, G, n, l):
	N = 2*n-1

	cst = []                # make each w_i,j,b binary
	for b in xrange(0, l):
		cst.append(0 <= W[b])
		cst.append(W[b] <= 1)

	X = {}                  # x_i,j,b == 0 iff binary cp# goes from 0 to 1 across edge (i,j). X_i,j,b > 0 otherwise
	for b in xrange(0, l):
		X[b] = cvx.Int(N, N)
		for i in xrange(0, N):
			for j in xrange(0, N):
				cst.append(X[b][i, j] == 2 + C_bin[i, b] - C_bin[j, b] - E[i, j])

	X_bin = {}              # x_bin_i,j,b == 0 iff binary cp# goes from 0 to 1 across edge (i,j). 1 otherwise
	for b in xrange(0, l):
		X_bin[b] = cvx.Int(N, N)
		cst += _bin(X[b], X_bin[b], N, N, 3)   # largest x_i,j,b could be is 3

		cst.append(W[b] == 1 - X_bin[b])       # w_i,j,b == 1 iff bp b appears on edge (i,j). 0 otherwise

	for b in xrange(0, l):                   # breakpoints only appear once in the tree
		cst.append(cvx.sum_entries(W[b]) == 1)

	for s in xrange(0, l):                   # breakpoint pairs must appear on same edge
		for t in xrange(0, l):
			cst.append(W[s] - W[t] <= 1 - G[s, t])
			cst.append(W[t] - W[s] <= 1 - G[s, t])
	
	return cst

#  input: C_bin (cvx.Int) [2n-1, l+r] binary copy number c_k,s of mutation s in clone k. is 0 if cp# is 0. 1 otherwise
#         W (dict of cvx.Int) key is bp_index (int). val is w_i,j (cvx.Int) == 1 iff bp b appears on edge (i,j)
#         U (np.array of float) [m, 2n-1] 0 <= u_p,k <= 1. percent of sample p made by clone k
#         m (int) number of samples
#         n (int) number of leaves in tree. total number of nodes is 2n-1
#         l (int) number of breakpoints
# output: cst (list of cvx.Constraint)
#   does: demands breakpoints occuring higher in the tree have larger percent of cells with this bp must he larger
def _get_sum_condition_constraints(C_bin, W, U, m, n, l):
	N = 2*n-1
	cst = []

	W_node = {}              # w_node_j,b == 1 iff bp b appears at node j. 0 otherwise
	for b in xrange(0, l):
		W_node[b] = cvx.Int(N)
		for j in xrange(0, N):
			cst.append(W_node[b][j] == sum([ W[b][i, j] for i in xrange(0, N) ]))

	Phi = cvx.Variable(m, l)          # phi_p,b is percent of cells in sample p with breakpoint b
	cst.append(Phi == U * C_bin[:, :l])

	Y = {}                        # y_i,j,s,t == 0 iff breakpoint s occurs at node i then t immediately occurs
	for i in xrange(0, N):        #   at node j. > 0 otherwise
		for j in xrange(0, N):
			Y[(i, j)] = cvx.Int(l, l)
			for s in xrange(0, l):
				for t in xrange(0, l):
					cst.append(Y[i, j][s, t] == 2 - W_node[s][i] - W[t][i, j])

	Y_bin = {}                    # y_bin_i,j,s,t is binary version of y. so now any > 0 is 1
	for i in xrange(0, N):
		for j in xrange(0, N):
			Y_bin[(i, j)] = cvx.Int(l, l)
			cst += _bin(Y[(i, j)], Y_bin[(i, j)], l, l, 2)

	D = cvx.Int(l, l)             # variant parent indicator d_s,t == 1 iff breakpoint s appears in parent
	for s in xrange(0, l):        #   of node where breakpoint t appears
		for t in xrange(0, l):
			cst.append(D[s, t] == sum([ (1 - Y_bin[(i, j)][s, t]) for i in xrange(0, N) for j in xrange(0, N) ]))

	for p in xrange(0, m):        # cell fraction of breakpoint in parent must be > cell frac in child
		for s in xrange(0, l):
			for t in xrange(0, l):
				cst.append(Phi[p, s] >= Phi[p, t] - 1 + D[s, t])

	return cst

#  input: Gam (cvx.Int) [2n-1, l] copy number g_k,b if segment containing breakpoint b in clone k
#                                 should have no constraints added so far!
#         C (cvx.Int) [2n-1, l+r] copy number c_k,s of mutation s in clone k
#         Q (np.array of 0 or 1) [l, r] q_b,s == 1 if breakpoint b is in segment s. 0 otherwise
#         W (dict of cvx.Int) key is bp_index (int). val is w_i,j (cvx.Int) == 1 iff bp b appears on edge (i,j)
# output: cst (list of cvx.Constraint)
#   does: constrains Gamma (segment copy number) to be at least 1 when bp appears and always greater than bp cp num
def _get_segment_copy_num_constraints(Gam, C, Q, W, m, n, l, r):
	N = 2*n-1
	cst = []

	for k in xrange(0, N):
		for b in xrange(0, l):
			cst.append(Gam[k, b] == sum([ Q[b, s] * C[k, l+s] for s in xrange(0, r) ]))

	cst.append(C[:, :l] <= Gam) # copy number of breakpoint cannot exceed copy number of segment containing bp

	for b in xrange(0, l):   # copy number of segment containing bp must be at least 1 if bp appears at node j
		for j in xrange(0, N):
			cst.append(Gam[j, b] >= sum([ W[b][i, j] for i in xrange(0, N) ]))

	return cst

#  input: C (cvx.Int) [2n-1, l+r] copy number c_k,s of mutation s in clone k
#         U (np.array of float) [m, 2n-1] 0 <= u_p,k <= 1. percent of sample p made by clone k
#         Q (np.array of 0 or 1) [l, r] q_b,s == 1 if breakpoint b is in segment s. 0 otherwise
#         A (np.array of int) [m, l] a_p,b is number of mated reads for breakpoint b in sample p
#         H (np.array of int) [m, l] h_p,b is number of total reads for breakpoint b in sample p
#         m (int) number of samples
#         n (int) number of leaves in tree. total number of nodes is 2n-1
#         l (int) number of breakpoints
#         r (int) number of copy number segments
#         alpha (float) number of standard deviations allowed for estimator of breakpoint frequency
# output: cst (list of cvx.Constraint)
#   does: forces the ratio copy number of breakpoint to copy number of segment containing breakpoint
#           (reguardless of whether it is mated or not) to be approximately the breakpoint frequency
def _get_bp_frequency_constraints(C, U, Q, A, H, m, n, l, r, alpha):
	N = 2*n-1
	cst = []

	Gam = cvx.Int(N, l)      # gam_k,b is segment copy number of bp b at clone k
	for k in xrange(0, N):
		for b in xrange(0, l):
			cst.append(Gam[k, b] == sum([ Q[b, s] * C[k, l+s] for s in xrange(0, r) ]))

	Pi = np.divide(A, H) # breakpoint frequency (bpf)
	Std_pi = np.sqrt(np.divide(np.multiply(Pi, 1-Pi), H)) # standard deviation in bpf. assume binomial
	
	F_bp, F_sg = cvx.Variable(m, l), cvx.Variable(m, l)
	cst.append(F_bp == U * Gam)      # f_bp_p,b is mixed copy number of breakpoint (with mate) b in sample p
	cst.append(F_sg == U * C[:, :l]) # f_sg_p,b is mixed copy number of segment containing bp b in sample p

	cst.append( cvx.mul_elemwise(Pi - alpha*Std_pi, F_bp) <= F_sg ) # ratio of f_bp / f_sg should be close
	cst.append( cvx.mul_elemwise(Pi + alpha*Std_pi, F_bp) >= F_sg ) #   to bpf Pi

	return cst


# # # # # # # # # # # # # # # # # # # #
#   H E L P E R   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # #

#  input: X (cvx.Int) [m, n] matrix to binarize
#         Y (cvx.Int) [m, n] matrix to add constraints to. this will be the binary matrix version of X
#         m, n (int) number of rows and cols respectively for both inputs X and Y
#         x_max (int) maximum allowed value in input X. this is used to create bit representation of X
# output: cst (list of cvx.Constraint) constraints making Y vals = 1 iff X > 0 and Y vals == 0 otherwise
def _bin(X, Y, m, n, x_max):
	num_bits = int(math.floor(math.log(x_max, 2))) + 1
	cst = [Y <= 1]

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
			for b in xrange(0, num_bits):
				cst.append( Z[b][i, j] <= Y[i, j] )

			# constrain Y to be 0 if all bits are 0
			cst.append( Y[i, j] <= sum([ Z[b][i, j] for b in xrange(0, num_bits) ]) )

	return cst

#  input: n (int) number of leaves in tree. total number of nodes is 2n-1
#         l (int) number of breakpoints
# output: W (dict of cvx.Int) key is bp_index (int). val is w_i,j (cvx.Int) bp appearance indicator variables
def _get_bp_appearance_var(n, l):
	N = 2*n-1
	W = {}
	for b in xrange(0, l):
		W[b] = cvx.Int(N, N)
	return W

# generate random U matrix with m rows and 2n-1 cols. vals are between 0.0 and 1.0 and rows sum to 1.0
def gen_U(m, n):
	U = np.random.rand(m, 2*n-1)
	rowsums = np.sum(U, 1)
	for i in xrange(m):
		U[i, :] = U[i, :] / rowsums[i]
	return U

def printnow(s):
	sys.stdout.write(s)
	sys.stdout.flush()
