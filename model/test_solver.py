import solver as sv
import random
import sys
import numpy as np

# F = f_scale * np.random.rand(m, l + r)
# U = np.random.rand(m, 2*n-1)
# Q = np.array([ np.arange(0, r) == random.randint(0, r) for bp in xrange(l) ])
# A = np.random.binomial(100, 0.25, [m, l])
# H = 100 * np.ones([m, l])

# lamb = 1.0
# alpha = 2.0

# C = sv.get_C(F, U, Q, A, H, n, c_max, lamb, alpha)
# print C

def printnow(s):
	sys.stdout.write(s)
	sys.stdout.flush()

def main(argv):
	random.seed(1)
	np.random.seed(1)
	m = 3  # samples
	n = 4  # leaves
	l = 6 # breakpoints
	r = 10 # segments
	c_max = 7 # maximum copy number
	f_scale = 5

	F = f_scale * np.random.rand(m, l + r)
	U = gen_U(m, n)
	Q = np.array([ np.arange(0, r) == random.randint(0, r) for bp in xrange(l) ], dtype = int)
	G = gen_G(l)
	A = np.random.binomial(100, 0.25, [m, l])
	H = 100 * np.ones([m, l])
	lamb = 1.0
	alpha = 2.0

	test_get_U(F, n, l, r)
	test_get_C(F, Q, G, A, H, n, c_max, lamb, alpha)
	exit()

	num_steps = 50
	for i in xrange(0, num_steps):
		if i > 0:
			prevC = C
			prevU = U
		C = sv.get_C(F, U, Q, A, H, n, c_max, lamb, alpha)
		U = sv.get_U(F, C, n)
		if i > 0:
			diffC = abs((C-prevC)).sum()
			diffU = abs((U-prevU)).sum()
			printnow('\n' + str(i) + '\n')
			printnow('difference in C is ' + str(diffC) + '\n')
			printnow('difference in U is ' + str(diffU) + '\n')
			printnow(str(C) + '\n')

# mated pair binary matrix
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

def gen_C(n, l, r):
	C = np.random.rand(2*n-1, l+r)
	C = (C * 4).round()
	C[2*n-2, :l] = 0
	C[2*n-2, l:] = 2
	return C

# generate random U matrix
def gen_U(m, n):
	U = np.random.rand(m, 2*n-1)
	rowsums = np.sum(U, 1)
	for i in xrange(m):
		U[i, :] = U[i, :] / rowsums[i]
	return U

# # # # # # # # #
#   T E S T S   #
# # # # # # # # #

def test_get_U(F, n, l, r):
	C = gen_C(n, l, r)
	printnow('\ntest_get_U starting\n')
	U = sv.get_U(F, C, n)
	printnow(str(U) + '\n')
	printnow('test_get_U complete\n')

def test_get_C(F, Q, G, A, H, n, c_max, lamb, alpha):
	m = len(F)
	U = gen_U(m, n)
	printnow('\ntest_get_C starting\n')
	C = sv.get_C(F, U, Q, G, A, H, n, c_max, lamb, alpha)
	printnow(str(C) + '\n')
	printnow('test_get_C complete\n')


#
#   CALL TO MAIN
#


if __name__ == "__main__":
	main(sys.argv[1:])

