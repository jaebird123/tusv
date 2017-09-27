import cvxpy as cvx
import numpy as np

m = 3  # samples
n = 4  # leaves
l = 10 # breakpoints
r = 20 # segments
c_max = 10 # maximum copy number

# inputs
np.random.seed(1)
F = np.random.randn(m, l+r) ** 2 # returns m by n array of values drawn from N(0, 1)

# variables
U = cvx.Variable(m, 2*n-1)
C = cvx.Int(2*n-1, l+r)

cs = [] # constraints

# U constraints
cs.append(U >= 0)
cs.append(U <= 1)

# C constraints
cs.append(C >= 0)
cs.append(C <= c_max)

# objective
F_hat = U * C
obj = cvx.Minimize(cvx.sum_entries(F - F_hat))

# solve problem
prob = cvx.Problem(obj, cs)
sol = prob.solve

print sol
print 'optimal U'
print U.value
print ''
print 'optimal C'
print C.value
