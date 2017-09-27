import cvxpy as cvx
import numpy as np

# data
m = 30
n = 20
np.random.seed(1)
A = np.random.randn(m, n) # returns m by n array of values drawn from N(0, 1)
b = np.random.randn(m)    # returns array of len m of values ...

# define problem
x = cvx.Variable(n)
obj = cvx.Minimize(cvx.sum_squares(A * x - b))

# add constraints
constraints = []
constraints.append(0 <= x)
constraints.append(x <= 1)

prob = cvx.Problem(obj, constraints)

print prob

