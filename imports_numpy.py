# Import numpy and common aliases.

import numpy as np
import numpy.random as npr
from numpy import dot, linspace, pi, diag, arange, reshape
from numpy import abs, sqrt, exp, log, log10, sin, cos, tan, arccos, arcsin, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh
from numpy.linalg import norm
from pandas import DataFrame

import numpy.polynomial.chebyshev  as cheb
import numpy.polynomial.legendre   as lege
import numpy.polynomial.polynomial as poly

def normalize(v):
    return v / norm(v)

def maxabs(x):
    return np.max(np.abs(x.ravel()))

def colvec(x):
    return reshape(np.array(x), [-1, 1])

def unit_vec(n, j):
    return 1*(arange(n)==j-1)

# vim: et sw=4 sts=4
