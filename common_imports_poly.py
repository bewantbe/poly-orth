#

import sys
sys.path.append('/home/xyy/code/py/dnn_hess/')      # for def_dnn

import numpy as np
import numpy.random as npr
from numpy import dot, linspace, pi, diag, arange
from numpy import abs, sqrt, exp, log, log10, sin, cos, tan, arccos, arcsin, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh
from numpy.linalg import norm
from pandas import DataFrame

import numpy.polynomial.chebyshev  as cheb
import numpy.polynomial.legendre   as lege
import numpy.polynomial.polynomial as poly

exec(open('common_imports_plot.py').read())

def normalize(v):
    return v / norm(v)

def maxabs(x):
    return np.max(np.abs(x.ravel()))

# vim: et sw=4 sts=4
