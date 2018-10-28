#

import numpy as np
import numpy.polynomial as nppoly
import numpy.polynomial.chebyshev  as cheb
import numpy.polynomial.legendre   as lege
from numpy import exp, cos

exec(open('common_imports_plot.py').read())

#from weighted_orthpoly_solver import get_orthpoly
exec(open('weighted_orthpoly_solver.py').read())

n_sampling = 100
x_sampling, int_weight = lege.leggauss(n_sampling)

n_deg = 100
#f_weighting = lambda x: 1
#f_weighting = lambda x: exp(-x*x)
f_weighting = lambda x: 1 + cos(pi*x)
basis_repr_coef = get_orthpoly(n_deg, f_weighting, 10)

ii = 99
jj = 3
ff = nppoly.chebyshev.Chebyshev(basis_repr_coef[:, ii])
gg = nppoly.chebyshev.Chebyshev(basis_repr_coef[:, jj])

ff1 = lege.Legendre.basis(ii) * sqrt(ii + 0.5)
gg1 = lege.Legendre.basis(jj) * sqrt(jj + 0.5)

print('ff (1) = ', ff(1))
print('ff1(1) = ', ff1(1))

s = sum([ff(x)*gg(x)* f_weighting(x) * w for x, w in zip(x_sampling, int_weight)])
print('The <ff1, gg1>rho =', s)

s = sum([ff(x)*gg(x)* w for x, w in zip(x_sampling, int_weight)])
print('The <ff1, gg1>    =', s)

s = sum([ff1(x)*gg1(x)* f_weighting(x) * w for x, w in zip(x_sampling, int_weight)])
print('The <ff, gg>rho   =', s)

# Check the Gramian matrix.
ss = np.zeros((20, 20))
for ii in range(ss.shape[0]):
    for jj in range(ss.shape[1]):
        ff = nppoly.chebyshev.Chebyshev(basis_repr_coef[:, ii])
        gg = nppoly.chebyshev.Chebyshev(basis_repr_coef[:, jj])
        ss[ii,jj] = sum([ff(x)*gg(x)* f_weighting(x) * w for x, w in zip(x_sampling, int_weight)])

figure(14); clf()
matshow(ss - np.eye(ss.shape[0], ss.shape[1]), fignum=False)
title('Gramian matrix: <pi, pj>rho')
colorbar()

