# Compare different orthogonal polynomials.

import numpy as np
import numpy.polynomial as nppoly
import matplotlib.pyplot as plt
from numpy import abs, sqrt, exp, log, log10, sin, cos, tan, arccos, arcsin, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh, pi
from matplotlib.pyplot import figure, clf, plot, scatter, xlabel, ylabel, xlim, ylim, legend

import matplotlib as mpl
# see mpl.rcParams for full list
mpl.rc('figure', figsize=(1600/240.0, 1200/240.0), dpi=240)
mpl.rc('lines', markersize=2.0)
mpl.rcParams['axes.formatter.limits'] = [-3,3]
plt.ion()

exec(open('/home/xyy/Documents/research/poly-orth/ex/weighted_orthpoly_solver.py').read())

#f = lambda x: abs(x)
f = lambda x: 1/(1+x**2)
#f = lambda x: sin(10*pi*x+0.1)

f_weighting = lambda x: 1
#f_weighting = lambda x: exp(-x*x)
#f_weighting = lambda x: 1 + cos(pi*x)   # bad sampling at the boundary

n_deg = 100

# Coefficients for normalized Legendre polynomials.
x_sample_leg = nppoly.legendre.legroots(1*(np.arange(n_deg+2)==n_deg+1))
y_sample_leg = f(x_sample_leg)
coef_l = nppoly.legendre.legfit(x_sample_leg, y_sample_leg, n_deg) \
         * sqrt(np.arange(n_deg+1)+0.5)

# Coefficients for weighted orthogonal polynomials.
#x_sample_che = sin(pi * (np.arange(n_deg+1)-(n_deg)/2)/(n_deg+1))
x_sample_che = y_sample_leg
y_sample_che = f(x_sample_che)
basis_repr_coef = get_orthpoly(n_deg, f_weighting, 10)
coef_w = np.linalg.solve(dot(nppoly.chebyshev.chebvander(x_sample_che, n_deg), basis_repr_coef), y_sample_che)

# Check accuracy of the generated polynomials.
n_b = 100
ll_c = nppoly.chebyshev.Chebyshev(basis_repr_coef[:, n_b])
ll   = nppoly.legendre.Legendre.basis(n_b) * sqrt(n_b+0.5)
ssx = np.linspace(-1,1,1000)
figure(323); clf()
plot(ssx, ll_c(ssx), ssx, ll(ssx))
print("diff at 1 = %g" % (ll_c(1) - ll(1)))
print("diff at 1 = %g" % (ll_c(0.1) - ll(0.1)))

coefs = np.vstack([
 nppoly.chebyshev.Chebyshev.interpolate(f, n_deg).coef,
 coef_l,
 coef_w,
 ])

poly_name = ['Cheb', 'Leg', 'weight']

figure(10); clf()
plot(np.arange(n_deg+1), log10(abs(coefs.T)), '-o')
xlabel('k')
ylabel('log10(c_k)')
legend(poly_name)
plt.show()
