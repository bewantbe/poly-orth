# Verify Vandermonde matrix is low conditional number matrix

exec(open('common_imports.py').read())

import numpy.polynomial as nppoly

n_deg = 1000

x_sample = nppoly.legendre.legroots(1*(np.arange(n_deg+2)==n_deg+1))

# Other sample points also works.
#n_sampling = n_deg + 1
#x_sample = sin( pi * (-(n_sampling-1)/2 + np.arange(n_sampling))/n_sampling )
#x_sample = nppoly.chebyshev.chebroots(1*(np.arange(n_deg+2)==n_deg+1))

V = nppoly.legendre.legvander(x_sample, n_deg)
printf('Conditional number: %g (should be a small number)' % np.linalg.cond(V))  # A low number
