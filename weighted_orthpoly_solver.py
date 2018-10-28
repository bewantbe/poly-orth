# Utilities for generating orthogonal polynomials.

import numpy as np
from numpy import pi, sqrt, sin, dot, sum, sign
from numpy.linalg import inv, qr, cholesky

import numpy.polynomial.chebyshev  as cheb
import numpy.polynomial.legendre   as lege

def get_orthpoly(n_deg, f_weighting, n_extra_point = 10):
    """
     Get orthogonal polynomials with respect to the weighting (f_weighting).
     The polynomials are represented by coefficients of Chebyshev polynomials.
     Evaluate it as: np.dot(cheb.chebvander(set_of_x, n_deg), repr_coef)
     See weighted_orthpoly1.py for validation of this algorithm.
    """
    n_sampling = n_deg + 1 + n_extra_point
    # Do the intergration by discrete sampling at good points.
    if 0:
        # Here we use (first kind) Chebyshev points.
        x_sampling = sin( pi * (-(n_sampling-1)/2 + np.arange(n_sampling))/n_sampling )
        # Put together the weighting and the weighting of the Gauss-Chebyshev quadrature.
        diag_weight = np.array([f_weighting(x)/cheb.chebweight(x) for x in x_sampling]) * (pi/n_sampling)
        V = cheb.chebvander(x_sampling, n_deg)
    else:
        # Use Gauss-Legendre quadrature.
        x_sampling, int_weight = lege.leggauss(n_sampling)
        diag_weight = np.array([f_weighting(x)*w for x, w in zip(x_sampling, int_weight)])
        V = cheb.chebvander(x_sampling, n_deg)
#        V = lege.legvander(x_sampling, n_deg)
    
    # The Gramian matrix.
#    inner_prod_matrix = dot(V.T, diag_weight[:, np.newaxis] * V)
#    repr_coef = inv(cholesky(inner_prod_matrix).T)
    # QR decomposition should give more accurate result.
    repr_coef = inv(qr(sqrt(diag_weight[:, np.newaxis]) * V, mode='r'))
    repr_coef = repr_coef * sign(sum(repr_coef,axis=0))
    return repr_coef
    
# vim: et sw=4 sts=4
