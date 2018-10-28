# Get orthogonal polynomials through eigen factorization.
# Using representation of powers of x, i.e. x^k.
# Demo that the problem is not scalable -- due to the Vandermonde matrix had huge condition number.

exec(open('common_imports_poly.py').read())

# eigen value for 
n_deg = 17
s_x = linspace(-1, 1, n_deg+1)  # For demo
#s_x = np.cos(linspace(0, pi, n_deg+1))  # Cheb 2-kind, helps generalization, but not the problem conditon number

B = poly.polyvander(s_x, n_deg)

figure(10); clf()
matshow(dot(B.T, B), fignum=False)
title('The inner product matrix, ie square of Vandermonde matrix.')
colorbar()

eigval, eigvec = np.linalg.eigh(dot(B.T, B))

print("Condition number %g."%(eigval[-1]/eigval[0]))

figure(12); clf()
matshow(eigvec, fignum=False)
title('The eigen vectors')
colorbar()

# Get normalized poly coef
repr_coef = eigvec / sqrt(eigval)

figure(14); clf()
matshow(repr_coef, fignum=False)
title('The coefficients for the orthogonal polynomial')
colorbar()

figure(100); clf()
title('Well behavior at the interpolated point')
plot(s_x, dot(poly.polyvander(s_x, n_deg), repr_coef))

figure(102); clf()
title('Bad behavior otherwise - kind of bad generalization')
s_x_show = linspace(-1, 1, 500)
plot(s_x_show, dot(poly.polyvander(s_x_show, n_deg), repr_coef))


# vim: et sw=4 sts=4
