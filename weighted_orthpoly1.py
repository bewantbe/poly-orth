# Demo that generate the orthogonal polynomials of a customized weighting.
# In domain [-1, 1].

exec(open('common_imports_poly.py').read())

#f_weighting = lambda x: 1/sqrt(1-x*x)
f_weighting = lambda x: log(2.01)-log(x+1)
#f_weighting = lambda x: 1
#f_weighting = lambda x: exp(-x*x)
#f_weighting = lambda x: 1+cos(pi*x)

# Method 1: discrete sampling at good points (cheb points or Leb points)
n_deg = 100
n_x_sampling = n_deg + 1 + 10
x_sampling = sin( pi * (-(n_x_sampling-1)/2 + np.arange(n_x_sampling))/n_x_sampling )
# compare to xx, _ = cheb.chebgauss(n_deg+1)

V = cheb.chebvander(x_sampling, n_deg)

diag_weight = np.array([f_weighting(x)/cheb.chebweight(x) for x in x_sampling]) * (pi/n_x_sampling)

inner_prod_matrix = dot(V.T, diag_weight[:, np.newaxis] * V)

figure(10); clf()
matshow(inner_prod_matrix, fignum=False)
title('The inner product matrix.')
colorbar()

repr_coef = np.linalg.inv(np.linalg.cholesky(inner_prod_matrix).T)
print("Condition number %g."%(np.linalg.cond(inner_prod_matrix)))

#repr_coef = get_orthpoly(n_deg, f_weighting)
#repr_coef = inv(qr(sqrt(diag_weight[:, np.newaxis]) * V, mode='r'))

#eigval, eigvec = np.linalg.eigh(inner_prod_matrix)
#print("Condition number %g."%(eigval[-1]/eigval[0]))

## Get normalized poly coef
#repr_coef = eigvec / sqrt(eigval)

figure(14); clf()
matshow(repr_coef, fignum=False)
title('The coefficients for the orthogonal polynomial')
colorbar()

figure(102); clf()
s_x_show = linspace(-1, 1, 500)
plot(s_x_show, dot(cheb.chebvander(s_x_show, n_deg), repr_coef))

# Verify that they are orthogonal w.r.t. the sampling point
print("Should be machine eps %g" % (norm(dot(dot(V, repr_coef).T, diag_weight[:, np.newaxis]*dot(V, repr_coef)) - np.eye(n_deg+1))))

# it seems that sqrt(cond(inner_prod_matrix)) \approx max(poly)

