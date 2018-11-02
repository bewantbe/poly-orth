# Compare different orthogonal polynomials.

exec(open('imports_numpy.py').read())
exec(open('imports_plot.py').read())
exec(open('weighted_orthpoly_solver.py').read())

import numpy.polynomial as nppoly

#f = lambda x: abs(x)
#f = lambda x: np.maximum(0, x)**10
#f = lambda x: 1/(1+(x+0.01)**2)
#f = lambda x: sin(10*pi*x+0.1)
#f = lambda x: x**30
#f = lambda x: exp(x)
#f = lambda x: tanh(3*x+0.01)
#f = nppoly.legendre.Legendre.basis(20)
#f = nppoly.chebyshev.Chebyshev.basis(20)
f = lambda x: (x+1)*sin(1/(x+1))

#f_weighting = lambda x: 1
#f_weighting = lambda x: exp(-x*x)
f_weighting = lambda x: 1 + cos(pi*x)   # bad sampling at the boundary

n_deg = 100

#coef_c = nppoly.chebyshev.Chebyshev.interpolate(f, n_deg).coef
x_sample = sin( pi * (-n_deg/2 + np.arange(n_deg+1))/(n_deg+1) )
y_sample = f(x_sample)
coef_c = np.linalg.solve(cheb.chebvander(x_sample, n_deg), y_sample)

# Coefficients for normalized Legendre polynomials.
x_sample_leg = nppoly.legendre.legroots(1*(np.arange(n_deg+2)==n_deg+1))
y_sample_leg = f(x_sample_leg)
#coef_l = nppoly.legendre.legfit(x_sample_leg, y_sample_leg, n_deg) \
#         / sqrt(np.arange(n_deg+1)+0.5)
coef_l = np.linalg.solve(sqrt(np.arange(n_deg+1)+0.5) \
            * lege.legvander(x_sample_leg, n_deg), y_sample_leg)

# Coefficients for weighted orthogonal polynomials.
#x_sample_che = sin(pi * (np.arange(n_deg+1)-(n_deg)/2)/(n_deg+1))
x_sample_che = x_sample_leg
y_sample_che = f(x_sample_che)
basis_repr_coef = get_orthpoly(n_deg, f_weighting, 100)
coef_w = np.linalg.solve(dot(nppoly.chebyshev.chebvander(x_sample_che, n_deg), basis_repr_coef), y_sample_che)

#dot(dot(nppoly.chebyshev.chebvander(0.4, n_deg), basis_repr_coef), coef_w)

# Check accuracy of the generated polynomials.
if f_weighting(1/pi) == 1 and f_weighting(exp(-pi)) == 1:
    # f_weighting(x) == 1
    R = dot(nppoly.chebyshev.chebvander(x_sample_che, n_deg), basis_repr_coef)
    S = sqrt(arange(n_deg+1)+0.5) * lege.legvander(x_sample_che, n_deg)
    print('Maxabs polyvander =', maxabs(R-S))
    print('Cond              =', np.linalg.cond(S))

    n_b = 100
    ll_c = nppoly.chebyshev.Chebyshev(basis_repr_coef[:, n_b])
    ll   = nppoly.legendre.Legendre.basis(n_b) * sqrt(n_b+0.5)
#    ssx = np.linspace(-1,1,1000)
#    figure(323); clf()
#    plot(ssx, ll_c(ssx), ssx, ll(ssx))
    print("diff at 1 = %g" % (ll_c(1) - ll(1)))
    print("diff at 1 = %g" % (ll_c(0.1) - ll(0.1)))

coefs = np.vstack([
 coef_c,
 coef_l,
 coef_w,
 ])

coefs_v = np.vstack([
 orthpoly_coef(f, 'chebyshev', n_deg),
 orthpoly_coef(f, 'legendre', n_deg),
 orthpoly_coef(f, f_weighting, n_deg),
 ])

print('Verify coef =', maxabs(coefs - coefs_v), '(Must be zero)')

poly_name = ['Cheb', 'Leg', 'weight']

figure(10); clf()
plot(np.arange(n_deg+1), log10(abs(coefs.T)), '-o')
xlabel('k')
ylabel('log10(c_k)')
legend(poly_name)
plt.show()
