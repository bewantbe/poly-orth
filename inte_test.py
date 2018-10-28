# 

def integral(f, deg = 50):
    s_xw = cheb.chebgauss(deg)
    return sum([f(x) / cheb.chebweight(x) * w for x, w in np.stack(s_xw).T])

f_weighting = lambda x: 1/sqrt(1-x*x)
integral(f_weighting)

f_weighting = lambda x: exp(-x*x)
integral(f_weighting, 500)

f_weighting = lambda x: 1
integral(f_weighting)

f_weighting = lambda x: x*x
integral(f_weighting)

