%
addpath('/home/xyy/Documents/research/ann/codes/chebfun-master');

%f = chebfun('exp(4*x)');
%f = chebfun('tanh(x)');
f = chebfun('sqrt(x+2)');
plot(roots(f, 'all'), 'o');

