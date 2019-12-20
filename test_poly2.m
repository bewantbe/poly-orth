%

x = chebfun('x');
%f = sin(1024*x);
f = exp(-1000*x^2);

% get(f, 'coeffs')

p = chebfun(f);

figure(1)
plot(f)

figure(2)
plot(f - p)

figure(3)
a = chebcoeffs(f, 'kind', 1);
semilogy(abs(a),'-o')

% A = CHEBCOEFFS(F, 'kind', 2)

