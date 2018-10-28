% Demo that generate the orthogonal polynomials of a customized weighting.
% In domain [-1, 1].
% Method 2: use Gauss-Chebyshev quadrature for the inner product.

addpath('/home/xyy/Documents/research/ann/codes/chebfun-master/')

%f_weighting = chebfun('1');
f_weighting = chebfun('exp(-x*x)');
%f_weighting = chebfun('1+cos(x*pi)');  % is this work???

n_deg = 10;
inner_prod_matrix = zeros(n_deg+1, n_deg+1);

for ii = 0:n_deg
  for jj = 0:n_deg
    if jj > ii
      continue
    end
    v = sum(f_weighting * chebpoly(ii) * chebpoly(jj));
    inner_prod_matrix(ii+1, jj+1) = v;
    inner_prod_matrix(jj+1, ii+1) = v;
  end
end

figure(34)
imagesc(inner_prod_matrix)

repr_coef = inv(chol(inner_prod_matrix));

%[eigvec, eigval] = eig(inner_prod_matrix);
%repr_coef = eigvec ./ sqrt(diag(eigval)).';  % this not work due to not ordered

poly_w = zeros(1, n_deg+1);
for ii = 0:n_deg
  poly_w = poly_w + chebpoly(ii) * repr_coef(ii+1, :);
end

verify_orth = zeros(n_deg+1);
for ii = 0:n_deg
  verify_orth(ii+1, :) = sum(f_weighting * poly_w(:,ii+1) * poly_w);
end
fprintf('Should be machine eps: %g\n', norm(verify_orth - eye(n_deg+1)));

figure(35);
plot(poly_w)
xlim([-1,1])

