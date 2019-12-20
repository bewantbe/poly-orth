%

K = 100;

rc = chebpts(K, 1);
rl = legpts(K);

figure(34);
plot(rl - rc, '-o');

geomean = @(x) exp(mean(log(x)));

u = [];
for k = 1 : length(rc)
  u(k) = geomean(abs(rc((1:K)~=k) - rc(k)));
end

figure(34);
plot(rc, u);

