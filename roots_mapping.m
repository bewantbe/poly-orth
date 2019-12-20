%

% exact mapping between roots of cheb and Legendre polynomials.

deg_max = 40;

%% all roots
%root_map = [];
%for k = 1:deg_max
%  rc = sort(cos(pi*(0.5+(0:k-1)')/k));
%  rl = sort(roots(LegendrePoly(k)));
%  root_map = [root_map; [rc rl]];
%end
%root_map = sortrows(root_map);

k = deg_max;
rc = sort(cos(pi*(0.5+(0:k-1)')/k));
rl = sort(roots(LegendrePoly(k)));
root_map = [rc rl];

figure(1);
plot(root_map(:,1), root_map(:,2), '-o')

% looks like this curve
f_ref_map = @(x, a) (1-x).^a .* log(1-x) - x.^a .* log(x);
a = 0.68;
ref_map = -0.97/(2^(2-a)*(a*log(2)-1)) * f_ref_map(root_map(:,1)/2+0.5, a);

figure(2);
plot(root_map(:,1), deg_max*(root_map(:,2) - root_map(:,1)), '-o',...
     root_map(:,1), ref_map, '-o')

figure(4);
plot(deg_max*(root_map(:,2) - root_map(:,1)), '-o')

%axis([-1 1 -0.5 0.5]*0.2)


%figure(3);
%xx = 0:0.001:1;
%a = 0.75;
%plot(xx, f_ref_map(xx, a) / (2^(2-a)+(a*log(2)-1)));
%axis([0.5+[-1 1]*0.1 [-1 1]*0.1])

%diff(f_ref_map(xx, a))(round(end/2)) / (xx(2)-xx(1))
