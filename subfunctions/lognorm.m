function y = lognorm(beta, x);

k = beta(1);
m = beta(2);
s = beta(3);
o = beta(4);
y = (k*exp(-(x-m).^2/(2*s.^2)))+o;
