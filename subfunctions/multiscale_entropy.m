function msce = multiscale_entropy(dat)

dat = dat(:);
scale = 1:20; 
se = zeros(1, length(scale));
for z2 = 1:length(scale)
    dum = conv(dat, ones(scale(z2),1), 'same')./scale(z2);
    se(z2) = SampEn(2, 0.2*std(dum), dum);
end
msce(1) = se(1);
msce(2) = sum(se);
msce(3) = max(se);
msce(4) = sum(diff(se(1:5)));