function [nf1, nf2] = noise_floor(data, anno, fs1);

nf1 = zeros(1,16); nf2 = nf1; wn = 32; wl = 8192; olap = 0.5;
f = linspace(0, fs1/2, wl/2);
for ii = 1:16
    dat = data(ii,:);
    dat(anno(ii,:)==1) = 0; 
    M = floor(length(dat)./wn); O = floor(M*olap);
    spc = zeros(2*wn-1, wl/2); 
    for jj = 1:2*wn-1
        r1 = (jj-1)*O+1; r2 = r1+M-1;
        X = fft([dat(r1:r2) zeros(1,wl-M)]);
        dum = X.*conj(X)./(wl-sum(anno(ii,r1:r2)));
        spc(jj,:) = dum(1:wl/2);
    end
    S = mean(spc);
    rf = find(f>70);
    nf1(ii) = mean(S(rf))./(fs1/2-70);  % noise region
    rf = find(f>0.5 & f<16);
    nf2(ii) = mean(S(rf))./(16-0.5); % cortical activity region

end