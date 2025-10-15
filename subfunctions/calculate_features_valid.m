function fts = calculate_features_valid(data, fs)

nni = diff(data)/fs;
nn = data./fs;
mnn = mean(nni);
mnni = nni-mean(nni);
sdnn = std(nni);  % Standard Devition of NN interval in seconds
[pxx,f] = plomb(mnni, nn(1:end-1), 1/(2*min(diff(nn)))); % lomb periodogram will estimate the PSD you can extract power from
pxx = 2*pxx*var(nni)/mean(pxx);
band0 = find(f <= 0.01); % define band
band1 = find(f > 0.01 & f <= 0.04); % define band
band2 = find(f > 0.04 & f <= 0.15); % define band
band3 = find(f > 0.15 & f <= 1); % define band
df1 = (0.04-0.01)/max(f);
df2 = (0.15-0.04)/max(f);
df3 = (1-0.15)/max(f);
df0 = 0.01;
bp0 = mean(pxx(band0))*df0/(mean(pxx))*sqrt(mean(pxx)/2); % SLF power
bp1 = mean(pxx(band1))*df1/(mean(pxx))*sqrt(mean(pxx)/2); % VLF band power
bp2 = mean(pxx(band2))*df2/(mean(pxx))*sqrt(mean(pxx)/2); % LF band power
bp3 = mean(pxx(band3))*df3/(mean(pxx))*sqrt(mean(pxx)/2); % HF power
se_samp = SampEn(2, 0.2*sdnn, nni); % Sample Entropy NN interval
[sd1, sd2] = poincare(nni); % Poincare Analysis NN interval
[fd,~,~] = Higuchi1Dn(mnni);
msce = multiscale_entropy(mnni);
[pse, ~] = pentropy(mnni,1:length(nni));
med_se = median(pse);
max_pse = max(pse);
k=PLmeasure(nni,1:length(nni));
[shannonEntr, renyiEntr] = shanRenEntropy(mnni,4,3);
%wt = cwt(nni,'amor'); inst_spec=trapz(median(abs(wt),1));
max_permEn = PermEn(mnni,50);
[lzc1, ~] = calc_lz_complexity(nni,'primitive', 1);
nni_sk = skewness(nni);
nni_ku = kurtosis(nni);
fuzz_en = FuzzyEn(mnni,2,0.2*sdnn,50);
maxmin_nni=max(nni)-min(nni);
area_nni=trapz(nni);
median_snleo = median(nlin_energy(mnni));
nni_centile=quantile(nni,[0.05 0.5 0.95]);

fts = [mnn*1000 sdnn*1000, bp0*1000, bp1*1000, bp2*1000, bp3*1000, se_samp, sd1, sd2, fd, msce, med_se, max_pse, k, shannonEntr, ...
    renyiEntr, max_permEn, lzc1, nni_sk,nni_ku, fuzz_en, maxmin_nni,area_nni,median_snleo, nni_centile]; 
end

