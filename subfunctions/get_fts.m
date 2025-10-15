function feat = get_fts(fnames, ii)

fs1 = 64;
frange = [0.5 2 ; 2 4 ; 4 8 ; 8 16 ; 16 32];
dr = 'L:\Lab_JamesR\Paediatric_FBA\Dataset\';

[dat, ~, label, fs, scle, ~] = read_edf([dr fnames{ii}]);
data_bp_mont = read_into_montage(dat, label, scle, fs, fs1, 1);
bn = floor(length(data_bp_mont)/(30*fs1))-1;
feat = zeros(bn,32);
for jj = 1:bn
  r1 = (jj-1)*30*fs1+1; r2 = r1+60*fs1-1;
  feat(jj,:) = median(single_channel_features(data_bp_mont(:, r1:r2), frange, fs1));
end
