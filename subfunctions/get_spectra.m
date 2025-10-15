function [spc3a, spc3b] = get_spectra(fnames, ii)

fs1 = 64;
dr = 'L:\Lab_JamesR\Paediatric_FBA\Dataset\';
[dat, ~, label, fs, scle, ~] = read_edf([dr fnames{ii}]);
data_bp_mont = read_into_montage(dat, label, scle, fs, fs1, 1);
A = size(data_bp_mont);
bn = floor(length(data_bp_mont)/(30*fs1))-1;
%Pre-N2
spc2 = zeros(10, 257);
for jj = 1:10
  r1 = (jj-1)*30*fs1+1; r2 = r1+60*fs1-1;
  spc1 = zeros(A(1), 257);
  for zz = 1:A(1)
    dat = data_bp_mont(zz, r1:r2);
    Pxx = pwelch(dat, hamming(512),256);
    spc1(zz,:) = Pxx';
  end
  spc2(jj,:) = mean(spc1);
end
spc3a = median(spc2);

% Post N2
spc2 = zeros(bn-10, 257);
for jj = 11:bn
  r1 = (jj-1)*30*fs1+1; r2 = r1+60*fs1-1;
  spc1 = zeros(A(1), 257);
  for zz = 1:A(1)
    dat = data_bp_mont(zz, r1:r2);
    Pxx = pwelch(dat, hamming(512),256);
    spc1(zz,:) = Pxx';
  end
  spc2(jj-10,:) = mean(spc1);
end
spc3b = median(spc2);