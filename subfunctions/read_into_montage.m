function [data_ref_mont, data_bp_mont, data_av_mont, data_csd_mont, stx, rr1] = read_into_montage(dat, label, scle, fs, fs1);
% Takes the data from the Read EDF file and then outputs 3 data types
% INPUTS
%   dat - EEG data read in by read_edf - cell array in integer16 format
%   label - labels of each channel recorded in dat - cell array of
%   characters
%   scle - the scaling parameter to convert 16-bit integers into microvolts
%   fs - sampling frequency of the data (as recorded)
%   fs1 - sampling frequency of the data (to be processed)
%   rr1 - is the reference channel, with respect to the referential montage
%
% OUTPUTS
% data_ref_mont = referential montage (reference is on the vertex somewhere Fz,Cz,Pz)
% data_bp_mont = bipolar montage (double banana)
% data_av_mont = average montage with a slight twist (corrected for phase lags - sort of)
% stx - the labels of the various montages (3 x cell array)
%
% Nathan Stevenson

str1{1} = 'Fp1'; str1{2} = 'Fp2'; str1{3} = 'F7'; str1{4} = 'F3'; 
str1{5} = 'Fz'; str1{6} = 'F4'; str1{7} = 'F8'; str1{8} = 'T3'; 
str1{9} = 'C3'; str1{10} = 'Cz'; str1{11} = 'C4'; str1{12} = 'T4'; 
str1{13} = 'T5'; str1{14} = 'P3'; str1{15} = 'Pz'; str1{16} = 'P4'; 
str1{17} = 'T6'; str1{18} = 'O1'; str1{19} = 'O2'; 
ref1 = zeros(1,length(dat));
for ii = 1:length(dat); ref1(ii) = length(strfind(label{ii}', str1{1})); end
qq0 = find(ref1==1,1);
data_ref_mont = zeros(length(str1), length(resample(zeros(1,length(dat{qq0(1)})), fs1, fs(qq0(1)))));
[B, A] = butter(4, 2*[45 55]/fs(1), 'stop'); % altered to have wider bandpass
for jj = 1:length(str1)
    ref1 = zeros(1,length(dat));
    for ii = 1:length(dat)
        ref1(ii) = length(strfind(label{ii}', str1{jj})); 
    end
    qq1 = find(ref1==1,1);
    if isempty(dat{qq1})==0
        data_ref_mont(jj,:) = resample(filter(B,A,double(dat{qq1})*scle(qq1)), fs1, fs(qq1)); 
    end
end

str = cell(18,2); 
str{1,1} = 'Fp2'; str{1,2} = 'F4';     % Fp2-F4, Fp2-F4
str{2,1} = 'F4'; str{2,2} = 'C4';     % F4-C4, F4-C4
str{3,1} = 'C4'; str{3,2} = 'P4';     % C4-P4, C4-P4
str{4,1} = 'P4'; str{4,2} = 'O2';     % P4-O2, P4-O2
str{5,1} = 'Fp1'; str{5,2} = 'F3';     % Fp1-F3, Fp1-F3
str{6,1} = 'F3'; str{6,2} = 'C3';     % F3-C3, F3-C3
str{7,1} = 'C3'; str{7,2} = 'P3';     % C3-P3, C3-P3
str{8,1} = 'P3'; str{8,2} = 'O1';     % P3-O1, P3-O1
str{9,1} = 'Fp2'; str{9,2} = 'F8';     % Fp2-F8, Fp2-F8
str{10,1} = 'F8'; str{10,2} = 'T4';     % F8-T4, F8-T8
str{11,1} = 'T4'; str{11,2} = 'T6';     % T4-T6, T8-P8
str{12,1} = 'T6'; str{12,2} = 'O2';     % T6-O2, P8-O2
str{13,1} = 'Fp1';  str{13,2} ='F7';     % Fp1-F7, Fp1-F7
str{14,1} = 'F7'; str{14,2} = 'T3';     % F7-T3, F7-T7
str{15,1} = 'T3'; str{15,2} = 'T5';     % T3-T5, T7-P7
str{16,1} = 'T5'; str{16,2} = 'O1';     % T5-O1, P7-O1
str{17,1} = 'Fz'; str{17,2} = 'Cz';     % Fz-Cz, Fz-Cz
str{18,1} = 'Cz';  str{18,2} ='Pz';     % Cz-Pz, Cz-Pz
A = size(data_ref_mont);
data_bp_mont = zeros(A(1)-1, A(2));
for jj = 1:18
    ref1 = zeros(1,length(str1));
    ref2 = zeros(1,length(str1));
    for ii = 1:length(str1)
        ref1(ii) = length(strfind(str1{ii}, str{jj,1})); 
        ref2(ii) = length(strfind(str1{ii}, str{jj,2}));
    end
    qq1 = find(ref1==1,1);
    qq2 = find(ref2==1,1);
    if length(qq1)==length(qq2)
        data_bp_mont(jj,:) = data_ref_mont(qq1,:)-data_ref_mont(qq2,:); 
    end
end

rr = sum(abs(data_ref_mont),2);
rr1 = find(rr==min(rr));
rr2 = find(rr~=min(rr));
mref = hilbert(mean(data_ref_mont(rr2,:))); % this assumption is that the average here is an estimate of ref
mref = mref./norm(mref);
data_av_mont = zeros(A);
for ii = 1:A(1)
    w = (data_ref_mont(ii,:)*(mref'));
    data_av_mont(ii,:) = data_ref_mont(ii,:)-2*real(mref'*w)';  
end

str2 = str1; % convert to EGI Geodesics 129 format if necessary
str2{8} = '46'; str2{12} = '109'; str2{13} = '58'; str2{17} = '97';
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',str2');
[G,H] = GetGH(M);
data_csd_mont = CSD([data_ref_mont(1:rr1-1,:) ; zeros(1,length(data_ref_mont)) ; data_ref_mont(rr1+1:end,:)], G, H);

stx = cell(A(1), 4);
for ii = 1:A(1); stx{ii,1} = [str1{ii} ' - REF']; stx{ii,3} = [str1{ii} ' - Av']; stx{ii,4} = str1{ii}; end
for ii = 1:A(1)-1; stx{ii,2} = [str{ii,1} ' - ' str{ii,2}]; end

% figure; hold on;
% for ii = 1:19
%     plot(data_ref_mont(ii,:)+200*ii);
% end
