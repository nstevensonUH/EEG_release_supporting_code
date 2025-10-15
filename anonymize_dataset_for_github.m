

% This m-file finds all EDF files in the current working directory and then
% exports them into a new EDF file with random ID filename (4 numbers and 2 letters)
% and ensures the AGE variable is quantised (nonlinearly) so that no less than 10
% subjects falls within each AGE value

% CHANGE TO ENSURE IT MATCHES YOUR EQUIVALENT PATHS
addpath(genpath('D:\'))

% THIS FILE CONTAINS THE FILENAME ID, AGE AND SEX OF EACH SUBJECT
load age_filenames_v1_sex
sex2(find(isnan(sex2)==1))=1; % checked subject is female

% Read all EDF files in a directory
dir_name = 'D:\Anonymization\Before_EDF'; %test location on local machine (Helsinki)
cd(dir_name)
aa = dir('**\*.edf');
aa = {aa.name};

%EXTRACT AND SYNCHRONIZE DATA AND DEMOGRAPHIC INFO
age = zeros(1, length(aa)); sex = zeros(1, length(aa)); len = cellfun(@length, fns2);
for ii = 1:length(aa)
    val = contains(fns2, aa{ii}(1:end-4)) & (len == length(aa{ii})-4);
    age(ii) = age2(val==1);
    sex(ii) = sex2(val==1);
end

% QUANTIZE AGE - nonuniform (the age resolution changes with population density);
yref = 10; rng1 = []; ageqd = []; nage = age; c1=1;
while length(nage)>11
    val = unique(nage); hst = zeros(1,length(val));
    for ii = 1:length(val); hst(ii) = sum(nage==val(ii)); end
    chst = cumsum(hst);
    rng1(c1,2) = val(find(chst>yref,1));
    ageqd(c1) = mean(nage(nage<=rng1(c1,2)));
    nage(nage<=rng1(c1,2))=NaN;
    nage = nage(~isnan(nage));
    c1 = c1+1;
end
rng1(2:end,1) = rng1(1:end-1,2); rng1(end,2) = max(val);

xage = age; c1 = 1; ageqd = []; val1 = [];
while length(xage)>10
    [nage, idx] = sort(xage, 'descend');
    val1(c1,1) = nage(1); val1(c1,2) = nage(10);
    ageqd(c1) = mean(xage(xage>=val1(c1,2)));
    xage = xage(xage<val1(c1,2));
    c1 = c1+1;
end
val1(c1-1,2)= min(xage);
ageqd(c1-1) = mean(age(age<val1(c1-1,1) & age>=val1(c1-1,2)));

% TEST
% b = zeros(1, length(rng1));
% for ii = 1:length(rng1)    
%     b(ii) = length(find(age>rng1(ii,1) & age<=rng1(ii,2)));
% end
% ASSIGN QUANTISED AGES TO DATA
ageq = zeros(1, length(age)); rng1 = val1;
for ii = 1:length(rng1)    
    %ref = find(age>rng1(ii,1) & age<=rng1(ii,2));
    ref = find(age<=val1(ii,1) & age>=val1(ii,2));
    ageq(ref) = ageqd(ii); 
end

% generate random ID - 4 numerals and 2 letters randomly distributed in the
% filename
rng('default')
val = round(rand(1)*1000)+1;
rng(val); fname = {};
while length(fname)<length(aa)   
    str = [num2str(round(rand(1)*9)) num2str(round(rand(1)*9)) num2str(round(rand(1)*9)) num2str(round(rand(1)*9)) char(round(rand(1,2)*25+65))];
    k = randsample(1:6, 6);
    str = str(k);
    if length(unique([fname ; str]))>length(fname)
        fname = [fname ; str];
    end
end

% These are the channels to save
str1{1} = 'Fp1'; str1{2} = 'Fp2'; str1{3} = 'F7'; str1{4} = 'F3'; 
str1{5} = 'Fz'; str1{6} = 'F4'; str1{7} = 'F8'; str1{8} = 'T3'; 
str1{9} = 'C3'; str1{10} = 'Cz'; str1{11} = 'C4'; str1{12} = 'T4'; 
str1{13} = 'T5'; str1{14} = 'P3'; str1{15} = 'Pz'; str1{16} = 'P4'; 
str1{17} = 'T6'; str1{18} = 'O1'; str1{19} = 'O2'; 
str2{1} = 'EKG'; str2{2} = 'ECG'; 

% CHECK BEFORE PROCEEDING
if length(unique(fname)) == length(aa)
    disp('reading EDFs and writing to de-identified format')

    % THIS LOOP GOES THROUGH ALL FILES AN RE-EXPORTS WITH NEW FILENAME AND
    % DE-IDENTIFIED HEADER (which they already had) and only the EEG and ECG
    % parts of the signal
    for ii = 1:length(aa)
        [dat, hdr_old, label, fs, scle, offs] = read_edf(aa{ii});
        % ONLY EXPORT USEFUL CHANNELS WHICH HERE WILL BE 10-20 ELECTRODES
        % PLUS EKG
        for jj = 1:length(label); label{jj} = label{jj}'; end
        ref = zeros(1,length(str1));
        for jj = 1:length(str1)
            dum = contains(label, str1{jj});
            ref(jj) = find(dum==1,1);
        end
        rval = contains(label, str2{1}) | contains(label, str2{2});  
        if sum(rval)>0
            ref(length(str1)+1) = find(rval==1, 1);
        end

        datf = int16(zeros(length(ref), length(dat{ref(1)})));
        for jj = 1:length(ref)
            datf(jj,:) = dat{ref(jj)};
        end
    
    % THIS PART OF THE HEADER MAY CONTAIN IDENTIFIABLE INFORMATION SO
    % RE-WRITE WITH GENERIC UNIDENTIFIABLE TEXT - this is painful but
    % procedural
    len = length(dat{1})/fs(1);
    hdr_new = ['0' ; char(32*ones(7,1)) ; ['FBA in Pediatrics Study']' ; char(32*ones(160-23,1)) ; ['01.01.01']' ; ['01.01.01']' ; ['1792    ']' ; char(32*ones(44,1)) ; ...
        num2str(len)' ; char(32*ones(8-length(num2str(len)),1)) ; '1' ; char(32*ones(7,1)) ; num2str(length(ref))' ; char(32*ones(4-length(num2str(length(ref))),1))];
    rx8 = []; rx16 = []; rx32 = []; rx80 = [];
    for jj = 1:length(ref)
        r1 = (ref(jj)-1)*8+1; r2 = r1+7;
        rx8 = [rx8 r1:r2];
        r1 = (ref(jj)-1)*16+1; r2 = r1+15;
        rx16 = [rx16 r1:r2];
        r1 = (ref(jj)-1)*32+1; r2 = r1+31;
        rx32 = [rx32 r1:r2];
        r1 = (ref(jj)-1)*80+1; r2 = r1+79;
        rx80 = [rx80 r1:r2];
    end
    hdr_new = [hdr_new ; hdr_old{2}(rx16)];
    hdr_new = [hdr_new ; hdr_old{3}(rx80)];
    hdr_new = [hdr_new ; hdr_old{4}(rx8)];
    hdr_new = [hdr_new ; hdr_old{5}(rx8)];
    hdr_new = [hdr_new ; hdr_old{6}(rx8)];
    hdr_new = [hdr_new ; hdr_old{7}(rx8)];
    hdr_new = [hdr_new ; hdr_old{8}(rx8)];
    hdr_new = [hdr_new ; hdr_old{9}(rx80)];
    hdr_new = [hdr_new ; hdr_old{10}(rx8)];
    hdr_new = [hdr_new ; hdr_old{11}(rx32)];

    h1 = [num2str(length(hdr_new))' ; char(32*ones(8-length(num2str(length(hdr_new))),1))];
    hdr_new(185:192) = h1;

    % WRITE NEW DATA with randomized filename
    fdum = ['D:\Anonymization\After_EDF\' fname{ii} '.edf'];
    dum1 = write_edf(fdum, datf, hdr_new, length(ref), fs(ref(1)));
    fclose('all');

    end

    % WRITE CODEX AND ADJUSTED AGES TO SPREADSHEET
    dum = cell(length(fname)+1, 3);
    dum{1,1} = 'Anon. Filename';
    dum{1,2} = 'Quantized Age (years)'; dum{1,3} = 'Quantized Sex';
    for ii = 1:length(fname)
        dum{ii+1,1} = [fname{ii} '.edf'];
        dum{ii+1,2} = ageq(ii); 
        dum{ii+1,3} = sex(ii);
    end
    xlswrite('D:\Anonymization\anon_codex.xls', dum);
    writecell(dum, 'D:\Anonymization\anon_codex.csv')
else
    disp('Not enough unique filenames')
end


