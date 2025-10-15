% Quality Assessment Open EEG Data


% READ IN CSV file with filenames, quantized age and sex
filename = 'anon_codex_age_sex.csv';
% Specify range and delimiter
dataLines = [2, Inf];
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = dataLines;
opts.Delimiter = ",";
opts.VariableNames = ["AnonFilename", "QuantizedAgeyears", "QuantizedSex"];
opts.VariableTypes = ["string", "double", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "AnonFilename", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["AnonFilename", "QuantizedSex"], "EmptyFieldRule", "auto");
% Import the data
anoncodex = readtable(filename, opts);
fnames = cellstr(anoncodex{:,1});
ages = anoncodex{:,2};
sex = anoncodex{:,3};
s1 = contains(sex, 'M'); % M - Male
sex = s1; 
%s2 = contains(sex, 'F'); % F - Female

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EEG ANALYSIS
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERFORM BASIC EVALUATIONS ON DATASET
fs1 = 250;
val1a = zeros(1, length(fnames)); val1b = val1a; val2 = val1a;
for ii = 1:length(fnames)
    [dat, hdr, label, fs, scle, offs] = read_edf(fnames{ii});
    data_bp_mont = read_into_montage(dat, label, scle, fs, fs1, 0);
    anno = per_ha_art(data_bp_mont, fs1);
    val2(ii) = sum(sum(anno))/prod(size(anno))*100;
    [nf1, nf2] = noise_floor(data_bp_mont, anno, fs1);
    val1a(ii) = mean(nf1); val1b(ii) = mean(nf2); 
end
save('results_qa_v1.mat', 'val1a', 'val1b', 'val2');

load results_qa_v1
X = [ages sex]; % convert sex to binary first
y1 = fitlm(X, val1a');
y2 = fitlm(X, val1b');

% SIMPLE FEATURE EXTRACTION
fts = cell(1, length(fnames)); 
%parpool(8) consider using parpool if possible
for ii = 1:length(fnames)
    ii
    fts{ii} = get_fts(fnames, ii);
end
%cd('L:\Lab_JamesR\nathanST\FBA_paeds_project\Documents\paper3_open_data')
save('results_qa_feats_v1.mat', 'fts', 'ages', 'sex')

load results_qa_feats_v1 
% CORRELATIONS WITH AGE
fv1 = zeros(length(fts), 32); fv2 = fv1;
for ii = 1:length(fts)
    fval = fts{ii};
    fv1(ii,:) = median(fval(1:10,:));
    A = size(fval);
    if A(1)<21
        fv2(ii,:) = median(fval(11:end,:)); % pre N2 onset
    else
        fv2(ii,:) = median(fval(11:21,:)); % post N2 onset
    end
end
p1 = corr(fv1, ages);
p2 = corr(fv2, ages);

val = unique(ages);
prs = [1:2:length(val) ; 2:2:length(val)];
x1 = ages(sex==0);
x2 = ages(sex==1);
xx = linspace(min(ages), max(ages), 100);
A = size(fv2);
for ii = 1:A(2)
    ii
    y1 = fv2(sex==0,ii);
    y2 = fv2(sex==1,ii);   
    y1v = zeros(1, length(val)); y2v = y1v;
    for jj = 1:length(val)
        y1v(jj) = median(y1(x1==val(jj)));
        y2v(jj) = median(y2(x2==val(jj)));
    end
    cohD = zeros(length(prs), 5);
    for jj = 1:length(prs)
        dum = meanEffectSize(y1(x1==val(prs(1,jj)) | x1==val(prs(2,jj))), y2(x2==val(prs(1,jj)) | x2==val(prs(2,jj))));
        cohD(jj,:) = [dum.Variables sum(x1==val(prs(1,jj)) | x1==val(prs(2,jj))) sum(x2==val(prs(1,jj)) | x2==val(prs(2,jj)))];
    end
end

% START PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,2,1);
f = zeros(1,17);
for ii = 1:17
    rf = find(ages>=ii-1 & ages<ii);
    f(ii) = median(fv1(rf,25));    
end
x = [1:17]-0.5;
B = polyfit(x([1:15 17]), f([1:15 17]),2);
res = fv1(:,25)-polyval(B,ages);
rx = find(abs(res)<5e-4);
plot(ages(rx), fv1(rx,25), '.');
ylabel('Hjorth 2')
xlabel('Age (y)')
axis([-1 17 0 2e-3])
set(gca, 'position', [0.1 0.575 0.375 0.375], 'FontName', 'times', 'FontSize', 12)
title('Pre- N2 onset')
text(0, 1.75e-3, 'A', 'FontName', 'times', 'FontSize', 14)
pp = corr(ages, fv1(:,25), 'type','Spearman');
text(10, 0.25e-3, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 12)

subplot(2,2,3);
f = zeros(1,17);
for ii = 1:17
    rf = find(ages>=ii-1 & ages<ii);
    f(ii) = median(fv1(rf,30));    
end
x = [1:17]-0.5;
B = polyfit(x([1:15 17]), f([1:15 17]),2);
res = fv1(:,30)-polyval(B,ages);
rx = find(abs(res)<0.02);
plot(ages(rx), fv1(rx,30), '.');
ylabel('Burst Sharpness')
xlabel('Age (y)')
axis([-1 17 -.97 -0.89])
set(gca, 'position', [0.1 0.1 0.375 0.375], 'FontName', 'times', 'FontSize', 12)
text(15, -0.9, 'B', 'FontName', 'times', 'FontSize', 14)
pp = corr(ages, fv1(:,30), 'type','Spearman');
text(0, -0.96, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 12)

subplot(2,2,2);
f = zeros(1,17);
for ii = 1:17
    rf = find(ages>=ii-1 & ages<ii);
    f(ii) = median(fv2(rf,16));    
end
x = [1:17]-0.5;
B = polyfit(x([1:15 17]), f([1:15 17]),3);
res = fv2(:,16)-polyval(B,ages);
rx = find(abs(res)<5);
plot(ages(rx), fv2(rx,16), '.');
ylabel('Relative Alpha* Power (%)')
xlabel('Age (y)')
axis([-1 17 0 17.5])
set(gca, 'position', [0.6 0.575 0.375 0.375], 'FontName', 'times', 'FontSize', 12)
title('Post- N2 onset')
text(0, 15, 'C', 'FontName', 'times', 'FontSize', 14)
pp = corr(ages, fv2(:,16), 'type','Spearman');
text(10, 2, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 12)

subplot(2,2,4);
f = zeros(1,17);
for ii = 1:17
    rf = find(ages>=ii-1 & ages<ii);
    f(ii) = median(fv2(rf,1));    
end
x = [1:17]-0.5;
B = polyfit(x([1:15 17]), f([1:15 17]),2);
res = fv2(:,1)-polyval(B,ages);
rx = find(abs(res)<4);
plot(ages, fv2(:,1), '.');
ylabel('Amplitude (5^{th} centile; \muV)')
xlabel('Age (y)')
axis([-1 17 0 15])
set(gca, 'position', [0.6 0.1 0.375 0.375], 'FontName', 'times', 'FontSize', 12)
text(15, 13, 'D', 'FontName', 'times', 'FontSize', 14)
pp = corr(ages, fv1(:,1), 'type','Spearman');
text(0, 2, ['r=' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 12)

% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRAIN GPR based age prediction - 5-fold CV 80:20 split

% Find outliers in features
rfs = zeros(1,length(fv1));
for jj = 1:32
    y = fv1(:,jj); x = ages;
    B = polyfit(x,y,2);
    res = y - polyval(B, x);
    dum = zeros(1, length(y)); dum(find(abs(res)>6*std(res)))=1;
    rfs = rfs+dum;
    y = fv2(:,jj); x = ages;
    B = polyfit(x,y,2);
    res = y - polyval(B, x);
    dum = zeros(1, length(y)); dum(find(abs(res)>6*std(res)))=1;
    rfs = rfs+dum;
end

% get rid of outliers for now
nref = find(rfs==0);
fv1 = fv1(nref,:);
fv2 = fv2(nref,:);
ages = ages(nref);

% split up data into 5 folds for cross-validation
pid_in1 = 1:length(fv1); 
pma_in1 = ages;
Mc = 5;
out = zeros(1,1000);
for qq = 1:1000
    rng(qq)
    pd = unique(pid_in1);
    y = cell(1,Mc);
    rx = ones(1,length(pd)); 
    K = floor(length(pd)/Mc); dum = rem(length(pd), Mc);
    K = K.*ones(1,Mc); K(1:dum) = K(1)+1;
    for ii = 1:Mc
       rz = find(rx==1);
       dum = randsample(length(rz), K(ii), false);
       rx(rz(dum)) = 0;
       y{ii} = rz(dum);
    end
    y{end} = [y{end}' ; find(rx==1)']';
    for ii = 1:Mc
        y{ii} = pd(y{ii});
    end
    ss = cell(1, Mc);
    for ii = 1:Mc
        rf = y{ii};
        ag = [];
        for jj = 1:length(rf); ag = [ag pma_in1(find(pid_in1==rf(jj)))]; end
        ss{ii} = ag;
    end
    pv = zeros(1,6); c1 = 1;
    for ii = 1:Mc
        for jj = ii+1:Mc
            [~, pv(c1)] = kstest2(ss{ii}, ss{jj});
            c1 = c1+1;
        end
    end
    out(qq) = mean(pv);
end
nr = find(out==max(out)); % nr = 7 
rng(nr)
yr = cell(1,Mc);
rx = ones(1,length(pd)); %rt = 1-rx; %rr = 1:length(pid);
for ii = 1:Mc
   rz = find(rx==1);
   dum = randsample(length(rz), K(ii), false);
   rx(rz(dum)) = 0;
   yr{ii} = rz(dum);
end
yr{end} = [yr{end}' ; find(rx==1)']';
for ii = 1:Mc; yr{ii} = pd(yr{ii}); end

% do training and testing of age prediction for both pre- and post- N2
% onset data
pred1 = []; pred2 = []; agx = [];
for ii = 1:Mc
    ii
    dum = zeros(1,length(fv1));
    dum(yr{ii})=1;
    r1 = find(dum==1); r2 = find(dum==0);

    rGP1 = fitrgp(fv1(r2,:), ages(r2), 'BasisFunction', 'constant', 'KernelFunction', 'matern52', 'Standardize', true);
    out1 = predict(rGP1, fv1(r1,:));
    pred1 = [pred1 ; out1]; 

    rGP2 = fitrgp(fv2(r2,:), ages(r2), 'BasisFunction', 'constant', 'KernelFunction', 'matern52', 'Standardize', true);
    out2 = predict(rGP1, fv2(r1,:));
    pred2 = [pred2 ; out2]; 

    agx = [agx ; ages(r1)];

end

% START PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(1,2,1)
plot(agx, pred1 ,'.')
title('Pre- N2 Onset')
set(gca, 'FontName', 'times', 'Fontsize', 16)
xlabel('Age (y)'); ylabel('Predicted Age (y)')
grid on; hold on;
plot([-2 18], [-2, 18], 'k')
axis([-1 17 -1, 17])
set(gca, 'position', [0.075 0.15 0.4 0.8], 'Xtick', [0:2:16], 'Ytick', [0:2:16])
%text(0, 16, 'A', 'FontName', 'times', 'Fontsize', 14)
pp = corr(agx, pred1);
text(10, 2, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 16)

subplot(1,2,2)
plot(agx, pred2 ,'.')
title('Post- N2 Onset')
set(gca, 'FontName', 'times', 'Fontsize', 16)
xlabel('Age (y)'); ylabel('Predicted Age (y)')
grid on; hold on;
plot([-2 18], [-2, 18], 'k')
axis([-1 17 -1, 17])
set(gca, 'position', [0.575 0.15 0.4 0.8], 'Xtick', [0:2:16], 'Ytick', [0:2:16])
%text(0, 16, 'B', 'FontName', 'times', 'Fontsize', 14)
pp = corr(agx, pred2);
text(10, 2, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 16)

% END PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DO SPECTRAL ANALYSIS
addpath('L:\Lab_JamesR\nathanST\FBA_paeds_project\Code\fba_code_107052022\subfunctions')
spa = cell(1, length(fnames)); spb = spa;
%parpool(6)
for ii = 1:length(fnames)
    ii
    [spa{ii}, spb{ii}] = get_spectra(fnames, ii);
end
save('results_qa_spec_v1.mat', 'spa', 'spb', 'ages', 'sex')

load results_qa_spec_v1 

rn = [0 1 ; 1 3 ; 3 7 ; 7 17];
d1 = 10*log10(cell2mat(spa')); d2 = 10*log10(cell2mat(spb'));
spec_pre_mean = zeros(4, 257); spec_pre_std = spec_pre_mean; spec_post_mean = spec_pre_mean; 
spec_post_std = spec_pre_mean;
for ii = 1:4
    rf = find(ages> rn(ii,1) & ages <= rn(ii,2));
    length(rf)
    spec_pre_mean(ii,:) = mean(d1(rf,:));
    spec_pre_std(ii,:) = std(d1(rf,:));
    spec_post_mean(ii,:) = mean(d2(rf,:));
    spec_post_std(ii,:) = std(d2(rf,:));
end

% START PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; set(gcf, 'Position', [1000 150 800 680])
subplot(4,4,[1 2 5 6]); hold on;
f = linspace(0, 32, 257); rf = 6:255;
for ii = 1:4
    h(ii) = plot(f(rf), spec_pre_mean(ii,rf), 'LineWidth', 2);
end
title('Pre- N2 Onset')
axis([0 32 -10 40])
legend(h, '0-1y','1-3y','3-7y','7-16y')
set(gca, 'Fontname', 'times', 'fontsize', 14, 'Xtick', [0:4:32], 'XTickLabel', [0:4:32])
xlabel('Frequency (Hz)', 'Fontname', 'times', 'fontsize', 16);
ylabel('Power (dB)', 'Fontname', 'times', 'fontsize', 16);
grid on
set(gca, 'Position', [0.075 0.575 0.405 0.4])
%text(4, 35, 'A', 'Fontname', 'times', 'fontsize', 18)

fr = [9 10 13 14];
pst1(1,:) = [0.075 0.285 0.19 0.2];
pst1(2,:) = [0.295 0.285 0.19 0.2];
pst1(3,:) = [0.075 0.075 0.19 0.2];
pst1(4,:) = [0.295 0.075 0.19 0.2];
for ii = 1:4
    subplot(4,4, fr(ii))
    cc = get(h(ii), 'color');
    g1 = spec_pre_mean(ii,:)+2*spec_pre_std(ii,:);
    g2 = spec_pre_mean(ii,:)-2*spec_pre_std(ii,:);
    j = patch([f(rf(1)) f(rf) f(rf(250:-1:1)) f(rf(1))], [g2(rf(1)) g1(rf) g2(rf(250:-1:1)) g1(rf(1))], cc);
    set(j, "EdgeColor", 'none')
    hold on;
    plot(f(rf), spec_pre_mean(ii,rf), 'color', [1 1 1]);
    axis([0 32 0 40])
    set(gca, 'Position', pst1(ii,:),'Xtick', [0:8:32], 'Fontname', 'times', 'fontsize', 14)
    switch ii
        case 1
            set(gca, 'Xticklabel', {}, 'Fontname', 'times', 'fontsize', 14)
            ylabel('Power (dB)', 'Fontname', 'times', 'fontsize', 16);
            text(20, 32, '0-1y', 'Fontname', 'times', 'fontsize', 14)
        case 2
            set(gca, 'Yticklabel', {}, 'Xticklabel', {}, 'Fontname', 'times', 'fontsize', 14)
            text(20, 32, '1-3y', 'Fontname', 'times', 'fontsize', 14)
        case 3
            set(gca, 'Fontname', 'times', 'fontsize', 14)
            xlabel('Frequency (Hz)', 'Fontname', 'times', 'fontsize', 16);
            ylabel('Power (dB)', 'Fontname', 'times', 'fontsize', 16);
            text(20, 32, '3-7y', 'Fontname', 'times', 'fontsize', 14)
        case 4    
            set(gca, 'Yticklabel', {}, 'Fontname', 'times', 'fontsize', 14)
            xlabel('Frequency (Hz)', 'Fontname', 'times', 'fontsize', 16);
            text(20, 32, '7-16y', 'Fontname', 'times', 'fontsize', 14)
    end
    grid on;
    axis([0 32 -10 40])
end

subplot(4,4,[3 4 7 8]); 
hold on;
f = linspace(0, 32, 257); rf = 6:255;
for ii = 1:4
    h1(ii) = plot(f(rf), spec_post_mean(ii,rf), 'LineWidth', 2);
end
title('Post- N2 Onset')
axis([0 32 -10 40])
legend(h1, '0-1y','1-3y','3-7y','7-16y')
set(gca, 'Fontname', 'times', 'fontsize', 14, 'Xtick', [0:4:32], 'XTickLabel', [0:4:32])
xlabel('Frequency (Hz)', 'Fontname', 'times', 'fontsize', 16);
ylabel('Power (dB)', 'Fontname', 'times', 'fontsize', 16);
grid on
set(gca, 'Position', [0.575 0.575 0.405 0.4])
%text(4, 35, 'B', 'Fontname', 'times', 'fontsize', 18)

fr = [9 10 13 14];
pst1(1,:) = [0.575 0.285 0.19 0.2];
pst1(2,:) = [0.795 0.285 0.19 0.2];
pst1(3,:) = [0.575 0.075 0.19 0.2];
pst1(4,:) = [0.795 0.075 0.19 0.2];
for ii = 1:4
    subplot(4,4, fr(ii))
    cc = get(h1(ii), 'color');
    g1 = spec_post_mean(ii,:)+2*spec_post_std(ii,:);
    g2 = spec_post_mean(ii,:)-2*spec_post_std(ii,:);
    j = patch([f(rf(1)) f(rf) f(rf(250:-1:1)) f(rf(1))], [g2(rf(1)) g1(rf) g2(rf(250:-1:1)) g1(rf(1))], cc);
    set(j, "EdgeColor", 'none')
    hold on;
    plot(f(rf), spec_post_mean(ii,rf), 'color', [1 1 1]);
    axis([0 32 0 40])
    set(gca, 'Position', pst1(ii,:),'Xtick', [0:8:32], 'Fontname', 'times', 'fontsize', 14)
    switch ii
        case 1
            set(gca, 'Xticklabel', {}, 'Fontname', 'times', 'fontsize', 14)
            ylabel('Power (dB)', 'Fontname', 'times', 'fontsize', 16);
            text(20, 32, '0-1y', 'Fontname', 'times', 'fontsize', 14)
        case 2
            set(gca, 'Yticklabel', {}, 'Xticklabel', {}, 'Fontname', 'times', 'fontsize', 14)
            text(20, 32, '1-3y', 'Fontname', 'times', 'fontsize', 14)
        case 3
            set(gca, 'Fontname', 'times', 'fontsize', 14)
            xlabel('Frequency (Hz)', 'Fontname', 'times', 'fontsize', 16);
            ylabel('Power (dB)', 'Fontname', 'times', 'fontsize', 16);
            text(20, 32, '3-7y', 'Fontname', 'times', 'fontsize', 14)
        case 4    
            set(gca, 'Yticklabel', {}, 'Fontname', 'times', 'fontsize', 14)
            xlabel('Frequency (Hz)', 'Fontname', 'times', 'fontsize', 16);
            text(20, 32, '7-16y', 'Fontname', 'times', 'fontsize', 14)
    end
    grid on;
    axis([0 32 -10 40])
end

% END PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot some example EEG recordings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename = 'anon_codex_age_sex.csv';
% Specify range and delimiter

dataLines = [2, Inf];
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = dataLines;
opts.Delimiter = ",";
opts.VariableNames = ["AnonFilename", "QuantizedAgeyears", "QuantizedSex"];
opts.VariableTypes = ["string", "double", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "AnonFilename", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["AnonFilename", "QuantizedSex"], "EmptyFieldRule", "auto");
% Import the data
anoncodex = readtable(filename, opts);
fnames = cellstr(anoncodex{:,1});
ages = anoncodex{:,2};
sex = anoncodex{:,3};

ref1 = find(ages>7 & ages<=20);

% BASIC EVALUATIONS
fs1 = 250;
[B,A] = butter(4, [1 64]./fs1, 'bandpass');

r1 = fs1*5*60+1; r2 = r1+8*fs1-1;
cc(1,:) = [1 0.25 0.25];
cc(2,:) = [1 0.25 0.25];
cc(3,:) = [1 0.25 0.25];
cc(4,:) = [1 0.25 0.25];
cc(5,:) = [0.1 0.6 1];
cc(6,:) = [0.1 0.6 1];
cc(7,:) = [0.1 0.6 1];
cc(8,:) = [0.1 0.6 1];
cc(9,:) = [0.85 0 0];
cc(10,:) = [0.85 0 0];
cc(11,:) = [0.85 0 0];
cc(12,:) = [0.85 0 0];
cc(13,:) = [0 0 0.96];
cc(14,:) = [0 0 0.96];
cc(15,:) = [0 0 0.96];
cc(16,:) = [0 0 0.96];
cc(17,:) = [0 0 0];
cc(18,:) = [0 0 0];
fr = [8 57 694 158];

str{1} = 'Fp2-F4';     % Fp2-F4, Fp2-F4
str{2} = 'F4-C4';     % F4-C4, F4-C4
str{3} = 'C4-P4';     % C4-P4, C4-P4
str{4} = 'P4-O2';     % P4-O2, P4-O2
str{5} = 'Fp1-F3';     % Fp1-F3, Fp1-F3
str{6} = 'F3-C3';     % F3-C3, F3-C3
str{7} = 'C3-P3';     % C3-P3, C3-P3
str{8} = 'P3-O1';     % P3-O1, P3-O1
str{9} = 'Fp2-F8';     % Fp2-F8, Fp2-F8
str{10} = 'F8-T4';     % F8-T4, F8-T8
str{11} = 'T4-T6';     % T4-T6, T8-P8
str{12} = 'T6-O2';     % T6-O2, P8-O2
str{13} = 'Fp1-F7';     % Fp1-F7, Fp1-F7
str{14} = 'F7-T3';     % F7-T3, F7-T7
str{15} = 'T3-T5';     % T3-T5, T7-P7
str{16} = 'T5-O1';     % T5-O1, P7-O1
str{17} = 'Fz-Cz';     % Fz-Cz, Fz-Cz
str{18} = 'Cz-Pz';     % Cz-Pz, Cz-Pz
str{19} = 'ECG';

% START PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1); hold on
[dat, hdr, label, fs, scle, offs] = read_edf(fnames{fr(1)});
[data_bp_mont, ecg] = read_into_montage(dat, label, scle, fs, fs1, 0);
dbp = filtfilt(B,A, data_bp_mont')';
t = linspace(0,8, length(r1:r2));
for jj = 1:18
    plot(t, dbp(jj,r1:r2)-jj*100, 'color', cc(jj,:))
end
plot(t, (ecg(r1:r2)-median(ecg(r1:r2)))./10-(jj+1)*100, 'Color', [0 0.6 0])
set(gca, 'Ytick', [-1900:100:-100], 'Yticklabel', str(19:-1:1), 'Xticklabel', {})
axis([0 8 -2000 0])
set(gca, 'Position', [0.075 0.535 0.435 0.45], 'fontname', 'times', 'fontsize', 10)
text(0.175, -150, 'A', 'BackgroundColor', [1 1 1], 'fontname', 'times', 'fontsize', 14)

subplot(2,2,2); hold on
[dat, hdr, label, fs, scle, offs] = read_edf(fnames{fr(2)});
[data_bp_mont, ecg] = read_into_montage(dat, label, scle, fs, fs1, 0);
dbp = filtfilt(B,A, data_bp_mont')';
t = linspace(0,8, length(r1:r2));
for jj = 1:18
    plot(t, dbp(jj,r1:r2)-jj*100, 'color', cc(jj,:))
end
plot(t, (ecg(r1:r2)-median(ecg(r1:r2)))./20-(jj+1)*100, 'Color', [0 0.6 0])
set(gca, 'Ytick', [-1900:100:-100], 'Yticklabel', {}, 'Xticklabel', {})
axis([0 8 -2000 0])
set(gca, 'Position', [0.55 0.535 0.435 0.45], 'fontname', 'times', 'fontsize', 10)
text(0.175, -150, 'B', 'BackgroundColor', [1 1 1], 'fontname', 'times', 'fontsize', 14)

subplot(2,2,3); hold on
[dat, hdr, label, fs, scle, offs] = read_edf(fnames{fr(3)});
[data_bp_mont, ecg] = read_into_montage(dat, label, scle, fs, fs1, 0);
dbp = filtfilt(B,A, data_bp_mont')';
t = linspace(0,8, length(r1:r2));
for jj = 1:18
    plot(t, dbp(jj,r1:r2)-jj*100, 'color', cc(jj,:))
end
plot(t, (ecg(r1:r2)-median(ecg(r1:r2)))./10-(jj+1)*100, 'Color', [0 0.6 0])
set(gca, 'Ytick', [-1900:100:-100], 'Yticklabel', str(19:-1:1))
axis([0 8 -2000 0])
set(gca, 'Position', [0.075 0.075 0.435 0.45], 'fontname', 'times', 'fontsize', 10)
xlabel('time (s)')
text(0.175, -150, 'C', 'BackgroundColor', [1 1 1], 'fontname', 'times', 'fontsize', 14)

subplot(2,2,4); hold on
[dat, hdr, label, fs, scle, offs] = read_edf(fnames{fr(4)});
[data_bp_mont, ecg] = read_into_montage(dat, label, scle, fs, fs1, 0);
dbp = filtfilt(B,A, data_bp_mont')';
t = linspace(0,8, length(r1:r2));
for jj = 1:18
    plot(t, dbp(jj,r1:r2)-jj*100, 'color', cc(jj,:))
end
plot(t, (ecg(r1:r2)-median(ecg(r1:r2)))./10-(jj+1)*100, 'Color', [0 0.6 0])
set(gca, 'Ytick', [-1900:100:-100], 'Yticklabel', {})
axis([0 8 -2000 0])
set(gca, 'Position', [0.55 0.075 0.435 0.45], 'fontname', 'times', 'fontsize', 10)
xlabel('time (s)')
text(0.175, -150, 'D', 'BackgroundColor', [1 1 1], 'fontname', 'times', 'fontsize', 14)
plot([1.5 1.5 2.5], [-950 -1050 -1050], 'k', 'LineWidth', 1)

% END PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ECG ANALYSIS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename = 'anon_codex_age_sex.csv';
% Specify range and delimiter
dataLines = [2, Inf];
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = dataLines;
opts.Delimiter = ",";
opts.VariableNames = ["AnonFilename", "QuantizedAgeyears", "QuantizedSex"];
opts.VariableTypes = ["string", "double", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "AnonFilename", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["AnonFilename", "QuantizedSex"], "EmptyFieldRule", "auto");
% Import the data
anoncodex = readtable(filename, opts);
fnames = cellstr(anoncodex{:,1});
ages = anoncodex{:,2};
sex = anoncodex{:,3};

% BASIC EVALUATIONS
load('artefact_detection_models.mat');
fs1 = 250; fs2 = 12; epl = 5*60*fs1; olap = 30*fs1;
[Bn0,An0] = butter(4, 2*2/fs1, 'high'); 
[Bn1,An1] = butter(2, 2*[48 52]/fs1 , 'stop');
[Bn2,An2] = butter(2, 2*[98 102]/fs1 , 'stop');

% CALCULATE HRV FEATURES FROM 5-minute epochs of ECG
tadj = fs1*10; % correct for any poential end effects
feat = cell(length(fnames),2);
for ii = 1:length(fnames)
    ii
   
        [dat, hdr, label, fs, scle, offs] = read_edf(fnames{ii});
        [ecg, fs] = read_into_montage_ecg(dat, label, scle, fs);
        ecg = resample(ecg, fs1, fs); 
        ecg = filtfilt(Bn1, An1, ecg); ecg = filtfilt(Bn2, An2, ecg); % notch at 50 notch at 100
        ecg = filtfilt(Bn0,An0, ecg); 
        block_no = floor(length(ecg)/olap)-epl/olap+1; fts = NaN*ones(block_no,38);
        MM = (floor(length(ecg)/fs1/2))*2;
        flag1 = 0; flag2 = 0;
        fts = NaN*ones(block_no, 38); ag = NaN*ones(block_no,1);
        try
        for jj = 1:block_no    
            r1 = (jj-1)*olap+1-tadj; r2 = r1+epl-1+2*tadj;
            if r1<1; r1 = 1; flag1 = 1; end; if r2>MM*fs1; r2 = MM*fs1; flag2 = 1; end
            ecg1 = ecg(r1:r2); ecg2 = ecg((jj-1)*olap+1:(jj-1)*olap+epl); 
            if flag1 == 1; ecg1 = [zeros(1,tadj) ecg1]; end
            if flag2 == 1; ecg1 = [ecg1 zeros(1,tadj)]; end
            qts = quantile(ecg1, [0.01 0.99]); if abs(qts(1))>abs(qts(2)); ecg1 = -ecg1; end
            qts = quantile(ecg2, [0.01 0.99]); if abs(qts(1))>abs(qts(2)); ecg2 = -ecg2; end
            afts = calculate_features_ecg(ecg2, fs1);
            [~, out1] = predict(Mdl1, afts); [~, out2] = predict(Mdl2, afts); [~, out3] = predict(Mdl3, afts); 
            if (out1(2) < 0) && (out2(2) < 0.5) && (out3(2) < 0)
                [~, rr, ~] = pan_tompkin_adapt(ecg1,fs1);
                rr(rr<=tadj)=0; rr(rr>=epl+tadj)=0;
                rr = rr(rr>0);
                if (length(rr)>200) && (length(rr)<900)
                    fts(jj,:) = calculate_features_valid(rr, fs);     
                    ag(jj) = ages(ii);
                end
            end
        end

        feat{ii,1} = fts;
        feat{ii,2} = ag;
        catch 

        end

end

save('ecg_features_online_data.mat', 'feat', '-v7.3');
load('ecg_features_online_data.mat');
% Find outliers
fv1 = NaN*ones(length(feat), 38); fv2 = fv1; age = zeros(length(feat),1);
for ii = 1:length(feat)
    fts = feat{ii,1};
    fv1(ii,:) = median(fts(1:2,:));
    fv2(ii,:) = median(fts(11:end,:));
    age(ii) = feat{ii,2}(1);
end
rf = find(~isnan(age) & ~isnan(sum(fv1'))' & ~isnan(sum(fv2'))');
fv1 = fv1(rf,:);
fv2 = fv2(rf,:);
age = age(rf,:);

rfs = zeros(1,length(fv1));
for jj = 1:38
    y = fv1(:,jj); x = age;
    B = polyfit(x,y,2);
    res = y - polyval(B, x);
    dum = zeros(1, length(y)); dum(find(abs(res)>6*std(res)))=1;
    rfs = rfs+dum;
    y = fv2(:,jj); x = age;
    B = polyfit(x,y,2);
    res = y - polyval(B, x);
    dum = zeros(1, length(y)); dum(find(abs(res)>6*std(res)))=1;
    rfs = rfs+dum;
end

% remove outliers
nref = find(rfs==0);
fv1 = fv1(nref,:);
fv2 = fv2(nref,:);
age = age(nref);

% TRAIN GPR age predictor - 5 fold CV 80:20 split
pid_in1 = 1:length(fv1); 
pma_in1 = age;
Mc = 5;
out = zeros(1,1000);
for qq = 1:1000
    rng(qq)
    pd = unique(pid_in1);
    y = cell(1,Mc);
    rx = ones(1,length(pd)); %rt = 1-rx; %rr = 1:length(pid);
    K = floor(length(pd)/Mc); dum = rem(length(pd), Mc);
    K = K.*ones(1,Mc); K(1:dum) = K(1)+1;
    for ii = 1:Mc
       rz = find(rx==1);
       dum = randsample(length(rz), K(ii), false);
       rx(rz(dum)) = 0;
       y{ii} = rz(dum);
    end
    y{end} = [y{end}' ; find(rx==1)']';
    for ii = 1:Mc
        y{ii} = pd(y{ii});
    end
    ss = cell(1, Mc);
    for ii = 1:Mc
        rf = y{ii};
        ag = [];
        for jj = 1:length(rf); ag = [ag pma_in1(find(pid_in1==rf(jj)))]; end
        ss{ii} = ag;
    end
    pv = zeros(1,6); c1 = 1;
    for ii = 1:Mc
        for jj = ii+1:Mc
            [~, pv(c1)] = kstest2(ss{ii}, ss{jj});
            c1 = c1+1;
        end
    end
    out(qq) = mean(pv);
end
nr = find(out==max(out)); % nr = 7 
rng(nr)
yr = cell(1,Mc);
rx = ones(1,length(pd)); %rt = 1-rx; %rr = 1:length(pid);
for ii = 1:Mc
   rz = find(rx==1);
   dum = randsample(length(rz), K(ii), false);
   rx(rz(dum)) = 0;
   yr{ii} = rz(dum);
end
yr{end} = [yr{end}' ; find(rx==1)']';
for ii = 1:Mc; yr{ii} = pd(yr{ii}); end

pred1 = []; pred2 = []; agx = [];
for ii = 1:Mc
    ii
    dum = zeros(1,length(fv1));
    dum(yr{ii})=1;
    r1 = find(dum==1); r2 = find(dum==0);

    rGP1 = fitrgp(fv1(r2,:), age(r2), 'BasisFunction', 'constant', 'KernelFunction', 'matern52', 'Standardize', true);
    out1 = predict(rGP1, fv1(r1,:));
    pred1 = [pred1 ; out1]; 

    rGP2 = fitrgp(fv2(r2,:), age(r2), 'BasisFunction', 'constant', 'KernelFunction', 'matern52', 'Standardize', true);
    out2 = predict(rGP1, fv2(r1,:));
    pred2 = [pred2 ; out2]; 

    agx = [agx ; age(r1)];

end

% START PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(1,2,1)
plot(agx, pred1 ,'.')
title('Pre- N2 Onset')
set(gca, 'FontName', 'times', 'Fontsize', 16)
xlabel('Age (y)'); ylabel('Predicted Age (y)')
grid on; hold on;
plot([-2 18], [-2, 18], 'k')
axis([-1 17 -1, 17])
set(gca, 'position', [0.075 0.15 0.4 0.8], 'Xtick', [0:2:16], 'Ytick', [0:2:16])
%text(0, 16, 'A', 'FontName', 'times', 'Fontsize', 14)
pp = corr(agx(~isnan(agx)), pred1(~isnan(agx)));
text(10, 2, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 16)

subplot(1,2,2)
plot(agx, pred2 ,'.')
title('Post- N2 Onset')
set(gca, 'FontName', 'times', 'Fontsize', 16)
xlabel('Age (y)'); ylabel('Predicted Age (y)')
grid on; hold on;
plot([-2 18], [-2, 18], 'k')
axis([-1 17 -1, 17])
set(gca, 'position', [0.575 0.15 0.4 0.8], 'Xtick', [0:2:16], 'Ytick', [0:2:16])
%text(0, 16, 'B', 'FontName', 'times', 'Fontsize', 14)
pp = corr(agx(~isnan(agx)), pred2(~isnan(agx)));
text(10, 2, ['r = ' num2str(pp, '%1.3f')], 'FontName', 'times', 'FontSize', 16)

% END PLOTTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% just artefact stuff

filename = 'anon_codex_age_sex.csv';
% Specify range and delimiter
dataLines = [2, Inf];
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = dataLines;
opts.Delimiter = ",";
opts.VariableNames = ["AnonFilename", "QuantizedAgeyears", "QuantizedSex"];
opts.VariableTypes = ["string", "double", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "AnonFilename", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["AnonFilename", "QuantizedSex"], "EmptyFieldRule", "auto");
% Import the data
anoncodex = readtable(filename, opts);
fnames = cellstr(anoncodex{:,1});
ages = anoncodex{:,2};
sex = anoncodex{:,3};

% BASIC EVALUATIONS
load('artefact_detection_models.mat');
fs1 = 250; fs2 = 12; epl = 5*60*fs1; olap = epl;
[Bn0,An0] = butter(4, 2*2/fs1, 'high'); 
[Bnx,Anx] = butter(4, 1/fs1, 'high'); 
[Bn1,An1] = butter(2, 2*[48 52]/fs1 , 'stop');
[Bn2,An2] = butter(2, 2*[98 102]/fs1 , 'stop');
W = floor(fs1/3);

tadj = fs1*10; % correct for any poential end effects
feat = cell(length(fnames),3);
fband = [0.5 1 ; 1 3 ; 3 25 ; 25 125]; df = diff(fband');
for ii = 1:length(fnames)
    ii
   
        [dat, hdr, label, fs, scle, offs] = read_edf(fnames{ii});
        [ecg, fs] = read_into_montage_ecg(dat, label, scle, fs);
        ecg = resample(ecg, fs1, fs); 
        ecg = filtfilt(Bn1, An1, ecg); ecg = filtfilt(Bn2, An2, ecg); % notch at 50 notch at 100
        ecg0 = filtfilt(Bnx, Anx, ecg);
        ecg = filtfilt(Bn0,An0, ecg); 
        block_no = floor(length(ecg)/olap)-epl/olap+1; fts = NaN*ones(block_no,38);
        MM = (floor(length(ecg)/fs1/2))*2;
        flag1 = 0; flag2 = 0;
        fts = NaN*ones(block_no, 38); ag = NaN*ones(block_no,1);
        try
        artx1 = []; artx2 = []; artx3 = []; c1 = 0;
        pow = zeros(block_no,4);
        for jj = 1:block_no    
            r1 = (jj-1)*olap+1-tadj; r2 = r1+epl-1+2*tadj;
            if r1<1; r1 = 1; flag1 = 1; end; if r2>MM*fs1; r2 = MM*fs1; flag2 = 1; end
            ecg1 = ecg(r1:r2); ecg2 = ecg((jj-1)*olap+1:(jj-1)*olap+epl); 
            %spectral_analysis
            [Pxx, f] = pwelch(ecg2, hamming(2^13), 2^12, 2^13, fs1);
            for qq = 1:4
                rf = find(f>=fband(qq,1) & f<fband(qq,2));
                pow(jj,qq) = sum(Pxx(rf))./df(qq);
            end           
            if flag1 == 1; ecg1 = [zeros(1,tadj) ecg1]; end
            if flag2 == 1; ecg1 = [ecg1 zeros(1,tadj)]; end
            qts = quantile(ecg1, [0.01 0.99]); if abs(qts(1))>abs(qts(2)); ecg1 = -ecg1; end
            qts = quantile(ecg2, [0.01 0.99]); if abs(qts(1))>abs(qts(2)); ecg2 = -ecg2; end
            afts = calculate_features_ecg(ecg2, fs1);
            [~, out1] = predict(Mdl1, afts); [~, out2] = predict(Mdl2, afts); [~, out3] = predict(Mdl3, afts); 
          
            if (out1(2) < 0) && (out2(2) < 0.5) && (out3(2) < 0)
                c1 = c1+1;
                [~, rr, ~, art1, art2, art3] = pan_tompkin_adapt_count_error(ecg1,fs1);
                artx1(c1,:) = art1;
                artx2(c1,:) = art2;
                if isempty(art3); artx3(c1,:) = [0 length(rr)]; else; artx3(c1,:) = [art3 length(rr)]; end
            end

        end

        feat{ii,1} = artx2; feat{ii,2} = artx3;
        feat{ii,3} = pow;
        
        catch 

        end

end

save('artefact_features_online_data.mat', 'feat', '-v7.3');

seg_out = []; nn_adj = []; flag_bad = zeros(length(feat),2); pow = [];
for ii = 1:length(feat)
    A1 = size(feat{ii,3}); A2 = size(feat{ii,2});
    if isempty(feat{ii,2})
        flag_bad(ii,:) = [A1(1) 0];
    else
        flag_bad(ii,:) = [A1(1) A2(1)];
    end
    dum = feat{ii,1};
    if ~isempty(dum)
    dum(isnan(dum(:,2)),2) = 0;
    seg_out = [seg_out ; dum];
    nn_adj = [nn_adj ; feat{ii,2}];
    end
    pow = [pow; feat{ii,3}];
end


