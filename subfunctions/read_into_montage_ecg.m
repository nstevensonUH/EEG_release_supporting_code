function [data_ref_mont, fs] = read_into_montage_ecg(dat, label, scle, fs);
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

str1{1} = 'EKG'; str1{2} = 'ECG'; 
for ii = 1:length(label)
    label{ii} = label{ii}';
end
val = contains(label, 'EKG') | contains(label, 'ECG');  
rx = find(val==1, 1);
data_ref_mont = double(dat{rx})*scle(rx);
fs = fs(rx);

