function dum1 = write_edf(filename, dat, hdr, ns, fs)
% [data, hdr, label] = write_edf(filename);
%
% This functions writes an EDF file as per the format outlined in
%  http://www.edfplus.info/specs/edf.html. Note this version uses a
%  pre-existing header so the function is limited.
%
% INPUT: filename - EDF file name
%                dat - a cell array containing the data in the file (int16 format)
%                 hdr - a cell array the header file information in ASCII format
%
% Nathan Stevenson


fid = fopen(filename, 'w');

% WRITE HEADER (see next commented section as to what each bits relates to with respect to the EDF specification)
c1 = fwrite(fid, hdr, 'char');
% c2 = fwrite(fid, hdr{2}, 'char');
% c3 = fwrite(fid, hdr{3}, 'char');
% c4 = fwrite(fid, hdr{4}, 'char');
% c5 = fwrite(fid, hdr{5}, 'char');
% c6 = fwrite(fid, hdr{6}, 'char');
% c7 = fwrite(fid, hdr{7}, 'char');
% c8 = fwrite(fid, hdr{8}, 'char');
% c9 = fwrite(fid, hdr{9}, 'char');
% c10 = fwrite(fid, hdr{10}, 'char');
% c11 = fwrite(fid, hdr{11}, 'char');

% CORRESPONDING COMPONNET OF HEADER IN EDF FORMAT
% hdr{1} = fread(fid, 256, 'char');         % CONTAINS PATIENT INFORMATION, RECORDING INFORMATION
% ns = char(hdr{1}(253:256))';              % NUMBER OF SIGNALS
% hdr{2} = fread(fid, ns*16, 'char');    % LABEL channel label, temp or HR
% hdr{3} = fread(fid, ns*80,'char');     % TRANSDUCER TYPE
% hdr{4} = fread(fid, ns*8,'char');       % PHYSICAL DIMENSION, voltage - temperature
% hdr{5} = fread(fid, ns*8,'char');       % PHYSICAL MIN
% hdr{6} = fread(fid, ns*8,'char');       % PHYSICAL MAX
% hdr{7} = fread(fid, ns*8,'char');       % DIGITAL MIN
% hdr{8} = fread(fid, ns*8,'char');       % DIGITAL MAX
% hdr{9} = fread(fid, ns*80,'char');     % PRE FILTERING
% hdr{10} = fread(fid, ns*8, 'char');    % SAMPLING NO rec
% hdr{11} = fread(fid, ns*32,'char');     % RESERVED    

len_s = str2num(char(hdr(235:244))');        % START DATE AND TIME and a RESERVED
for ii = 1:len_s
    r1 = fs*(ii-1)+1; r2 = ii*fs;
    dum = dat(:,r1:r2);
    d1 = reshape(dum', 1, fs*ns);
    fwrite(fid, d1, 'short');    
end
dum1 = fclose(fid);

