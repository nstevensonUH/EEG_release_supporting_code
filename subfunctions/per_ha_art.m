function anno = per_ha_art(data, fs1);

anno = zeros(size(data));
offs = 4*fs1; N = length(anno);
for ii = 1:length(anno)/fs1
    r1 = (ii-1)*fs1+1-offs; r2 = r1+fs1-1+offs;
    if r1<1; r1=1; end; if r2>N; r2=N; end
    anno(:, r1:r2) = repmat(max(abs(data(:,r1:r2)),[],2), 1, length(r1:r2));
end
anno(anno<500) = 0; anno(anno>=500)=1;