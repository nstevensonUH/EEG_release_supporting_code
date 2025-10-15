function [c, trr] = rr_count(trr, test_ecg, lim1)

N1 = 1; N2 = 2; c1 =1; N = length(test_ecg);
while N1~=N2
    for ii = 1:length(trr)
        dum = trr(ii);
        t1 = dum-lim1; t2 = dum+lim1;
        if t1<1; t1 = 1; end; if t2>N; t2 = N; end
        dat = test_ecg(t1:t2);
        qq = find(dat == max(dat),1);
        trr(ii) = t1+qq-1;
    end
    N1 = length(trr); N2 = length(unique(trr));
    trr = unique(trr);c1 = c1+1;
end
c = length(trr);