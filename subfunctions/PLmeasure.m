function k=PLmeasure(hr,T);
fs=256;
t1 = linspace(T(1),T(end),(T(end)-T(1))*fs);
HR=pchip(T,hr,t1);%polynomial fit curve of heart rate


hilbt_HR=abs(hilbert(HR-mean(HR)));%apply hilbert transform
%find threshold value which will return the most bursts
qi=quantile(hilbt_HR,200);
brst_num=zeros(1,200);
for k=1:200;
ru=find([0,hilbt_HR]<qi(k)&[hilbt_HR,0]>=qi(k));
brst_num(k)=length(ru);
end
[BN,n]=max(brst_num);
rd=find([0,hilbt_HR]>qi(n)& [hilbt_HR,0]<=qi(n));
ru=find([0,hilbt_HR]<qi(n)& [hilbt_HR,0]>=qi(n));

%Calculate powerlaw parameters
BA=zeros(1,length(rd));
BS=BA; BH=BA;BK=BA;BW=BA;
for i=1:length(rd)-1;
HRburst=hilbt_HR(ru(i):rd(i));
BA(i)=trapz(HRburst); %burst area
BS(i)=skewness(HRburst); %burst skewness
BH(i)=unique(max(HRburst)); %burst peak (amplitude)
BW(i)=rd(i)-ru(i);
BK(i)=kurtosis(HRburst);%burst kurtosis
end
%filtering out oscillations from hilbert transform
% Y=30;
% md=mod(length(BA),Y);
% o=BA;
% lo=length(o);
% o(lo-md+1:lo)=[];
% rBA=reshape(o,[],Y);
% A=mean(min(rBA.'));
% fO=find(BA<=(A+0.05*A));
% BA(fO)=[];BS(fO)=[];BH(fO)=[];BK(fO)=[];BW(fO)=[]; %delete values identified as hilbert oscillations
[xmin,alpha]=PLxmin(BA);% find the parameters for the powerlaw distribution

FO=find(BA>=xmin);
BW=BW(FO);BH=BH(FO);BS=BS(FO);BK=BK(FO);

BSK=fitlm(BS,BK);
BSKc=BSK.Coefficients.Estimate(1);
BSKm=BSK.Coefficients.Estimate(2);
BWH=fitlm(BW,BH);
BWHc=BWH.Coefficients.Estimate(1);
BWHm=BWH.Coefficients.Estimate(2);
BAavg=mean(BA(FO));
k=[xmin,alpha,BSKc,BSKm,BAavg,BWHc,BWHm,BN,BN*BAavg];
end

function [xmin,alpha]=PLxmin(x)
z= sort(x);
xmins = unique(z); % search over all unique values of x
dat = zeros(size(xmins));

for xm=1:length(xmins)
    xmin = xmins(xm); % choose next xmin candidate
    z_h = z(z>=xmin); % truncate data below this xmin value
    N = length(z); %
    a = 1 + N ./ sum( log(z_h./xmin) ); % estimate alpha using direct MLE
    [fe,n]=ecdf(z_h); % construct the empirical CDF
    ft = (n./xmin).^(1-a); % construct the fitted theoretical CDF
    dat(xm) = unique(max(abs(ft-fe))); % compute the KS statistic
end
D = min(dat); % find smallest D value
xmin = xmins(find(dat<=D,1,'first')); % find corresponding xmin value
z = x(x>=xmin); %
N = length(z); %
alpha = 1 + N ./ sum( log(z./xmin) ); % get corresponding alpha estimate
L = N*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin)); % log-likelihood at estimate
end