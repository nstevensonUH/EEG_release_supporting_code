function [rsq]=linplot(x,y,col)
%linear plot line
% idx = isnan(x);
% p = polyfit(x(~idx),y(~idx),1);
p=polyfit(x,y,1); 
f=polyval(p,x); 
% hold on; plot(x,f,'Color',col,'linewidth',2); 

yresid=y-f;

SSresid =sum(yresid.^2);
SStotal = (length(y)-1)*var(y);

plot([min(x) max(x)], [min(f) max(f)]-2, 'k--'); 
plot([min(x) max(x)], [min(f) max(f)]+2, 'k--');

hold on, 
plot([min(x) max(x)], [min(f) max(f)]-1, '--','Color',[0.7 0.7 0.7]); 
plot([min(x) max(x)], [min(f) max(f)]+1, '--','Color',[0.7 0.7 0.7]); 


rsq= 1 - SSresid/SStotal; 
end

