function plotshaded(x,y,fstr,levels,outliers)
% x: x coordinates
% y: either just one y vector, or 2xN or 3xN matrix of y-data
% fstr: format ('r' or 'b--' etc)
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');

if nargin < 5
    if exist('outliers','var') == 0
        outliers = 0;
    else
        error('plotshaded:argChk','Wrong number of input arguments');
    end
end

if nargin < 4
    if exist('levels','var') == 0
        levels = 1;
    else
        error('plotshaded:argChk','Wrong number of input arguments');
    end
end

elements = size(y,1);

Datmean = mean(y);
Datminimum = min(y);
Datmaximum = max(y);

Datmin = min(y);
Datmax = max(y);
px=[x,fliplr(x)];
py=[Datmin, fliplr(Datmax)];
patch(px,py,1,'FaceAlpha',0.7*(1/(levels+1)),'FaceColor',fstr,'EdgeColor','none');
y = sort(y);
y([1:floor(elements*outliers/2),size(y,1)-floor(elements*outliers/2)+1:size(y,1)],:) = [];

for i = 1:levels
    Datmin = min(y);
    Datmax = max(y);
    px=[x,fliplr(x)];
    py=[Datmin, fliplr(Datmax)];
    patch(px,py,1,'FaceAlpha',0.7*(1/(levels+1)),'FaceColor',fstr,'EdgeColor','none');
    y = sort(y);
    y([1:floor(elements/levels/2),size(y,1)-floor(elements/levels/2)+1:size(y,1)],:) = [];
    hold on
end

plot(x,Datmean,fstr,'LineWidth',1.5);
spread = max(Datmaximum) - min(Datminimum);
axis([min(x), max(x), min(Datminimum) - 0.1 * spread, max(Datmaximum) + 0.1 * spread])
end