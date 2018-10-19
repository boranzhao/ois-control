function goodplot(papersize, margin, fontsize)
% function which produces a nice-looking plot
% and sets up the page for nice printing
if nargin == 0
    papersize = [7 4];
    margin = [0.01 0.01];
    fontsize = 18;
elseif nargin == 1
    margin = [0.01 0.01];
    fontsize = 18;
elseif nargin == 2
    fontsize = 18;
end
set(get(gca,'xlabel'),'FontName','Arial','FontSize', fontsize);
set(get(gca,'ylabel'),'FontName','Arial','FontSize', fontsize);
set(get(gca,'title'),'FontName','Arial','FontSize', fontsize-3);
% box off; axis square;
% set(gca,'LineWidth',2);
% set(gca,'DefaultLineLineWidth',2)
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
% set(gca,'FontWeight','Bold');
% set(gcf,'color','w');

h = findobj(gcf,'Tag','legend');
set(h,'FontSize',13) % 15 is default

set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', papersize);
set(gcf,'PaperPosition',[margin(1) margin(2) papersize(1)+0.2-2*margin(1) papersize(2)-2*margin(2)]);
set(gcf,'PaperPositionMode','Manual');
% for print 
%print -painters -dpdf -r150 ModelSetFR.pdf
end