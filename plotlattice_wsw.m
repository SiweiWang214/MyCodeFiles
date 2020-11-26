function plotlattice_wsw(THERING,points)
% 'points' can be 1:length(THERING) 
% Written by S.W.Wang on Nov.15,2019

if nargin < 2
    points = 1:length(THERING);
end
TD = twissring(THERING,0,1:length(THERING)+1,'chrom');
beta = cat(1,TD.beta);
s = cat(1,TD.SPos);
eta = cat(2,TD.Dispersion);
L = s(points(end));

FontNa='Helvetica';
FontSi=18;

figure;
hold on;
box on;
drawlattice_wsw(2,2,points(end));
h1 = plot(s(points),beta(points,1),'linewidth',1.5,'color','k','linestyle','-');
h2 = plot(s(points),beta(points,2),'linewidth',1.5,'color','r','linestyle','-.');
ylabel('\beta functions (m)','FontName',FontNa,'fontsize',FontSi);
yyaxis right
h3 = plot(s(points),eta(1,points),'linewidth',1.5,'color','b','linestyle','--');
l1=legend([h1,h2,h3],'\beta_x','\beta_y','\eta_x');
set(l1,'FontName',FontNa,'FontSize',FontSi,'location','best');
xlabel('s (m)');
ylabel('Horizontal \eta function (m)','FontName',FontNa,'fontsize',FontSi,'color','r');
set(gcf,'position',[400 150 800 600]);
set(gca,'FontName',FontNa,'FontSize',FontSi);
xlim([0 L]);