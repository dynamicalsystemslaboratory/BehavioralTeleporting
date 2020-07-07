%%iScience Supplementary Figure for time analysis
clearvars,clc;
load('Supplement_2minAnalysis.mat');
figure('Renderer','Painters','Units','Inches','Position',[2 2 6 4]);
tt = tiledlayout(2,1);
tt.Padding = 'Compact';
% Polarization
nexttile
errorbar(1:5,pol_means(1:5),pol_std(1:5).*sqrt(12),'ko--',...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','white','CapSize',0);
hold on
errorbar(.95:1:4.95,pol_means(6:10),pol_std(6:10).*sqrt(11),'k^-',...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','blue','CapSize',0);
errorbar(1.05:1:5.05,pol_means(11:15),pol_std(11:15).*sqrt(11),'ks-',...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','red','CapSize',0);
hold off
xticks(1:5);
xticklabels(2:2:10);
yticks(0.5:0.1:1)
xlabel('Time (min)');
ylabel('Polarization');
box off
set(gca,'FontName','Arial','FontSize',11,'tickdir','out','LineWidth',1);
axis([0.5 5.5 0.5 1])
legend('Control','Match','Mismatch')
% Distance Between
nexttile
errorbar(1:5,dist_means(1:5),dist_std(1:5).*sqrt(12),'ko--',...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','white','CapSize',0);
hold on
errorbar(.95:1:4.95,dist_means(6:10),dist_std(6:10).*sqrt(11),'k^-',...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','blue','CapSize',0);
errorbar(1.05:1:5.05,dist_means(11:15),dist_std(11:15).*sqrt(11),'ks-',...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','red','CapSize',0);
hold off
xticks(1:5);
yticks(0:5:25)
xticklabels(2:2:10);
% xlabel('Time (min)');
ylabel('Distance between\newlinetwo fish (cm)');
box off
set(gca,'FontName','Arial','FontSize',11,'tickdir','out','LineWidth',1);
axis([0.5 5.5 0 25])

% % Angular momentum
% nexttile
% errorbar(1:5,ang_means(1:5),ang_std(1:5).*sqrt(12),'ko--',...
%     'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','white','CapSize',0);
% hold on
% errorbar(.95:1:4.95,ang_means(6:10),ang_std(6:10).*sqrt(11),'k^-',...
%     'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','blue','CapSize',0);
% errorbar(1.05:1:5.05,ang_means(11:15),ang_std(11:15).*sqrt(11),'ks-',...
%     'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','red','CapSize',0);
% hold off
% xticks(1:5);
% yticks(0:0.2:1)
% xticklabels(2:2:10);
% xlabel('Time (min)');
% ylabel('Angular momentum');
% box off
% set(gca,'FontName','Arial','FontSize',11,'tickdir','out','LineWidth',1);
% axis([0.5 5.5 0 1])