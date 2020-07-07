clearvars,clc;
figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 7 2.8])
tiledlayout(1,2)
nexttile
% load('TE_DS_Entropy_02242020.mat');
load('TE_DS_Entropy_All.mat');
% bar(0.1:0.05:1,mean([EntSmall;EntLarge],1),'FaceColor',[150 150 150]./255,'EdgeColor','k','LineWidth',1)
axis([0.05 1.05 1.99 1.9971])
hold on
for i=2:20
errorbar(i*0.05,mean([EntSmall(:,i-1);EntLarge(:,i-1)],1),std([EntSmall(:,i-1);EntLarge(:,i-1)],1)/sqrt(68),'k.','LineWidth',1,'CapSize',0,'MarkerSize',15)
end
hold off
ylabel('Entropy (bits)')
xlabel('Resolution (s)')
set(gca,'FontName','Arial','FontSize',11,'LineWidth',1,'box','off','tickdir','out')

clearvars,clc;
load('Match_TE_Delay.mat')
nexttile
bar(0:4,mean(teSmalltoLarge_Control(13:23,:)+teLargetoSmall_Control(13:23,:),1),'FaceColor',[150 150 150]./255,'EdgeColor','k','LineWidth',1)
axis([-1 5 0.01 0.013])
hold on
for i=0:4
errorbar(i,mean(teSmalltoLarge_Control(13:23,i+1)+teLargetoSmall_Control(13:23,i+1),1),...
    [],std(teSmalltoLarge_Control(13:23,i+1)+teLargetoSmall_Control(13:23,i+1),1)/sqrt(11),'k','LineWidth',1,'CapSize',0)
end
hold off
ylabel('Sum of transfer entropy (bits)')
xlabel('Time-delay')
set(gca,'FontName','Arial','FontSize',11,'LineWidth',1,'box','off','tickdir','out')