% %  Figure 2
load('Figure2a.mat')
figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 7 2.8])
tl=tiledlayout(2,3)
tl.Padding='compact';
tl.TileSpacing='compact';
nexttile([2 1])
plot(posFish(:,1),posFish(:,2),'b','LineWidth',3)
hold on
drawcircle('Center',[18,18],'Radius',18,'StripeColor','green','InteractionsAllowed','none');
xticks([])
yticks([])
quiver(0,0,10,0,1,'Color','black','LineWidth',1,'MaxHeadSize',.5)
quiver(0,0,0,10,1,'Color','black','LineWidth',1,'MaxHeadSize',.5)
box on
axis([0 36 0 36])
pbaspect([1 1 1])
% title('Fish in right tank')
box off
set(gca,'FontName','Arial','FontSize',14,'TickDir','out')

nexttile([2 1])
plot(posRobot(:,1),posRobot(:,2),'r','LineWidth',3)
hold on
drawcircle('Center',[18,18],'Radius',18,'StripeColor','green','InteractionsAllowed','none');
xticks([])
yticks([])
box on
axis([0 36 0 36])
pbaspect([1 1 1])
% title('Robot in left tank')
box off
set(gca,'FontName','Arial','FontSize',11,'TickDir','out')

% nexttile
% plot(posFish(:,1),posFish(:,2),'b','LineWidth',3)
% hold on
% plot(posRobot(:,1),posRobot(:,2),'r','LineWidth',3)
% drawcircle('Center',[18,18],'Radius',18,'StripeColor','green','InteractionsAllowed','none');
% xticks([])
% yticks([])
% box on
% axis([0 36 0 36])
% pbaspect([1 1 1])
% title('Robot and fish trajectories superimposed')
% box off
% set(gca,'FontName','Arial','FontSize',14,'TickDir','out')


% X-Y trajectories
% figure
% tl=tiledlayout(2,3)
% tl.Padding='compact';
% tl.TileSpacing='compact';
% nexttile
% plot(0.05:0.05:20,posFish(:,1),'b','LineWidth',3)
% box on
% axis([0 20 0 36])
% ylabel('X position (cm)')
% title('Fish in right tank')
% xticks([])
% box off
% set(gca,'FontName','Arial','FontSize',14,'TickDir','out')
% 
% nexttile
% plot(0.05:0.05:20,posRobot(:,1),'r','LineWidth',3)
% box on
% axis([0 20 0 36])
% title('Robot in left tank')
% xticks([])
% yticks([])
% box off
% set(gca,'FontName','Arial','FontSize',14,'TickDir','out')


nexttile
plot(0.05:0.05:20,posFish(1:end,1),'b','LineWidth',3)
hold on
plot(0.05:0.05:20,posRobot(1:end,1),'--r','LineWidth',3)
box on
axis([0 20 0 36])
% title('Robot and fish trajectories superimposed')
ylabel('X position (cm)')
xticks([])
box off
set(gca,'FontName','Arial','FontSize',11,'TickDir','out')

% nexttile
% plot(0.05:0.05:20,posFish(:,2),'b','LineWidth',3)
% box on
% axis([0 20 0 36])
% ylabel('Y position (cm)')
% xlabel('Time (s)')
% box off
% set(gca,'FontName','Arial','FontSize',14,'TickDir','out')
% 
% 
% nexttile
% plot(0.05:0.05:20,posRobot(:,2),'r','LineWidth',3)
% box on
% axis([0 20 0 36])
% xlabel('Time (s)')
% yticks([])
% box off
% set(gca,'FontName','Arial','FontSize',14,'TickDir','out')

nexttile
plot(0.05:0.05:20,posFish(1:end,2),'b','LineWidth',3)
hold on
plot(0.05:0.05:20,posRobot(1:end,2),'--r','LineWidth',3)
box on
axis([0 20 0 36])
ylabel('Y position (cm)')
xlabel('Time (s)')
box off
set(gca,'FontName','Arial','FontSize',11,'TickDir','out')


