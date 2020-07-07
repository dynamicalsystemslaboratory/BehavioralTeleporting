%% Figure 4
clearvars,clc;
load('Data.mat')

speed=fillmissing(downsample(ProxyTrajectories(1).fish1Speed(1000:1040,1),2),'linear');
turnrate=fillmissing(downsample(ProxyTrajectories(1).fish1TurnRate(1000:1040,1),2),'linear');

symbX1 = diff(speed,1,1)>0;
symbX2 = diff(turnrate,1,1)>0;
symbX = zeros(numel(symbX1),1);
symbX(intersect(find(symbX1==0),find(symbX2==0))) = 1;
symbX(intersect(find(symbX1==0),find(symbX2==1))) = 2;
symbX(intersect(find(symbX1==1),find(symbX2==0))) = 3;
symbX(intersect(find(symbX1==1),find(symbX2==1))) = 4;

ff=figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 7 4.6])
tt=tiledlayout(3,1)
tt.TileSpacing='compact';
nexttile
for ct=1:20
    hold on
    if symbX1(ct)==0
        plot((ct-1)*0.1:0.1:ct*0.1,speed(ct:ct+1),'--k','LineWidth',1)
    elseif symbX1(ct)==1
        plot((ct-1)*0.1:0.1:ct*0.1,speed(ct:ct+1),'-k','LineWidth',1)
    end
    scatter(0:0.1:2,speed,36,'k')
    hold off
    xticklabels([])
    ylabel('Speed (cm/s)')
    axis([0 2 0 12])
    set(gca,'FontName','Ariel','FontSize',10,'TickDir','out');
end
nexttile
for ct=1:20
    hold on
    if symbX2(ct)==0
        plot((ct-1)*0.1:0.1:ct*0.1,turnrate(ct:ct+1),'--k','LineWidth',1)
    elseif symbX2(ct)==1
        plot((ct-1)*0.1:0.1:ct*0.1,turnrate(ct:ct+1),'-k','LineWidth',1)
    end
    scatter(0:0.1:2,turnrate,36,'k')
    hold off
    xticklabels([])
    ylabel('Turn rate (rad/s)')
    axis([0 2 -20 20])
    set(gca,'FontName','Ariel','FontSize',10,'TickDir','out');
end
nexttile
for ct=1:20
    hold on
    if symbX(ct)==1
        scatter(ct*0.1,symbX(ct),100,'ks','filled')
    elseif symbX(ct)==2
        scatter(ct*0.1,symbX(ct),100,'ks','filled')
    elseif symbX(ct)==3
        scatter(ct*0.1,symbX(ct),100,'ks','filled')
    elseif symbX(ct)==4
        scatter(ct*0.1,symbX(ct),100,'ks','filled')
    end
    hold off
    xlabel('Time (s)')
    ylabel('Symbolized time-series')
    axis([0 2 0 5])
    yticks([1 2 3 4])
    yticklabels({'--','-+','+-','++'})
    set(gca,'FontName','Ariel','FontSize',10,'TickDir','out');
end