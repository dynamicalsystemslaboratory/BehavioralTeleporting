%% Heatmaps
clearvars,clc;
load('Proxy_1226_cleaned.mat')
spt = 0;
figure('Renderer','Painters','Units','Inches','Position',[2 2 10 6])
for smallCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
    spt = ProxyTrajectories(smallCt).PairID;
    if isequal(TrialFishIds.Small(smallCt),'1')
        small = ProxyTrajectories(smallCt).fish1Pos;
        large = ProxyTrajectories(smallCt).fish2Pos;
    elseif isequal(TrialFishIds.Small(smallCt),'2')
        small = ProxyTrajectories(smallCt).fish2Pos;
        large = ProxyTrajectories(smallCt).fish1Pos;
    end
    subplot(6,12,spt)
    hist3(small,'CDataMode','auto','CDataMapping','direct','Edges',{0:3:33 0:3:33})
    xticks([]);
    yticks([]);
    box on
%     colorbar('Ticks',[0:60:300],'TickLabels',[0:3:15])
    axis square
    view(2)
    if spt==1
        ylabel('Small')
    end
    set(gca,'FontName','Arial','FontSize',11);
end
for smallCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
    spt = ProxyTrajectories(smallCt).PairID+12;
    if isequal(TrialFishIds.Small(smallCt),'1')
        small = ProxyTrajectories(smallCt).fish1Pos;
        large = ProxyTrajectories(smallCt).fish2Pos;
    elseif isequal(TrialFishIds.Small(smallCt),'2')
        small = ProxyTrajectories(smallCt).fish2Pos;
        large = ProxyTrajectories(smallCt).fish1Pos;
    end
    subplot(6,12,spt)
    hist3(large,'CDataMode','auto','CDataMapping','direct','Edges',{0:3:33 0:3:33})
    xticks([]);
    yticks([]);
    box on
    axis square
    view(2)
    if spt==13
        ylabel('Large')
    end
    set(gca,'FontName','Arial','FontSize',11);
end
for smallCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
    spt = ProxyTrajectories(smallCt).PairID+24;
    if isequal(TrialFishIds.Small(smallCt),'L')
        small = ProxyTrajectories(smallCt).fishLPos;
        large = ProxyTrajectories(smallCt).fishRPos;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        small = ProxyTrajectories(smallCt).fishRPos;
        large = ProxyTrajectories(smallCt).fishLPos;
    end
    subplot(6,12,spt)
    hist3(small,'CDataMode','auto','CDataMapping','direct','Edges',{0:3:33 0:3:33})
    xticks([]);
    yticks([]);
    box on
    axis square
    view(2)
    if spt==25
        ylabel('Small')
    end
    set(gca,'FontName','Arial','FontSize',11);
end
for smallCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
    spt = ProxyTrajectories(smallCt).PairID+36;
    if isequal(TrialFishIds.Small(smallCt),'L')
        small = ProxyTrajectories(smallCt).fishLPos;
        large = ProxyTrajectories(smallCt).fishRPos;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        small = ProxyTrajectories(smallCt).fishRPos;
        large = ProxyTrajectories(smallCt).fishLPos;
    end
    subplot(6,12,spt)
    hist3(large,'CDataMode','auto','CDataMapping','direct','Edges',{0:3:33 0:3:33})
    xticks([]);
    yticks([]);
    box on
    axis square
    view(2)
    if spt==37
        ylabel('Large')
    end
    set(gca,'FontName','Arial','FontSize',11);
end
for smallCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
    spt = ProxyTrajectories(smallCt).PairID+48;
    if isequal(TrialFishIds.Small(smallCt),'L')
        small = ProxyTrajectories(smallCt).fishLPos;
        large = ProxyTrajectories(smallCt).fishRPos;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        small = ProxyTrajectories(smallCt).fishRPos;
        large = ProxyTrajectories(smallCt).fishLPos;
    end
    subplot(6,12,spt)
    hist3(small,'CDataMode','auto','CDataMapping','direct','Edges',{0:3:33 0:3:33})
    xticks([]);
    yticks([]);
    box on
    axis square
    view(2)
    if spt==49
        ylabel('Small')
    end
    set(gca,'FontName','Arial','FontSize',11);
end
for smallCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
    spt = ProxyTrajectories(smallCt).PairID+60;
    if isequal(TrialFishIds.Small(smallCt),'L')
        small = ProxyTrajectories(smallCt).fishLPos;
        large = ProxyTrajectories(smallCt).fishRPos;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        small = ProxyTrajectories(smallCt).fishRPos;
        large = ProxyTrajectories(smallCt).fishLPos;
    end
    subplot(6,12,spt)
    hist3(large,'CDataMode','auto','CDataMapping','direct','Edges',{0:3:33 0:3:33})
    xticks([]);
    yticks([]);
    box on
    axis square
    view(2)
    if spt==61
        ylabel('Large')
    end
    xlabel(num2str(spt-60));
    set(gca,'FontName','Arial','FontSize',11);
end
subplot(6,12,31)
xticks([]);
yticks([]);
box on
axis square
subplot(6,12,43)
xticks([]);
yticks([]);
box on
axis square
subplot(6,12,58)
xticks([]);
yticks([]);
box on
axis square
subplot(6,12,70)
xticks([]);
yticks([]);
box on
axis square
xlabel('10')
set(gca,'FontName','Arial','FontSize',11);
