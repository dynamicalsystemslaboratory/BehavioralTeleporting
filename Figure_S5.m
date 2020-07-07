%% Polarization for supplemental mat
clearvars,clc;
load('Proxy_1226_cleaned.mat');
ctS = 0;
pol=struct;
DistBet=struct;
pol.Control=[];
DistBet.Control=[];
AngMoment.Control=[];

pol.Match=[];
DistBet.Match=[];
AngMoment.Match=[];

pol.Mismatch=[];
DistBet.Mismatch=[];
AngMoment.Mismatch=[];

pol.Chance=[];
DistBet.Chance=[];
AngMoment.Chance=[];

for smallCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
    ctS = ctS+1
    if isequal(TrialFishIds.Small(smallCt),'1')
        small = ProxyTrajectories(smallCt).fish1Pos;
        large = ProxyTrajectories(smallCt).fish2Pos;
    elseif isequal(TrialFishIds.Small(smallCt),'2')
        small = ProxyTrajectories(smallCt).fish2Pos;
        large = ProxyTrajectories(smallCt).fish1Pos;
    end
%     temp=[];
    for tt=1:30
        wndwS = (tt-1)*400 + 1;
        wndwe = (tt)*400;
        pol.Control= [pol.Control;nanmean(vecnorm((diff(small(wndwS:wndwe,:),1,1)./vecnorm(diff(small(wndwS:wndwe,:),1,1),2,2))+...
            (diff(large(wndwS:wndwe,:),1,1)./vecnorm(diff(large(wndwS:wndwe,:),1,1),2,2)),2,2)./2)];
        DistBet.Control = [DistBet.Control;nanmean(sqrt((small(wndwS:wndwe,1)-large(wndwS:wndwe,1)).^2+(small(wndwS:wndwe,2)-large(wndwS:wndwe,2)).^2))];
        sp1 = [0,0;diff(small(wndwS:wndwe,:),1,1)].*20;
        sp2 = [0,0;diff(large(wndwS:wndwe,:),1,1)].*20;
        temp = vecnorm(cross([(small(wndwS:wndwe,:)-large(wndwS:wndwe,:)).*0.5,zeros(400,1)],[sp1,zeros(400,1)])+...
            cross([(-small(wndwS:wndwe,:)+large(wndwS:wndwe,:)).*0.5,zeros(400,1)],[sp2,zeros(400,1)]),2,2)./...
            (vecnorm([(small(wndwS:wndwe,:)-large(wndwS:wndwe,:)).*0.5,ones(400,1)],2,2).*...
            vecnorm([(-small(wndwS:wndwe,:)+large(wndwS:wndwe,:)).*0.5,ones(400,1)],2,2).*...
            vecnorm([sp1,zeros(400,1)],2,2).*vecnorm([sp2,zeros(400,1)],2,2));
        temp(temp>1)=NaN;
        AngMoment.Control = [AngMoment.Control;nanmean(temp,'all')];
    end  
%     pol(end+1)=nanmean(temp,'all');
    
end

ctS=0;
for smallCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
    ctS = ctS+1
    if isequal(TrialFishIds.Small(smallCt),'L')
        small = ProxyTrajectories(smallCt).fishLPos;
        large = ProxyTrajectories(smallCt).fishRPos;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        small = ProxyTrajectories(smallCt).fishRPos;
        large = ProxyTrajectories(smallCt).fishLPos;
    end
%     temp=[];
    for tt=find(delayAnl(smallCt).Simultaneous==1)'
        wndwS = (tt-1)*400 + 1;
        wndwe = (tt)*400;
        pol.Match = [pol.Match;nanmean(vecnorm((diff(small(wndwS:wndwe,:),1,1)./vecnorm(diff(small(wndwS:wndwe,:),1,1),2,2))+...
            (diff(large(wndwS:wndwe,:),1,1)./vecnorm(diff(large(wndwS:wndwe,:),1,1),2,2)),2,2)./2)];
        DistBet.Match = [DistBet.Match;nanmean(sqrt((small(wndwS:wndwe,1)-large(wndwS:wndwe,1)).^2+(small(wndwS:wndwe,2)-large(wndwS:wndwe,2)).^2))];       
        sp1 = [0,0;diff(small(wndwS:wndwe,:),1,1)].*20;
        sp2 = [0,0;diff(large(wndwS:wndwe,:),1,1)].*20;
        temp = vecnorm(cross([(small(wndwS:wndwe,:)-large(wndwS:wndwe,:)).*0.5,zeros(400,1)],[sp1,zeros(400,1)])+...
            cross([(-small(wndwS:wndwe,:)+large(wndwS:wndwe,:)).*0.5,zeros(400,1)],[sp2,zeros(400,1)]),2,2)./...
            (vecnorm([(small(wndwS:wndwe,:)-large(wndwS:wndwe,:)).*0.5,ones(400,1)],2,2).*...
            vecnorm([(-small(wndwS:wndwe,:)+large(wndwS:wndwe,:)).*0.5,ones(400,1)],2,2).*...
            vecnorm([sp1,zeros(400,1)],2,2).*vecnorm([sp2,zeros(400,1)],2,2));
        temp(temp>1)=NaN;
        AngMoment.Match = [AngMoment.Match;nanmean(temp,'all')];
    end   
%     pol(end+1)=nanmean(temp,'all');
    
end

ctS=0;
for smallCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
    ctS = ctS+1
    if isequal(TrialFishIds.Small(smallCt),'L')
        small = ProxyTrajectories(smallCt).fishLPos;
        large = ProxyTrajectories(smallCt).fishRPos;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        small = ProxyTrajectories(smallCt).fishRPos;
        large = ProxyTrajectories(smallCt).fishLPos;
    end
%     temp=[];
    for tt=find(delayAnl(smallCt).Simultaneous==1)'
        wndwS = (tt-1)*400 + 1;
        wndwe = (tt)*400;
        pol.Mismatch = [pol.Mismatch;nanmean(vecnorm((diff(small(wndwS:wndwe,:),1,1)./vecnorm(diff(small(wndwS:wndwe,:),1,1),2,2))+...
            (diff(large(wndwS:wndwe,:),1,1)./vecnorm(diff(large(wndwS:wndwe,:),1,1),2,2)),2,2)./2)];
        DistBet.Mismatch = [DistBet.Mismatch;nanmean(sqrt((small(wndwS:wndwe,1)-large(wndwS:wndwe,1)).^2+(small(wndwS:wndwe,2)-large(wndwS:wndwe,2)).^2))];     
        sp1 = [0,0;diff(small(wndwS:wndwe,:),1,1)].*20;
        sp2 = [0,0;diff(large(wndwS:wndwe,:),1,1)].*20;
        temp = vecnorm(cross([(small(wndwS:wndwe,:)-large(wndwS:wndwe,:)).*0.5,zeros(400,1)],[sp1,zeros(400,1)])+...
            cross([(-small(wndwS:wndwe,:)+large(wndwS:wndwe,:)).*0.5,zeros(400,1)],[sp2,zeros(400,1)]),2,2)./...
            (vecnorm([(small(wndwS:wndwe,:)-large(wndwS:wndwe,:)).*0.5,ones(400,1)],2,2).*...
            vecnorm([(-small(wndwS:wndwe,:)+large(wndwS:wndwe,:)).*0.5,ones(400,1)],2,2).*...
            vecnorm([sp1,zeros(400,1)],2,2).*vecnorm([sp2,zeros(400,1)],2,2));
        temp(temp>1)=NaN;
        AngMoment.Mismatch = [AngMoment.Mismatch;nanmean(temp,'all')];
    end
%     pol(end+1)=nanmean(temp,'all');
    
end

ctS=0;
for smallCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
    ctS=ctS+1;
    if isequal(TrialFishIds.Small(smallCt),'1')
        small = ProxyTrajectories(smallCt).fish1Pos;
        large = ProxyTrajectories(smallCt).fish2Pos;
    elseif isequal(TrialFishIds.Small(smallCt),'2')
        small = ProxyTrajectories(smallCt).fish2Pos;
        large = ProxyTrajectories(smallCt).fish1Pos;
    end
    randIdx = round(rand*11+1);
    largeCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1);
        while smallCt==largeCt(randIdx)
            randIdx = round(rand*11+1);
        end
        if isequal(TrialFishIds.Small(largeCt(randIdx)),'1')
            large = ProxyTrajectories(largeCt(randIdx)).fish2Pos;
        elseif isequal(TrialFishIds.Small(largeCt(randIdx)),'2')
            large = ProxyTrajectories(largeCt(randIdx)).fish1Pos;
        end
%         temp=[];
        for tt=1:30
            dice = round(rand(2,1)*29+1);
            wndwS = (dice(1)-1)*400 + 1;
            wndwe = (dice(1))*400;
            wndwS2 = (dice(2)-1)*400 + 1;
            wndwe2 = (dice(2))*400;
            pol.Chance = [pol.Chance;nanmean(vecnorm((diff(small(wndwS:wndwe,:),1,1)./vecnorm(diff(small(wndwS:wndwe,:),1,1),2,2))+...
            (diff(large(wndwS2:wndwe2,:),1,1)./vecnorm(diff(large(wndwS2:wndwe2,:),1,1),2,2)),2,2)./2)];
            DistBet.Chance = [DistBet.Chance;nanmean(sqrt((small(wndwS:wndwe,1)-large(wndwS2:wndwe2,1)).^2+(small(wndwS:wndwe,2)-large(wndwS2:wndwe2,2)).^2))];
            sp1 = [0,0;diff(small(wndwS:wndwe,:),1,1)].*20;
            sp2 = [0,0;diff(large(wndwS2:wndwe2,:),1,1)].*20;
            temp = vecnorm(cross([(small(wndwS:wndwe,:)-large(wndwS2:wndwe2,:)).*0.5,zeros(400,1)],[sp1,zeros(400,1)])+...
                cross([(-small(wndwS:wndwe,:)+large(wndwS2:wndwe2,:)).*0.5,zeros(400,1)],[sp2,zeros(400,1)]),2,2)./...
                (vecnorm([(small(wndwS:wndwe,:)-large(wndwS2:wndwe2,:)).*0.5,ones(400,1)],2,2).*...
                vecnorm([(-small(wndwS:wndwe,:)+large(wndwS2:wndwe2,:)).*0.5,ones(400,1)],2,2).*...
                vecnorm([sp1,zeros(400,1)],2,2).*vecnorm([sp2,zeros(400,1)],2,2));
            temp(temp>1)=NaN;
            AngMoment.Chance = [AngMoment.Chance;nanmean(temp,'all')];
        end   
%         pol(end+1)=nanmean(temp,'all');        
end
clearvars,clc;
load('FigS3data.mat')
figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 8 6])
tt=tiledlayout(2,2)
tt.Padding = 'compact'

nexttile
tbl = table(fillmissing(pol.Control(:),'constant',0),fillmissing(DistBet.Control(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
% pbaspect([1 1 1])
% xlabel('Polarization')
ylabel('Distance (cm)')
title('Control');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(pol.Match(:),'constant',0),fillmissing(DistBet.Match(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
% pbaspect([1 1 1])
% xlabel('Polarization')
ylabel('')
title('Match');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(pol.Mismatch(:),'constant',0),fillmissing(DistBet.Mismatch(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
% pbaspect([1 1 1])
xlabel('Polarization')
ylabel('Distance (cm)')
title('Mismatch');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(pol.Chance(:),'constant',0),fillmissing(DistBet.Chance(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
% pbaspect([1 1 1])
xlabel('Polarization')
ylabel('')
title('Chance');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

% figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 7 5])
% tt=tiledlayout(2,2)
% tt.Padding = 'compact'

nexttile
tbl = table(fillmissing(pol.Control(:),'constant',0),fillmissing(AngMoment.Control(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 1])
pbaspect([1 1 1])
xlabel('Polarization')
ylabel('Angular momentum')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(pol.Match(:),'constant',0),fillmissing(AngMoment.Match(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 1])
pbaspect([1 1 1])
xlabel('Polarization')
ylabel('')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(pol.Mismatch(:),'constant',0),fillmissing(AngMoment.Mismatch(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 1])
pbaspect([1 1 1])
xlabel('Polarization')
ylabel('')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(pol.Chance(:),'constant',0),fillmissing(AngMoment.Chance(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 1])
pbaspect([1 1 1])
xlabel('Polarization')
ylabel('')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

% figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 7 5])
% tt=tiledlayout(2,2)
% tt.Padding = 'compact'

nexttile
tbl = table(fillmissing(AngMoment.Control(:),'constant',0),fillmissing(DistBet.Control(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
pbaspect([1 1 1])
xlabel('Angular momentum')
ylabel('Distance (cm)')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(AngMoment.Match(:),'constant',0),fillmissing(DistBet.Match(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
pbaspect([1 1 1])
xlabel('Angular momentum')
ylabel('')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(AngMoment.Mismatch(:),'constant',0),fillmissing(DistBet.Mismatch(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
pbaspect([1 1 1])
xlabel('Angular momentum')
ylabel('')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off

nexttile
tbl = table(fillmissing(AngMoment.Chance(:),'constant',0),fillmissing(DistBet.Chance(:),'constant',0));
mdl = fitlm(tbl,'linear')
h = plot(mdl);
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
axis([0 1 0 30])
pbaspect([1 1 1])
xlabel('Angular momentum')
ylabel('')
title('');
box off
set(gca,'TickDir','out','LineWidth',1,'FontName','Arial','FontSize',11)
legend off