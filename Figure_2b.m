%% Time delay analysis
clearvars;
load('Figure2b.mat');
clc
cntr = 0;
delayAnl = struct;
for ct=1:34  
    if isequal(ProxyTrajectories(ct).Condition,{'Control'})
        continue
    end
    cntr = cntr+1;
    for tt=1:30
        wndwS = (tt-1)*400 + 1;
        wndwe = (tt)*400;
        %% Left fish right robot
        [delayAnl(cntr).XrCorr(tt,1),delayAnl(cntr).Lag(tt,1)]=xcorr_windowed(ProxyTrajectories(ct).fishRPos(wndwS:wndwe,:),ProxyTrajectories(ct).robotLPos(wndwS:wndwe,:),10,0.05);
        delayAnl(cntr).DistBet(tt,1) = nanmean(sqrt((ProxyTrajectories(ct).fishLPos(wndwS:wndwe,1)-ProxyTrajectories(ct).robotLPos(wndwS:wndwe,1)).^2+(ProxyTrajectories(ct).fishLPos(wndwS:wndwe,2)-ProxyTrajectories(ct).robotLPos(wndwS:wndwe,2)).^2));
        
        %% Right fish left robot
        [delayAnl(cntr).XrCorr(tt,2),delayAnl(cntr).Lag(tt,2)]=xcorr_windowed(ProxyTrajectories(ct).fishLPos(wndwS:wndwe,:),ProxyTrajectories(ct).robotRPos(wndwS:wndwe,:),10,0.05);
        delayAnl(cntr).DistBet(tt,2) = nanmean(sqrt((ProxyTrajectories(ct).fishRPos(wndwS:wndwe,1)-ProxyTrajectories(ct).robotRPos(wndwS:wndwe,1)).^2+(ProxyTrajectories(ct).fishRPos(wndwS:wndwe,2)-ProxyTrajectories(ct).robotRPos(wndwS:wndwe,2)).^2));
        
    end
end

for ii = 1:22
    cc_lag(ii,:,1) = [mean(delayAnl(ii).XrCorr(:,1)),mean(delayAnl(ii).Lag(:,1))];
    cc_lag(ii,:,2) = [mean(delayAnl(ii).XrCorr(:,2)),mean(delayAnl(ii).Lag(:,2))];
    cc_lag_err(ii,:,1) = [std(delayAnl(ii).XrCorr(:,1))/sqrt(30),std(delayAnl(ii).Lag(:,1))/sqrt(30)];
    cc_lag_err(ii,:,2) = [std(delayAnl(ii).XrCorr(:,2))/sqrt(30),std(delayAnl(ii).Lag(:,2))/sqrt(30)];
end
figure('Renderer', 'painters', 'Units', 'Inches', 'Position', [0.5 0.5 7 4])
rectangle('Position',[-0.2 .95 0.2 0.05],'FaceColor','#7ede7e','LineStyle','none')
hold on
errorbar(cc_lag(:,2,1),cc_lag(:,1,1),-cc_lag_err(:,1,1),cc_lag_err(:,1,1),-cc_lag_err(:,2,1),cc_lag_err(:,2,1),'ok','MarkerSize',5,'CapSize',0,'MarkerFaceColor','k');
errorbar(cc_lag(:,2,2),cc_lag(:,1,2),-cc_lag_err(:,1,2),cc_lag_err(:,1,2),-cc_lag_err(:,2,2),cc_lag_err(:,2,2),'ok','MarkerSize',5,'CapSize',0,'MarkerFaceColor','k');
hold off
axis([-0.4 0 0.75 1])
pbaspect([0.4 .25 1])
box off
xlabel('Lag (s)')
ylabel('Similarity index')
% title('Cross correlation vs lags according to 20s windows')
% legend({'Left robot','right robot'},'Location','southwest');
set(gca,'FontName','Arial','FontSize',11,'TickDir','out','LineWidth',1)