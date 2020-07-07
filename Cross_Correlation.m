%% Cross Correlation testing
clearvars
clc
load('Data.mat');
ds=2

    ctS = 0;
    for smallCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
        ctL=0;
        ctS = ctS+1;
        if isequal(TrialFishIds.Small(smallCt),'1')
            small = ProxyTrajectories(smallCt).fish1TurnRate;
        elseif isequal(TrialFishIds.Small(smallCt),'2')
            small = ProxyTrajectories(smallCt).fish2TurnRate;            
        end
        small(isnan(small))=0;
        for largeCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
            ctL=ctL+1;
            if isequal(TrialFishIds.Small(largeCt),'1')
                large = ProxyTrajectories(largeCt).fish2TurnRate;
            elseif isequal(TrialFishIds.Small(largeCt),'2')
                large = ProxyTrajectories(largeCt).fish1TurnRate; 
            end
            large(isnan(large))=0;
            tot_ccorr = 0;
            tot_lag = 0;
            tot_ccorr2 = 0;
            tot_lag2 = 0;
            n_windows = 0;
            if smallCt==largeCt        
                for tt=1:30
                    wndwSsm = (tt-1)*400 + 1;
                    wndwesm = (tt)*400;
                    wndwSl = (tt-1)*400 + 1;
                    wndwel = (tt)*400;
                    n_windows = n_windows+1; 
                    [temp_ccor,temp_lag] = xcorr(downsample(small(wndwSsm:wndwesm,:),ds),...
                        downsample(large(wndwSl:wndwel,:),ds),round(20/ds),'normalized');
                    [temp_mx_ccor,idx] = nanmax(temp_ccor);
                    if isnan(temp_mx_ccor)
                        n_windows = n_windows-1;
                        continue
                    end
                    temp_mx_lag = temp_lag(idx);
                    tot_ccorr = nansum([tot_ccorr, temp_mx_ccor]);
                    tot_lag = tot_lag + temp_mx_lag;
                    clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                    [temp_ccor,temp_lag] = xcorr(downsample(large(wndwSl:wndwel,:),ds),...
                        downsample(small(wndwSsm:wndwesm,:),ds),...
                        round(20/ds),'normalized');
                    [temp_mx_ccor,idx] = nanmax(temp_ccor);
                    temp_mx_lag = temp_lag(idx);
                    tot_ccorr2 = nansum([tot_ccorr2, temp_mx_ccor]);
                    tot_lag2 = tot_lag2 + temp_mx_lag;
                    clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                end
            else
                for tt=1:30
                    wndwSsm = (tt-1)*400 + 1;
                    wndwesm = (tt)*400;
                    for ttt=1:30
                        wndwSl = (ttt-1)*400 + 1;
                        wndwel = (ttt)*400;
                        n_windows = n_windows+1; 
                        [temp_ccor,temp_lag] = xcorr(downsample(small(wndwSsm:wndwesm,:),ds),...
                            downsample(large(wndwSl:wndwel,:),ds),round(20/ds),'normalized');
                        [temp_mx_ccor,idx] = nanmax(temp_ccor);
                        if isnan(temp_mx_ccor)
                            n_windows = n_windows-1;
                            continue
                        end
                        temp_mx_lag = temp_lag(idx);
                        tot_ccorr = nansum([tot_ccorr, temp_mx_ccor]);
                        tot_lag = tot_lag + temp_mx_lag;
                        clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                        [temp_ccor,temp_lag] = xcorr(downsample(large(wndwSl:wndwel,:),ds),...
                            downsample(small(wndwSsm:wndwesm,:),ds),...
                            round(20/ds),'normalized');
                        [temp_mx_ccor,idx] = nanmax(temp_ccor);
                        temp_mx_lag = temp_lag(idx);
                        tot_ccorr2 = nansum([tot_ccorr2, temp_mx_ccor]);
                        tot_lag2 = tot_lag2 + temp_mx_lag;
                        clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                    end
                end  
            end
            CCorr(ctS,ctL,1) = tot_ccorr/n_windows;
            Lags(ctS,ctL,1) = tot_lag/n_windows;
            CCorr2(ctS,ctL,1) = tot_ccorr2/n_windows;
            Lags2(ctS,ctL,1) = tot_lag2/n_windows;
        end
    end

    ctS=0;
    for smallCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
        ctL=0;
        ctS = ctS+1;
        if isequal(TrialFishIds.Small(smallCt),'L')
            small = ProxyTrajectories(smallCt).fishLTurnRate;
        elseif isequal(TrialFishIds.Small(smallCt),'R')
            small = ProxyTrajectories(smallCt).fishRTurnRate;
        end
        small(isnan(small))=0;
        for largeCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
            ctL=ctL+1;
            if isequal(TrialFishIds.Small(largeCt),'L')
                large = ProxyTrajectories(largeCt).fishRTurnRate;
            elseif isequal(TrialFishIds.Small(largeCt),'R')
                large = ProxyTrajectories(largeCt).fishLTurnRate; 
            end
            large(isnan(large))=0;
            tot_ccorr = 0;
            tot_lag = 0;
            tot_ccorr2 = 0;
            tot_lag2 = 0;
            n_windows = 0;
            smallDelay = find(delayAnl(smallCt).Simultaneous==1);
            largeDelay = find(delayAnl(largeCt).Simultaneous==1);
            min_bin = min(sum(delayAnl(smallCt).Simultaneous),sum(delayAnl(largeCt).Simultaneous));
            if smallCt==largeCt        
                for tt=1:min_bin
                    wndwSsm = (smallDelay(tt)-1)*400 + 1;
                    wndwesm = (smallDelay(tt))*400;
                    wndwSl = (largeDelay(tt)-1)*400 + 1;
                    wndwel = (largeDelay(tt))*400;
                    n_windows = n_windows+1;
                    sml = downsample(small(wndwSsm:wndwesm,:),ds);
                    lrg = downsample(large(wndwSl:wndwel,:),ds);
                    [temp_ccor,temp_lag] = xcorr(sml(1:end-2),...
                        lrg(3:end),round(20/ds),'normalized');
                    [temp_mx_ccor,idx] = nanmax(temp_ccor);
                    if isnan(temp_mx_ccor)
                        n_windows = n_windows-1;
                        continue
                    end
                    temp_mx_lag = temp_lag(idx);
                    tot_ccorr = nansum([tot_ccorr, temp_mx_ccor]);
                    tot_lag = tot_lag + temp_mx_lag;
                    clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                    [temp_ccor,temp_lag] = xcorr(lrg(1:end-2),...
                            sml(3:end),...
                            round(20/ds),'normalized');
                    [temp_mx_ccor,idx] = nanmax(temp_ccor);
                    temp_mx_lag = temp_lag(idx);
                    tot_ccorr2 = nansum([tot_ccorr2, temp_mx_ccor]);
                    tot_lag2 = tot_lag2 + temp_mx_lag;
                    clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                end
            else
                for tt=1:numel(smallDelay)
                    wndwSsm = (smallDelay(tt)-1)*400 + 1;
                    wndwesm = (smallDelay(tt))*400;
                    for ttt=1:numel(largeDelay)
                        wndwSl = (largeDelay(ttt)-1)*400 + 1;
                        wndwel = (largeDelay(ttt))*400;
                        n_windows = n_windows+1;
                        sml = downsample(small(wndwSsm:wndwesm,:),ds);
                        lrg = downsample(large(wndwSl:wndwel,:),ds);
                        [temp_ccor,temp_lag] = xcorr(sml(1:end-2),...
                            lrg(3:end),round(20/ds),'normalized');
                        [temp_mx_ccor,idx] = nanmax(temp_ccor);
                        if isnan(temp_mx_ccor)
                            n_windows = n_windows-1;
                            continue
                        end
                        temp_mx_lag = temp_lag(idx);
                        tot_ccorr = nansum([tot_ccorr, temp_mx_ccor]);
                        tot_lag = tot_lag + temp_mx_lag;
                        clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                        [temp_ccor,temp_lag] = xcorr(lrg(1:end-2),...
                            sml(3:end),...
                            round(20/ds),'normalized');
                        [temp_mx_ccor,idx] = nanmax(temp_ccor);
                        temp_mx_lag = temp_lag(idx);
                        tot_ccorr2 = nansum([tot_ccorr2, temp_mx_ccor]);
                        tot_lag2 = tot_lag2 + temp_mx_lag;
                        clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                    end
                end
            end
            CCorr(ctS,ctL,2) = tot_ccorr/n_windows;
            Lags(ctS,ctL,2) = tot_lag/n_windows;
            CCorr2(ctS,ctL,2) = tot_ccorr2/n_windows;
            Lags2(ctS,ctL,2) = tot_lag2/n_windows;
        end
    end

    ctS=0;
    for smallCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
        ctL=0;
        ctS = ctS+1;
        if isequal(TrialFishIds.Small(smallCt),'L')
            small = ProxyTrajectories(smallCt).fishLTurnRate;
        elseif isequal(TrialFishIds.Small(smallCt),'R')
            small = ProxyTrajectories(smallCt).fishRTurnRate;
        end
        small(isnan(small))=0;
        for largeCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
            ctL=ctL+1;
            if isequal(TrialFishIds.Small(largeCt),'L')
                large = ProxyTrajectories(largeCt).fishRTurnRate;
            elseif isequal(TrialFishIds.Small(largeCt),'R')
                large = ProxyTrajectories(largeCt).fishLTurnRate; 
            end
            large(isnan(large))=0;
            tot_ccorr = 0;
            tot_lag = 0;
            tot_ccorr2 = 0;
            tot_lag2 = 0;
            n_windows = 0;
            smallDelay = find(delayAnl(smallCt).Simultaneous==1);
            largeDelay = find(delayAnl(largeCt).Simultaneous==1);
            min_bin = min(sum(delayAnl(smallCt).Simultaneous),sum(delayAnl(largeCt).Simultaneous));
            if smallCt==largeCt        
                for tt=1:min_bin
                    wndwSsm = (smallDelay(tt)-1)*400 + 1;
                    wndwesm = (smallDelay(tt))*400;
                    wndwSl = (largeDelay(tt)-1)*400 + 1;
                    wndwel = (largeDelay(tt))*400;
                    n_windows = n_windows+1;
                    sml = downsample(small(wndwSsm:wndwesm,:),ds);
                    lrg = downsample(large(wndwSl:wndwel,:),ds);
                    [temp_ccor,temp_lag] = xcorr(sml(1:end-1),...
                        lrg(2:end),round(20/ds),'normalized');
                    [temp_mx_ccor,idx] = nanmax(temp_ccor);
                    if isnan(temp_mx_ccor)
                        n_windows = n_windows-1;
                        continue
                    end
                    temp_mx_lag = temp_lag(idx);
                    tot_ccorr = nansum([tot_ccorr, temp_mx_ccor]);
                    tot_lag = tot_lag + temp_mx_lag;
                    clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                    [temp_ccor,temp_lag] = xcorr(lrg(1:end-1),...
                            sml(2:end),...
                            round(20/ds),'normalized');
                    [temp_mx_ccor,idx] = nanmax(temp_ccor);
                    temp_mx_lag = temp_lag(idx);
                    tot_ccorr2 = nansum([tot_ccorr2, temp_mx_ccor]);
                    tot_lag2 = tot_lag2 + temp_mx_lag;
                    clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                end
            else
                for tt=1:numel(smallDelay)
                    wndwSsm = (smallDelay(tt)-1)*400 + 1;
                    wndwesm = (smallDelay(tt))*400;
                    for ttt=1:numel(largeDelay)
                        wndwSl = (largeDelay(ttt)-1)*400 + 1;
                        wndwel = (largeDelay(ttt))*400;
                        n_windows = n_windows+1;
                        sml = downsample(small(wndwSsm:wndwesm,:),ds);
                        lrg = downsample(large(wndwSl:wndwel,:),ds);
                        [temp_ccor,temp_lag] = xcorr(sml(1:end-1),...
                            lrg(2:end),round(20/ds),'normalized');
                        [temp_mx_ccor,idx] = nanmax(temp_ccor);
                        temp_mx_lag = temp_lag(idx);
                        if isnan(temp_mx_ccor)
                            n_windows = n_windows-1;
                            continue
                        end
                        tot_ccorr = nansum([tot_ccorr, temp_mx_ccor]);
                        tot_lag = tot_lag + temp_mx_lag;
                        clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                        [temp_ccor,temp_lag] = xcorr(lrg(1:end-1),...
                            sml(2:end),...
                            round(20/ds),'normalized');
                        [temp_mx_ccor,idx] = nanmax(temp_ccor);
                        temp_mx_lag = temp_lag(idx);
                        tot_ccorr2 = nansum([tot_ccorr2, temp_mx_ccor]);
                        tot_lag2 = tot_lag2 + temp_mx_lag;
                        clear temp_ccor temp_lag temp_mx_ccor temp_mx_lag idx
                    end
                end
            end
            CCorr(ctS,ctL,3) = tot_ccorr/n_windows;
            Lags(ctS,ctL,3) = tot_lag/n_windows;
            CCorr2(ctS,ctL,3) = tot_ccorr2/n_windows;
            Lags2(ctS,ctL,3) = tot_lag2/n_windows;
        end
    end
    Lags = Lags.*(ds/20);
    Lags2 = Lags2.*(ds/20);
    %% Figures
    ymax = .2;
    figure
    tiledlayout(3,2)
    %Control
    nexttile
    xlabl = 'Normalized cross correlation';
    [real, surr, ~] = ...
        func_perm_stats_CCorr(20000,CCorr(:,:,1), xlabl, ymax,1);
    title('Control')
%     [h,pvalue,~] = ztest(mean(real),nanmean(surr,'all'),nanstd(surr))
    'Control'
    1-invprctile(surr,mean(real))/100

    nexttile
    xlabl = 'Time Lags (s)';
    [real, surr, ~] = ...
        func_perm_stats_CCorr(20000,Lags(:,:,1), xlabl, ymax,2);
%     [h,pvalue,~] = ztest(mean(real),nanmean(surr,'all'),nanstd(surr))
    1-invprctile(surr,mean(real))/100
    title('Control')
%     scatter(diag(Lags(:,:,1)),diag(CCorr(:,:,1)),75,'k.');
%     hold on
%     plot([0,0],[0.5,1],'g--')
%     xlabel('Time lags (s)');
%     ylabel('Cross correlation');
%     axis([-1 1 0 1])

    %Match
    nexttile
    xlabl = 'Normalized cross correlation';
    [real, surr, ~] = ...
        func_perm_stats_CCorr(20000,CCorr(1:11,1:11,2), xlabl, ymax,1);
%     [h,pvalue,~] = ztest(mean(real(1:11)),nanmean(surr,'all'),nanstd(surr))
    'Match'
    1-invprctile(surr,mean(real(1:11)))/100
    title('Match')

    nexttile
    xlabl = 'Time Lags (s)';
    [real, surr, ~] = ...
        func_perm_stats_CCorr(20000,Lags(1:11,1:11,2), xlabl, ymax,2);
%     [h,pvalue,~] = ztest(mean(real(1:11)),nanmean(surr,'all'),nanstd(surr))
    1-invprctile(surr,mean(real(1:11)))/100
    title('Match')
%     scatter(diag(Lags(1:11,1:11,2)),diag(CCorr(1:11,1:11,2)),75,'k.');
%     hold on
%     plot([0,0],[0.5,1],'g--')
%     xlabel('Time lags (s)');
%     ylabel('Cross correlation');
%     axis([-1 1 0 1])

    %Mismatch
    nexttile
    xlabl = 'Normalized cross correlation';
    [real, surr, ~] = ...
        func_perm_stats_CCorr(20000,CCorr(1:11,1:11,3), xlabl, ymax,1);
%     [h,pvalue,~] = ztest(mean(real(1:11)),nanmean(surr,'all'),nanstd(surr))
    'Mismatch'
    1-invprctile(surr,mean(real(1:11)))/100
    title('Mismatch')

    nexttile
    xlabl = 'Time Lags (s)';
    [real, surr, ~] = ...
        func_perm_stats_CCorr(20000,Lags(1:11,1:11,3), xlabl, ymax,2);
%     [h,pvalue,~] = ztest(mean(real(1:11)),nanmean(surr,'all'),nanstd(surr))
    1-invprctile(surr,mean(real(1:11)))/100
    title('Mismatch')
%     scatter(diag(Lags(1:11,1:11,3)),diag(CCorr(1:11,1:11,3)),75,'k.');
%     hold on
%     plot([0,0],[0.5,1],'g--')
%     xlabel('Time lags (s)');
%     ylabel('Cross correlation');
%     axis([-1 1 0 1])

    suptitle(['Cross correlation and lag on TurnRate with sampling interval ',num2str(0.05*ds),' seconds'])




