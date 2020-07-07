clearvars
clc
load('Data.mat');
ds=2;
ctS = 0;
tau = 0;
for smallCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
    ctL=0;
    ctS = ctS+1
    if isequal(TrialFishIds.Small(smallCt),'1')
        smallSp = ProxyTrajectories(smallCt).fish1Speed;
        smallTR = ProxyTrajectories(smallCt).fish1TurnRate;
    elseif isequal(TrialFishIds.Small(smallCt),'2')
        smallSp = ProxyTrajectories(smallCt).fish2Speed;
        smallTR = ProxyTrajectories(smallCt).fish2TurnRate;
    end
    for largeCt=find(ismember([ProxyTrajectories.Condition],{'Control'})==1)
        ctL=ctL+1;
        if isequal(TrialFishIds.Small(largeCt),'1')
            largeSp = ProxyTrajectories(largeCt).fish2Speed;
            largeTR = ProxyTrajectories(largeCt).fish2TurnRate;
        elseif isequal(TrialFishIds.Small(largeCt),'2')
            largeSp = ProxyTrajectories(largeCt).fish1Speed; 
            largeTR = ProxyTrajectories(largeCt).fish1TurnRate; 
        end
        pmfXpX=zeros(1,16);
        pmfXpXY=zeros(1,64);
        pmfX=zeros(1,4);
        pmfXY=zeros(1,16);
        pmfYpY=zeros(1,16);
        pmfYpYX=zeros(1,64);
        pmfY=zeros(1,4);
        pmfYX=zeros(1,16);
        n_windows = 0;
        if smallCt==largeCt        
            for tt=1:30
                wndwSsm = (tt-1)*400 + 2;
                wndwesm = (tt)*400;
                wndwSl = (tt-1)*400 + 2;
                wndwel = (tt)*400;
                n_windows = n_windows+1;
                % Y to X
                [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp(wndwSl:wndwel,:),...
                    largeTR(wndwSl:wndwel,:),smallSp(wndwSsm:wndwesm,:),smallTR(wndwSsm:wndwesm,:),ds,tau);
%                 [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp,...
%                     largeTR,smallSp,smallTR,ds,tau);
                pmfXpX = pmfXpX + temp_pmfXpX;
                pmfXpXY = pmfXpXY + temp_pmfXpXY;
                pmfX = pmfX + temp_pmfX;
                pmfXY = pmfXY + temp_pmfXY;
                % X to Y
                [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp(wndwSsm:wndwesm,:),...
                    smallTR(wndwSsm:wndwesm,:),largeSp(wndwSl:wndwel,:),largeTR(wndwSl:wndwel,:),ds,tau);
%                 [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp,...
%                     smallTR,largeSp,largeTR,ds,tau);
                pmfYpY = pmfYpY + temp_pmfYpY;
                pmfYpYX = pmfYpYX + temp_pmfYpYX;
                pmfY = pmfY + temp_pmfY;
                pmfYX = pmfYX + temp_pmfYX;
            end
        else
            for tt=1:30
                wndwSsm = (tt-1)*400 + 2;
                wndwesm = (tt)*400;
                for ttt=1:30
                    wndwSl = (ttt-1)*400 + 2;
                    wndwel = (ttt)*400;
                    n_windows = n_windows+1;
                    % Y to X
                    [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp(wndwSl:wndwel,:),...
                    largeTR(wndwSl:wndwel,:),smallSp(wndwSsm:wndwesm,:),smallTR(wndwSsm:wndwesm,:),ds,tau);
%                     [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp,...
%                         largeTR,smallSp,smallTR,ds,tau);
                    pmfXpX = pmfXpX + temp_pmfXpX;
                    pmfXpXY = pmfXpXY + temp_pmfXpXY;
                    pmfX = pmfX + temp_pmfX;
                    pmfXY = pmfXY + temp_pmfXY;
                    % X to Y
                    [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp(wndwSsm:wndwesm,:),...
                    smallTR(wndwSsm:wndwesm,:),largeSp(wndwSl:wndwel,:),largeTR(wndwSl:wndwel,:),ds,tau);
% [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp,...
%                     smallTR,largeSp,largeTR,ds,tau);
                    pmfYpY = pmfYpY + temp_pmfYpY;
                    pmfYpYX = pmfYpYX + temp_pmfYpYX;
                    pmfY = pmfY + temp_pmfY;
                    pmfYX = pmfYX + temp_pmfYX;
                end
            end
        end
        pmfXpX = pmfXpX./n_windows;
        pmfXpXY = pmfXpXY./n_windows;
        pmfX = pmfX./n_windows;
        pmfXY = pmfXY./n_windows;

        pmfYpY = pmfYpY./n_windows;
        pmfYpYX = pmfYpYX./n_windows;
        pmfY = pmfY./n_windows;
        pmfYX = pmfYX./n_windows;
        
        [teSmalltoLarge_Control(ctS,ctL,1)] = symbolicTE_pmf(pmfXpX,pmfXpXY,pmfX,pmfXY);
        [teLargetoSmall_Control(ctS,ctL,1)] = symbolicTE_pmf(pmfYpY,pmfYpYX,pmfY,pmfYX);
    end
end

ctS = 0;
tau = 2;
for smallCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
    ctL=0;
    ctS = ctS+1
    if isequal(TrialFishIds.Small(smallCt),'L')
        smallSp = ProxyTrajectories(smallCt).fishLSpeed;
        smallTR = ProxyTrajectories(smallCt).fishLTurnRate;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        smallSp = ProxyTrajectories(smallCt).fishRSpeed;
        smallTR = ProxyTrajectories(smallCt).fishRTurnRate;
    end
    for largeCt=find(ismember([ProxyTrajectories.Condition],{'SameSize'})==1)
        ctL=ctL+1;
        if isequal(TrialFishIds.Small(largeCt),'L')
            largeSp = ProxyTrajectories(largeCt).fishRSpeed;
            largeTR = ProxyTrajectories(largeCt).fishRTurnRate;
        elseif isequal(TrialFishIds.Small(largeCt),'R')
            largeSp = ProxyTrajectories(largeCt).fishLSpeed; 
            largeTR = ProxyTrajectories(largeCt).fishLTurnRate;
        end
        pmfXpX=zeros(1,16);
        pmfXpXY=zeros(1,64);
        pmfX=zeros(1,4);
        pmfXY=zeros(1,16);
        pmfYpY=zeros(1,16);
        pmfYpYX=zeros(1,64);
        pmfY=zeros(1,4);
        pmfYX=zeros(1,16);
        n_windows = 0;
        smallDelay = find(delayAnl(smallCt).Simultaneous==1);
        largeDelay = find(delayAnl(largeCt).Simultaneous==1);
        min_bin = min(sum(delayAnl(smallCt).Simultaneous),sum(delayAnl(largeCt).Simultaneous));
        if smallCt==largeCt        
            for tt=1:min_bin
                wndwSsm = (smallDelay(tt)-1)*400 + 2;
                wndwesm = (smallDelay(tt))*400;
                wndwSl = (largeDelay(tt)-1)*400 + 2;
                wndwel = (largeDelay(tt))*400;
                n_windows = n_windows+1;
                % Y to X
                [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp(wndwSl:wndwel,:),...
                    largeTR(wndwSl:wndwel,:),smallSp(wndwSsm:wndwesm,:),smallTR(wndwSsm:wndwesm,:),ds,tau);
                pmfXpX = pmfXpX + temp_pmfXpX;
                pmfXpXY = pmfXpXY + temp_pmfXpXY;
                pmfX = pmfX + temp_pmfX;
                pmfXY = pmfXY + temp_pmfXY;
                % X to Y
                [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp(wndwSsm:wndwesm,:),...
                    smallTR(wndwSsm:wndwesm,:),largeSp(wndwSl:wndwel,:),largeTR(wndwSl:wndwel,:),ds,tau);
                pmfYpY = pmfYpY + temp_pmfYpY;
                pmfYpYX = pmfYpYX + temp_pmfYpYX;
                pmfY = pmfY + temp_pmfY;
                pmfYX = pmfYX + temp_pmfYX;
            end
        else
            for tt=1:numel(smallDelay)
                wndwSsm = (smallDelay(tt)-1)*400 + 2;
                wndwesm = (smallDelay(tt))*400;
                for ttt=1:numel(largeDelay)
                    wndwSl = (largeDelay(ttt)-1)*400 + 2;
                    wndwel = (largeDelay(ttt))*400;
                    n_windows = n_windows+1;
                    % Y to X
                    [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp(wndwSl:wndwel,:),...
                    largeTR(wndwSl:wndwel,:),smallSp(wndwSsm:wndwesm,:),smallTR(wndwSsm:wndwesm,:),ds,tau);
                    pmfXpX = pmfXpX + temp_pmfXpX;
                    pmfXpXY = pmfXpXY + temp_pmfXpXY;
                    pmfX = pmfX + temp_pmfX;
                    pmfXY = pmfXY + temp_pmfXY;
                    % X to Y
                    [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp(wndwSsm:wndwesm,:),...
                    smallTR(wndwSsm:wndwesm,:),largeSp(wndwSl:wndwel,:),largeTR(wndwSl:wndwel,:),ds,tau);
                    pmfYpY = pmfYpY + temp_pmfYpY;
                    pmfYpYX = pmfYpYX + temp_pmfYpYX;
                    pmfY = pmfY + temp_pmfY;
                    pmfYX = pmfYX + temp_pmfYX;
                end
            end
        end
        pmfXpX = pmfXpX./n_windows;
        pmfXpXY = pmfXpXY./n_windows;
        pmfX = pmfX./n_windows;
        pmfXY = pmfXY./n_windows;

        pmfYpY = pmfYpY./n_windows;
        pmfYpYX = pmfYpYX./n_windows;
        pmfY = pmfY./n_windows;
        pmfYX = pmfYX./n_windows;
        
        [teSmalltoLarge(ctS,ctL,2)] = symbolicTE_pmf(pmfXpX,pmfXpXY,pmfX,pmfXY);
        [teLargetoSmall(ctS,ctL,2)] = symbolicTE_pmf(pmfYpY,pmfYpYX,pmfY,pmfYX);
    end
end        
        
        
ctS = 0;
tau=1;
for smallCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
    ctL=0;
    ctS = ctS+1
    if isequal(TrialFishIds.Small(smallCt),'L')
        smallSp = ProxyTrajectories(smallCt).fishLSpeed;
        smallTR = ProxyTrajectories(smallCt).fishLTurnRate;
    elseif isequal(TrialFishIds.Small(smallCt),'R')
        smallSp = ProxyTrajectories(smallCt).fishRSpeed;
        smallTR = ProxyTrajectories(smallCt).fishRTurnRate;
    end
    for largeCt=find(ismember([ProxyTrajectories.Condition],{'DiffSize'})==1)
        ctL=ctL+1;
        if isequal(TrialFishIds.Small(largeCt),'L')
            largeSp = ProxyTrajectories(largeCt).fishRSpeed;
            largeTR = ProxyTrajectories(largeCt).fishRTurnRate;
        elseif isequal(TrialFishIds.Small(largeCt),'R')
            largeSp = ProxyTrajectories(largeCt).fishLSpeed; 
            largeTR = ProxyTrajectories(largeCt).fishLTurnRate;
        end
        pmfXpX=zeros(1,16);
        pmfXpXY=zeros(1,64);
        pmfX=zeros(1,4);
        pmfXY=zeros(1,16);
        pmfYpY=zeros(1,16);
        pmfYpYX=zeros(1,64);
        pmfY=zeros(1,4);
        pmfYX=zeros(1,16);
        n_windows = 0;
        smallDelay = find(delayAnl(smallCt).Simultaneous==1);
        largeDelay = find(delayAnl(largeCt).Simultaneous==1);
        min_bin = min(sum(delayAnl(smallCt).Simultaneous),sum(delayAnl(largeCt).Simultaneous));
        if smallCt==largeCt        
            for tt=1:min_bin
                wndwSsm = (smallDelay(tt)-1)*400 + 2;
                wndwesm = (smallDelay(tt))*400;
                wndwSl = (largeDelay(tt)-1)*400 + 2;
                wndwel = (largeDelay(tt))*400;
                n_windows = n_windows+1;
                % Y to X
                [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp(wndwSl:wndwel,:),...
                    largeTR(wndwSl:wndwel,:),smallSp(wndwSsm:wndwesm,:),smallTR(wndwSsm:wndwesm,:),ds,tau);
                pmfXpX = pmfXpX + temp_pmfXpX;
                pmfXpXY = pmfXpXY + temp_pmfXpXY;
                pmfX = pmfX + temp_pmfX;
                pmfXY = pmfXY + temp_pmfXY;
                % X to Y
                [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp(wndwSsm:wndwesm,:),...
                    smallTR(wndwSsm:wndwesm,:),largeSp(wndwSl:wndwel,:),largeTR(wndwSl:wndwel,:),ds,tau);
                pmfYpY = pmfYpY + temp_pmfYpY;
                pmfYpYX = pmfYpYX + temp_pmfYpYX;
                pmfY = pmfY + temp_pmfY;
                pmfYX = pmfYX + temp_pmfYX;
            end
        else
            for tt=1:numel(smallDelay)
                wndwSsm = (smallDelay(tt)-1)*400 + 2;
                wndwesm = (smallDelay(tt))*400;
                for ttt=1:numel(largeDelay)
                    wndwSl = (largeDelay(ttt)-1)*400 + 2;
                    wndwel = (largeDelay(ttt))*400;
                    n_windows = n_windows+1;
                    % Y to X
                    [temp_pmfXpX,temp_pmfXpXY,temp_pmfX,temp_pmfXY] = probability_mass_paired(largeSp(wndwSl:wndwel,:),...
                    largeTR(wndwSl:wndwel,:),smallSp(wndwSsm:wndwesm,:),smallTR(wndwSsm:wndwesm,:),ds,tau);
                    pmfXpX = pmfXpX + temp_pmfXpX;
                    pmfXpXY = pmfXpXY + temp_pmfXpXY;
                    pmfX = pmfX + temp_pmfX;
                    pmfXY = pmfXY + temp_pmfXY;
                    % X to Y
                    [temp_pmfYpY,temp_pmfYpYX,temp_pmfY,temp_pmfYX] = probability_mass_paired(smallSp(wndwSsm:wndwesm,:),...
                    smallTR(wndwSsm:wndwesm,:),largeSp(wndwSl:wndwel,:),largeTR(wndwSl:wndwel,:),ds,tau);
                    pmfYpY = pmfYpY + temp_pmfYpY;
                    pmfYpYX = pmfYpYX + temp_pmfYpYX;
                    pmfY = pmfY + temp_pmfY;
                    pmfYX = pmfYX + temp_pmfYX;
                end
            end
        end
        pmfXpX = pmfXpX./n_windows;
        pmfXpXY = pmfXpXY./n_windows;
        pmfX = pmfX./n_windows;
        pmfXY = pmfXY./n_windows;

        pmfYpY = pmfYpY./n_windows;
        pmfYpYX = pmfYpYX./n_windows;
        pmfY = pmfY./n_windows;
        pmfYX = pmfYX./n_windows;
        
        [teSmalltoLarge(ctS,ctL,3)] = symbolicTE_pmf(pmfXpX,pmfXpXY,pmfX,pmfXY);
        [teLargetoSmall(ctS,ctL,3)] = symbolicTE_pmf(pmfYpY,pmfYpYX,pmfY,pmfYX);
    end
end
ymax = .4;
figure
tiledlayout(1,3)
%Control
nexttile
xlabl = 'TE_{Small\rightarrowLarge}';
[RealTE(:,1,1), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge_Control(:,:,1), xlabl, ymax,1);
title('Control')
% [h,pvalue,~] = ztest(mean(RealTE(:,1,1)),nanmean(surr,'all'),nanstd(surr))
'Control'
1-invprctile(surr,mean(RealTE(:,1,1)))/100

nexttile
xlabl = 'TE_{Large\rightarrowSmall}';
[RealTE(:,2,1), surr, ~] = ...
    func_perm_stats(20000,teLargetoSmall_Control(:,:,1), xlabl, ymax,1);
% [h,pvalue,~] = ztest(mean(RealTE(:,2,1)),nanmean(surr,'all'),nanstd(surr))
title('Control')
1-invprctile(surr,mean(RealTE(:,2,1)))/100

nexttile
xlabl = 'TE_{Small\rightarrowLarge} - TE_{Large\rightarrowSmall}';
[RealTE(:,3,1), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge_Control(:,:,1)-teLargetoSmall_Control(:,:,1), xlabl, ymax,2);
% [h,pvalue,~] = ztest(mean(RealTE(:,3,1)),nanmean(surr,'all'),nanstd(surr))
title(['Control DS ',num2str(ds)])      

ymax = .4;
figure
tiledlayout(3,3)
%Control
nexttile
xlabl = 'TE_{Small\rightarrowLarge}';
[RealTE(:,1,1), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge_Control(:,:,1), xlabl, ymax,1);
title('Control')
% [h,pvalue,~] = ztest(mean(RealTE(:,1,1)),nanmean(surr,'all'),nanstd(surr))
'Control'
1-invprctile(surr,mean(RealTE(:,1,1)))/100

nexttile
xlabl = 'TE_{Large\rightarrowSmall}';
[RealTE(:,2,1), surr, ~] = ...
    func_perm_stats(20000,teLargetoSmall_Control(:,:,1), xlabl, ymax,1);
% [h,pvalue,~] = ztest(mean(RealTE(:,2,1)),nanmean(surr,'all'),nanstd(surr))
title('Control')
1-invprctile(surr,mean(RealTE(:,2,1)))/100

nexttile
xlabl = 'TE_{Small\rightarrowLarge} - TE_{Large\rightarrowSmall}';
[RealTE(:,3,1), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge_Control(:,:,1)-teLargetoSmall_Control(:,:,1), xlabl, ymax,2);
% [h,pvalue,~] = ztest(mean(RealTE(:,3,1)),nanmean(surr,'all'),nanstd(surr))
title('Control')
'Net'
1-invprctile(surr,mean(RealTE(:,3,1)))/100
invprctile(surr,mean(RealTE(:,3,1)))/100
        
%Match
nexttile
xlabl = 'TE_{Small\rightarrowLarge}';
[RealTE(1:11,1,2), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge(:,:,2), xlabl, ymax,1);
% [h,pvalue,~] = ztest(mean(RealTE(1:11,1,2)),nanmean(surr,'all'),nanstd(surr))
title('Match')
'Match'
1-invprctile(surr,mean(RealTE(1:11,1,2)))/100

nexttile
xlabl = 'TE_{Large\rightarrowSmall}';
[RealTE(1:11,2,2), surr, ~] = ...
    func_perm_stats(20000,teLargetoSmall(:,:,2), xlabl, ymax,1);
% [h,pvalue,~] = ztest(mean(RealTE(1:11,2,2)),nanmean(surr,'all'),nanstd(surr))
1-invprctile(surr,mean(RealTE(1:11,2,2)))/100
title('Match')

nexttile
xlabl = 'TE_{Small\rightarrowLarge} - TE_{Large\rightarrowSmall}';
[RealTE(1:11,3,2), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge(:,:,2)-teLargetoSmall(:,:,2), xlabl, ymax,2);
% [h,pvalue,~] = ztest(mean(RealTE(1:11,3,2)),nanmean(surr,'all'),nanstd(surr))
'Net'
1-invprctile(surr,mean(RealTE(1:11,3,2)))/100
invprctile(surr,mean(RealTE(1:11,3,2)))/100
title('Match')

%Mismatch
nexttile
xlabl = 'TE_{Small\rightarrowLarge}';
[RealTE(1:11,1,3), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge(:,:,3), xlabl, ymax,1);
% [h,pvalue,~] = ztest(mean(RealTE(1:11,1,3)),nanmean(surr,'all'),nanstd(surr))
'Mismatch'
1-invprctile(surr,mean(RealTE(1:11,1,3)))/100
title('Mismatch')

nexttile
xlabl = 'TE_{Large\rightarrowSmall}';
[RealTE(1:11,2,3), surr, ~] = ...
    func_perm_stats(20000,teLargetoSmall(:,:,3), xlabl, ymax,1);
% [h,pvalue,~] = ztest(mean(RealTE(1:11,2,3)),nanmean(surr,'all'),nanstd(surr))
1-invprctile(surr,mean(RealTE(1:11,2,3)))/100
title('Mismatch')

nexttile
xlabl = 'TE_{Small\rightarrowLarge} - TE_{Large\rightarrowSmall}';
[RealTE(1:11,3,3), surr, ~] = ...
    func_perm_stats(20000,teSmalltoLarge(:,:,3)-teLargetoSmall(:,:,3), xlabl, ymax,2);
% [h,pvalue,~] = ztest(mean(RealTE(1:11,3,3)),nanmean(surr,'all'),nanstd(surr))
'Net'
1-invprctile(surr,mean(RealTE(1:11,3,3)))/100
invprctile(surr,mean(RealTE(1:11,3,3)))/100
title('Mismatch')

% suptitle(['TE on speed and turn rate combine with sampling interval ',num2str(0.05*ds),' seconds and tau=0 for control, tau=2 for mis/match'])
% end

% ttestTE