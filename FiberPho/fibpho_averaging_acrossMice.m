%Written 01/2024 by MRT
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%Function to combine fiber photometry data across animals

%input:
%   loads fiber photometry data from multiple mice
%output:
%   combines to get mean and SEM analysis and across time plots

%Manually run mode for each condition:
%  1. Save Figure 4 (top 6 animals) as .svg
%  2. copy/paste values from command line to TwoGroupsMeanDiff
%  3. save TwoGroupsMeanDiff as .svg
%        note: ymin/max=  [-6, +11]
%        symbol sizes all max

%% LOAD DATA

%NOTE: IF ANY OF THESE ARE CHANGED, clear data to reload GRP
clear;  %force a clearing of data
fcutoff = 3; %3 Hz   %lowpass each trial over time

GraphicMode=2; %Fit2m1 PDF, with 95CI of MixingTest
DS = 12; %dot size
TRIM = 0; %33;  %trimmean -- don't use it

%FLAGS
bExludeMice=0; %change to 0 to look at all mice
KEEPTOP=6; %change based on burst cutoff
bBootStrapSpectogram = 1;
 
%%CHANGE THIS BASED ON WHICH DAY/TRIAL TYPE YOU WANT TO ANALYZE
%See below switch/case section for information on available
%ExperimentalStage cases
ExperimentalStage = 'B1' 

%Analysis Windows
AnalysisWindows = { ...
    [-1.5  -0.5   ]  ...
    [ 0     1.5  ] ...
    };

TMINMAX = [-3 3]; %time window
switch ExperimentalStage
    case 'XCRIT'
        %all trials across training days 6-10
        Day= [6:10];
        Type=1;
        TRIALS = [];        
    case {'A01' 'A02' 'A03' 'A04' 'A05' 'A06' 'A07' 'A08' 'A09'  'A10'  'A11'}
        %one of day 1-11, training, first 20 trials
        Day= NumbersInString(ExperimentalStage);
        Type=1;
        TRIALS = [1:20];
    case {'AA01' 'AA02' 'AA03' 'AA04' 'AA05' 'AA06' 'AA07' 'AA08' 'AA09'  'AA10'}
        %one of day 1-10, all trials
        Day= NumbersInString(ExperimentalStage);
        Type=1;
        TRIALS = [];
    case 'B0'
        %first 1 extinction trial
        Day=11;
        Type=2;
        TRIALS = [1];
    case 'B1'
        %first 4 extinction trials
        Day=11;
        Type=2;
        TRIALS = [1:4];
    case 'B2'
        %first #37-40 extinction trials
        Day=11;
        Type=2;
        TRIALS = [37:40];
    case 'B3'
        %day 12 extinction 1-50
        Day=12;
        Type=1;
        TRIALS = [1:50];        
    case 'C1'
        %day 10 probe trails
        Day=10;
        Type=2;
        TRIALS = [];        
    case 'D1'
        %first 4 conditioning trials
        Day=11;
        Type=3;
        TRIALS = [1:4];
    case 'D2'
        %first #37-40 conditioning trials
        Day=11;
        Type=3;
        TRIALS = [37:40];
    case 'D3'
        %day 12 conditioning
        Day=12;
        Type=2;
        TRIALS = [1:50];
end

try
    GRP = SAVED.GRP;
catch
    for j=1:2
        if j==1
            %CTL mice
            cd '\Fiber Photometry\ddHTP controls'; %change to your full path
        else
            %EXPT mice
            cd '\Fiber Photometry\HTP'; %change to your full path
        end
        if bExludeMice
            %exlude files that start with _VA...
            D = dir('VA*dFF*.mat');
        else
            %include all files that start with VA and _VA...
            D = dir('*VA*dFF*.mat');
        end
        NumTrials = [];
        GRP(j).FIB = [];
        for q=1:length(Day)
            for k=1:length(D)
                fprintf([num2str(k) ' .. ']);
                load(D(k).name);
                
                GRP(j).ts = fibdat(1).ts; %load time window information
                ts = GRP(j).ts;
                fs = 1/mean(diff(ts)); %sampling frequency in Hz
                TMP = fibdat(Day(q)).trialType(Type).F465_dFF; %fiber pho dFF data
                for n=1:size(TMP,1)
                    %each trial
                    TMP(n,:) = 100*LowpassFilter(ts, TMP(n,:), fcutoff, 1);                    
                end
                
                %add data to larger structure
                if q==1
                    GRP(j).FIB(1:size(TMP,1),1:size(TMP,2),k) = TMP;
                    GRP(j).NumTrials(k) = size(TMP,1);
                else
                    GRP(j).FIB(GRP(j).NumTrials(k)+(1:size(TMP,1)),1:size(TMP,2),k) = TMP;
                    GRP(j).NumTrials(k) = GRP(j).NumTrials(k) + size(TMP,1);
                end
            end
            for k=1:length(D)
                GRP(j).FIB(GRP(j).NumTrials(k)+1:end,:,k) = NaN;
            end
        end
        fprintf(['\n']);
    end
    SAVED.GRP = GRP;
end

if isempty(TRIALS)
    %all trials
    MaxTrial = max([GRP(:).NumTrials]);
    for k=1:length(GRP)
        GRP(k).FIB(end+1:MaxTrial,:,:) = NaN;  %pad with NaNs so all GRP have the same # trials for later cat(3,...)
    end
    TRIALS = 1:MaxTrial;
else
    %subset of trials was specified
    MaxTrial = max(TRIALS);
    for k=1:length(GRP)
        GRP(k).FIB(MaxTrial+1:end,:,:) = [];  %truncate
    end
end

ALLFIB = cat(3, GRP(:).FIB);
TIME = GRP(1).ts;

I = find(TIME>TMINMAX(1) & TIME<TMINMAX(2));

ALLFIB = ALLFIB(:,I,:);
TIME = TIME(I);

%% Bootstrap
if bBootStrapSpectogram
    NumBoot = 1e3;
    BOOTDIFF = zeros(size(ALLFIB,1), size(ALLFIB,2), NumBoot);
    
    for k=1:NumBoot
        P = randperm(size(ALLFIB,3));
        a = P(1:size(GRP(1).FIB,3));
        b = P(size(GRP(1).FIB,3)+1:end);
        BOOTDIFF(:,:,k) = mean(ALLFIB(:,:,b),3) - mean(ALLFIB(:,:,a),3);
        if mod(k,5000)==0
            disp(k);
        end
    end
    P = 1:size(ALLFIB,3);
    a = P(1:size(GRP(1).FIB,3));
    b = P(size(GRP(1).FIB,3)+1:end);
    ACTUALDIFF = mean(ALLFIB(:,:,b),3) - mean(ALLFIB(:,:,a),3);
else
    P = 1:size(ALLFIB,3);
    a = P(1:size(GRP(1).FIB,3));
    b = P(size(GRP(1).FIB,3)+1:end);
end


%% PLOT heatmap IMAGES  (x=time, y=trial)
set(0,'DefaultFigureWindowStyle','docked');
figure(1); clf;
for j=1:2
    subplot(2,1,j);
    if j==1
        mice = a;
        title('CTL');
    else
        mice = b;
        title('GBZ');
    end
    imagesc(TIME,1:size(GRP(j).FIB,1), mean(ALLFIB(:,:,mice),3,'omitnan') );
    xlabel('t (sec)');
    ylabel('trial#');
    caxis([-0.02 0.08]*100);
    colormap(hot(300))
    %ylim([0 20]);
    colorbar;
    drawnow
end

%% Plot bootstrap difference
if bBootStrapSpectogram
    figure(2); clf;
    PVAL = mean(ACTUALDIFF<BOOTDIFF,3);
    imagesc(TIME, 1:size(PVAL,1), -log10(PVAL));
    caxis(-log10([0.1 0.001]))
    h=colorbar;
    TICK=[0.1 0.03 0.01 0.003 0.001]'; set(h,'Ticks',-log10(TICK), 'TickLabels',TICK)
end


%NOTE: in case bINCLUDE_ALL_MICE ends as 0,so that TMPA and TMPB are complete
TMPA = squeeze(mean(ALLFIB(TRIALS,:,a), 1, 'omitnan'))';  %average over trials for CTL mice
TMPB = squeeze(mean(ALLFIB(TRIALS,:,b), 1, 'omitnan'))';  %average over trials for GBZ mice

%%analysis of time windows
for k=1:length(AnalysisWindows)
    W = AnalysisWindows{k};
    I = find(TIME>=W(1) & TIME<=(W(2)));
    CTL_GBZ{k} = [mean(TMPA(:,I),2)  mean(TMPB(:,I),2)];
    
    %ALL mice
    [ALL.Fit1(k), ALL.Fit2(k), ALL.Fit2m1(k), ALL.MixingTest(k)] = BootFitDifferencePlot(TRIM, zeros(size(CTL_GBZ{k},1),1), CTL_GBZ{k}(:,1), ones(size(CTL_GBZ{k},1),1), CTL_GBZ{k}(:,2), [0 1]); 

    %TOP# mice are the ones we keep
    [TOP.Fit1(k), TOP.Fit2(k), TOP.Fit2m1(k), TOP.MixingTest(k)] = BootFitDifferencePlot(TRIM, zeros(KEEPTOP,1), CTL_GBZ{k}(1:KEEPTOP,1), ones(KEEPTOP,1), CTL_GBZ{k}(1:KEEPTOP,2), [0 1]); 
end

%%
%Plot two figures: 
%Figure 3: all mice
%Figure 4: KEEPTOP mice
for bINCLUDE_ALL_MICE = [1 0]
    figure(4 - bINCLUDE_ALL_MICE); clf;
    
    if bINCLUDE_ALL_MICE
        TMPA = squeeze(mean(ALLFIB(TRIALS,:,a), 1, 'omitnan'))';  %average over trials
        TMPB = squeeze(mean(ALLFIB(TRIALS,:,b), 1, 'omitnan'))';  %average over 1st 4 trials
        STAT = ALL;
        
    else
        %only top #
        TMPA = squeeze(mean(ALLFIB(TRIALS,:,a(1:KEEPTOP)), 1, 'omitnan'))';  %average over trials
        TMPB = squeeze(mean(ALLFIB(TRIALS,:,b(1:KEEPTOP)), 1, 'omitnan'))';  %average over 1st trials
        STAT = TOP;
    end
      
    for k=1:length(AnalysisWindows)
        subplot(1, length(AnalysisWindows)+2, k); 
        
        if bINCLUDE_ALL_MICE
            %ALL mice
            CTL = CTL_GBZ{k}(:,1);
            GBZ = CTL_GBZ{k}(:,2);
        else
            %TOP# mice are the ones we keep
            CTL = CTL_GBZ{k}(1:KEEPTOP,1);
            GBZ = CTL_GBZ{k}(1:KEEPTOP,2);
        end
    
        axis([-1 3 -10 20]); hold on;
        xbee = beeswarm([CTL*0; GBZ*0+1], [CTL; GBZ], 'use_current_axes', 1, 'dot_size', DS*1.2, 'sort_style', 'fan'); 
        cla;
        plot([-1 3],[0 0],    ':k'); hold on;
        plot([-1 3],[0 0]+10, ':k'); hold on;
        plot(xbee(1:length(CTL)), CTL, 'o', 'markersize', DS, 'MarkerFaceColor',     [0 0 0]*0.9 + [1 1 1]*0.1, 'MarkerEdgeColor', 'none', 'MarkerEdgeColor', 'none'); hold on;
        plot(xbee(length(CTL)+1:end), GBZ, 'o', 'markersize', DS, 'MarkerFaceColor', [1 0 0]*0.9 + [1 1 1]*0.1, 'MarkerEdgeColor', 'none', 'MarkerEdgeColor', 'none'); hold on;
        axis off;
                
        for j=1:length(GraphicMode)
            if abs(GraphicMode(j))==1
                %Fit2m1 mode -- for both pdf and 95CI
                [N,E] = histcounts(STAT.Fit2m1(k).TermDiffBoot, 30, 'Normalization', 'pdf');
                C = (E(1:end-1) + E(2:end))/2; %bin centers
                C = C + STAT.Fit1(k).Term; %offset relative to control
                
                CI95 = STAT.Fit1(k).Term + STAT.Fit2m1(k).TermDiffCI95;
                CImid = STAT.Fit2(k).Term;
                
                XPDF  = 1.8; %x-position of PDF distribution               
                X95CI = 1.8; %x-position for CI
                clr1 = 'r';
                clr2 = 'k';
                
            elseif abs(GraphicMode(j))==2
                %Hybrid Mode -- PDF of Fit2m1, with 95CI of Mixing Test
                
                %pdf of Fit2m1
                [N,E] = histcounts(STAT.Fit2m1(k).TermDiffBoot, 30, 'Normalization', 'pdf');
                C = (E(1:end-1) + E(2:end))/2; %bin centers
                C = C + STAT.Fit1(k).Term; %offset relative to control
                
                %95 CI from Mixing Test
                CI95 = STAT.Fit2(k).Term + STAT.MixingTest(k).TermDiffNullCI95;
                CImid = STAT.Fit2(k).Term;
                
                X95CI = 1.8; %x-position for CI
                XPDF  = 2.0; %x-position of PDF distribution               
                clr1 = 'r';
                clr2 = 'k';
                
            elseif abs(GraphicMode(j))==3
                %MixingTest -- for both pdf and 95CI
                [N,E] = histcounts(STAT.MixingTest(k).TermDiffNullBoot, 31, 'Normalization', 'pdf');
                C = (E(1:end-1) + E(2:end))/2; %bin centers
                C = C + STAT.Fit2(k).Term; %offset relative to positive
                CI95 = STAT.Fit2(k).Term + STAT.MixingTest(k).TermDiffNullCI95;
                CImid = STAT.Fit2(k).Term;
                
                X95CI = 1.8; %x-position for CI
                XPDF  = 1.8; %x-position of PDF distribution               
                clr1 = 'r';
                clr2 = 'k';
            end
            %add points to start and end of histogram
            dC = mean(diff(C));
            C = [C(1)-dC C  C(end)+dC];
            N = [0 N 0];
            N = 1.5*0.3*N/max(N);
            h=patch(N + XPDF+(j-1)/1.5, C, clr1);
            set(h,'LineStyle', 'none', 'FaceAlpha', 1');
            
            %95% confidence interval and line
            plot([0 0]+X95CI+(j-1)/1.5,  CI95, clr2, 'linewidth', 5);
            %plot(X95CI+(j-1)/2, CImid,   'o', 'markersize', DS, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerEdgeColor', 'none');
        end
        plot([1 X95CI], [0 0]+STAT.Fit2(k).Term, '-', 'color', 'k', 'linewidth', 1);
        plot([0 XPDF+0.3], [0 0]+STAT.Fit1(k).Term, '-', 'color', 'k', 'linewidth', 1);
        grid on;
        
        %text on figure
        text(-0.8,18, [...
            'Time Window =  ' num2str(AnalysisWindows{k})  newline ...
            'Permutation Test:'  newline ...
            '     P = ' num2str(STAT.MixingTest(k).P_TermDiffVsNull) newline ...
            '     Null 95CI = (' num2str(STAT.MixingTest(k).TermDiffNullCI95') ')' newline ...
            'Effect Size: Mean Diff = ' num2str(STAT.Fit2m1(k).TermDiff) newline ...
            '     Boot 95CI = (' num2str(STAT.Fit2m1(k).TermDiffCI95') ')' newline  ...
            ], 'FontName', 'Courier');
    end
    
    %Waveform Plot
    subplot(1, length(AnalysisWindows)+2, length(AnalysisWindows)+[1:2]);
    plot(TIME, TIME*0, 'k'); hold on;
    
    plot(TIME, mean(TMPA),'k', 'linewidth', 3); hold on;  %average over mice
    SEM = std(TMPA)/sqrt(size(TMPA,1));  %SEM over mice
    h=patch([TIME fliplr(TIME)], [mean(TMPA)-SEM  fliplr(mean(TMPA)+SEM)], 'k');
    set(h,'LineStyle', 'none', 'FaceAlpha', 0.2');
    xlim([TIME(1) TIME(end)]);
    
    plot(TIME, mean(TMPB),'r', 'linewidth', 3); hold on;%average over mice
    SEM = std(TMPB)/sqrt(size(TMPB,1));%SEM over mice
    h=patch([TIME fliplr(TIME)], [mean(TMPB)-SEM  fliplr(mean(TMPB)+SEM)], 'r');
    set(h,'LineStyle', 'none', 'FaceAlpha', 0.2');
    grid on;
    ylim([-0.1 0.2]*100);
    title(['Day=' num2str(Day) ', Type=' num2str(Type), ', Trials=' num2str(TRIALS)]);

end

%print data into the command window 
ExperimentalStage
for k=1:length(AnalysisWindows)
    W = AnalysisWindows{k};
    disp(['Analysis Window #' num2str(k) ':  '  num2str(W)]);
    
    fprintf(['(ALL) Perm_P = ' num2str(ALL.MixingTest(k).P_TermDiffVsNull,'%.10f') '\n']);
    fprintf(['(TOP' num2str(KEEPTOP) ') Perm_P = ' num2str(TOP.MixingTest(k).P_TermDiffVsNull,'%.10f') '\n']);
    
    fprintf('%+.10f\t%+.10f\n', CTL_GBZ{k}');    
end
