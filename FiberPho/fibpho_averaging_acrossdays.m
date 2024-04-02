%Written 03/2023 by SB, with edits 12/2023, and input from MRT
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%Modified from https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/fiber-photometry-epoch-averaging-example/

%Requires TDT functions that can be downloaded at: https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/getting-started/

%Function to extract trial timepoint 415nm and 465nm signals across behavior days for one
%mouse, and to calculate the time-locked dF/F
%To run, be in folder with one animal's fiber pho data from multiple days in individual subfolders
%each day's folder needs to include the trialData variable from that day

%inputs:
%       Cohort: cohort string
%       Mouse: mouse number string
%       Bank: string "upper" for top 415A/465A bank, "lower" for bottom 415C/465C
%       bank, of the RZ10X
%
%outputs:
%       fibdat: structure that contains 1 row per day of analysis. Within
%       each day, there are multiple structures relative to each trial type
%       of interest (ex probe, Tone A rewarded, Tone A extinction, etc)
%       NOTE: this fibdat includes the raw 415/465 and the calculated dF/F for each trial for ease of
%       subsequent analysis

function fibpho_averaging_acrossdays(Cohort, Mouse, Bank)

path = cd;

%load the data
folderDir = dir(path);
analysisFiles = [];
for q = 1:size(folderDir)
    if contains(folderDir(q).name, Mouse) && folderDir(q).isdir
        subDirPath = [path filesep folderDir(q).name filesep];
        analysisFiles = [analysisFiles; subDirPath];
    end
end

%extract data from the correct Bank (upper or lower)
if isequal(Bank, 'upper')
    REF_EPOC = 'U11_'; % event store name. This holds solenoid activation times
    STREAM_STORE1 = 'x415A'; % name of the 415 store
    STREAM_STORE2 = 'x465A'; % name of the 465 store
    solval = 1; % the set solval by the User Input GUI in Synapse
elseif isequal (Bank, 'lower')
    REF_EPOC = 'U12_'; % event store name. This holds solenoid activation times
    STREAM_STORE1 = 'x415C'; % name of the 415 store
    STREAM_STORE2 = 'x465C'; % name of the 465 store
    solval = 2; % the set solval by the User Input GUI in Synapse
else
    error('not sure which bank to use');
end

fibdat = []; %initialize variable to store data

TRANGE = [-5 10]; % window size [start time relative to epoc onset, window duration]
BASELINE_PER = [-5 -1.5]; % baseline period within our window

for n=1:size(analysisFiles) 
    path = analysisFiles(n, :);
    data = TDTbin2mat(path, 'TYPE', {'epocs', 'scalars', 'streams'}); %load TDT data

    %load the trialData
    subFolderDir = dir(path);
    trialFiles = [];
    for q = 1:size(subFolderDir)
        if contains(subFolderDir(q).name, [Cohort '_' Mouse]) && contains(subFolderDir(q).name, 'trialData')
           subDirPath = [path subFolderDir(q).name];
           trialFiles = [trialFiles; subDirPath];
        end
    end
    
    if size(trialFiles, 1) == 1
        load(trialFiles(1, :));
    else
        error('too many trial files');
    end
    
    Values = []; %initialize temporary variable

    %calculate the time offset between matlab solenoid and tdt solenoid to
    %align
    rewTrials = trialData.trialInfo(:, trialData.trialInfo(3, :) > 0); %cue start times as recorded by the behavior code
    solTimes = (rewTrials(6, :) + 1.5)'; %add 1.5 seconds to get the solenoid activation times
    issueWithTrialDat = 1;
    try
        %calculate the offset for each trial
        if Bank == 'upper'
            offset = data.epocs.U11_.onset - solTimes;
        elseif Bank == 'lower'
            offset = data.epocs.U12_.onset - solTimes;
        end
        meanOffset = mean(offset); %calculate the mean of all the offsets
        stdOffset = std(offset); %std can be used to indicate if there is an issue - ex, offsets are very different trial to trial
        
        issueWithTrialDat = 0;      
        %use the calculated offset to get the fiber pho time points for all
        %trials (both rewarded and unrewarded)
        allTrialTimes = trialData.trialInfo(6, :) + 1.5 + meanOffset;      
    catch
        if issueWithTrialDat %this is generally if the TDT onset time points and trial number aren't the same size
            fprintf("Couldn't match fiber indices to trial data");
            if Bank == 'upper'
                Values(1).timeStamps = data.epocs.U11_.onset';
            elseif Bank == 'lower'
                Values(1).timeStamps = data.epocs.U12_.onset';
            end
            %if this issue arises, can still calculate the fibpho signals
            %around the solenoid times for the rewarded trials, but won't
            %be able to align to the behavior to get time points of other
            %unrewarded trials
            allTrialTimes = Values(1).timeStamps;
            Values(1).timeStamps = logical(ones(size(allTrialTimes)));
            Values(1).assetName = 'Tone A Training';
        end
    end
        
    %raw data all trials
    data_trials = TDTfilter_timestamp(data, 'VALUES', allTrialTimes', 'TIME', TRANGE); %filter the data around the time stamps for all trials (rewarded and unrewarded)
    
    %"Applying a time filter to a uniformly sampled signal means that the length of each segment could vary by one sample. 
    %Let's find the minimum length so we can trim the excess off before calculating the mean."
    minLength1 = min(cellfun('prodofsize', data_trials.streams.(STREAM_STORE1).filtered));
    minLength2 = min(cellfun('prodofsize', data_trials.streams.(STREAM_STORE2).filtered));
    data_trials.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data_trials.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
    data_trials.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data_trials.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);
    
    % downsample 10x and average 415 signal
    allSignals = cell2mat(data_trials.streams.(STREAM_STORE1).filtered');
    N = 10;
    F415_all = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F415_all(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength1 = size(F415_all,2);
    
    % downsample 10x and average 465 signal
    allSignals = cell2mat(data_trials.streams.(STREAM_STORE2).filtered');
    F465_all = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F465_all(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength2 = size(F465_all,2);
    
    %%%%%%%%dFF calculation%%%%%%%%%%%%%
    ts2 = TRANGE(1) + (1:minLength2) / data_trials.streams.(STREAM_STORE2).fs*N;
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    F0_raw = mean(F465_all(:, ind), 2); %F0 is calculated from the 3.5s prior to cue initiation per trial 
    x = 1:size(F0_raw);
    
    %calculate a double term exponential fit to all baseline trial F0s
    %this enables a less noisy approximation of F0 changes from bleaching
    %across the 1hr recording sessions
    ft = fittype('c+a1*exp(-x/tau1)+a2*exp(-x/tau2)');
    c = F0_raw(end);
    a = F0_raw(1)-F0_raw(end);
    tau1 = 10;
    tau2 = 100;
    ftoptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower',[0,0,0,0,0],'Upper',[Inf,Inf,Inf,Inf,Inf],'StartPoint',[a, a, c, tau1, tau2]);
    E = fit(x', F0_raw, ft, ftoptions);
    F0 = feval(E, x); %F0 per trial is the estimated F0 from the exponential fit
    err = immse(F0, F0_raw);
    
    %calculate (F-F0)/F0
    F465_dFF_all = (F465_all-F0)./F0;
    
    %save data
    fibdat(n).F0_raw = F0_raw;
    fibdat(n).F0 = F0;
    fibdat(n).F0_err = err;
        
    %save trial time stamps, split for each trial type for that behavior
    %session
    if ~issueWithTrialDat
        %training days
        if trialData.dayType < 2
            Values(1).timeStamps = trialData.trialInfo(2, :)==2500 & trialData.trialInfo(3, :)>0;
            Values(1).assetName = 'Training Tone A';
            %probes of unexpected tone B, unrewarded
            if trialData.numProbeTrials > 0
                Values(2).timeStamps = trialData.trialInfo(8,:) >0;
                Values(2).assetName = 'Probe trials';
            end
            if isfield(trialData, 'numOmissionTrials')   
                if trialData.numOmissionTrials > 0 %probes of learned tone A, unexpectedly unrewarded
                    Values(3).timeStamps = trialData.trialInfo(9,:) >0;
                    Values(3).assetName = 'Omission trials';
                end
            end
        %reward contingency flip session (testing session 1)
        elseif trialData.dayType == 2
            %baseline before the rule change (tone A still rewarded)
            Values(1).timeStamps = trialData.trialInfo(2, :)==2500 & trialData.trialInfo(3, :)>0;
            Values(1).assetName = 'Tone A preflip';
            %extinction
            Values(2).timeStamps = trialData.trialInfo(2, :)==2500 & trialData.trialInfo(3, :)==0;
            Values(2).assetName = 'Tone A post flip extinction';
            %conditioning
            Values(3).timeStamps = trialData.trialInfo(2, :)==11000;
            Values(3).assetName = 'Tone B post flip conditioning';
        %testing session 2    
        elseif trialData.dayType == 3
            %extinction
            Values(1).timeStamps = trialData.trialInfo(2, :)==2500;
            Values(1).assetName = 'Tone A test day extinction';
            %conditioning
            Values(2).timeStamps = trialData.trialInfo(2, :)==11000;
            Values(2).assetName = 'Tone B test day conditioning';
        end
    end
    
    for r=1:size(Values, 2) 
        %extract only the trials we care about        
        F415 = F415_all(Values(r).timeStamps', :);
        F465 = F465_all(Values(r).timeStamps', :);
        F465_dFF = F465_dFF_all(Values(r).timeStamps', :);

        minLength1 = size(F415,2);
        minLength2 = size(F465,2);
        
        %--Subplot 1: average and plot the raw signals--
        % Create mean signal, standard error of signal, and DC offset of 415 signal
        meanSignal1 = mean(F415);
        stdSignal1 = std(double(F415))/sqrt(size(F415,1));
        dcSignal1 = mean(meanSignal1);
                
        % Create mean signal, standard error of signal, and DC offset of 465 signal
        meanSignal2 = mean(F465);
        stdSignal2 = std(double(F465))/sqrt(size(F465,1));
        dcSignal2 = mean(meanSignal2);
        
        % Create the time vector for each stream store
        ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
        ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
      
        %save the raw data and timestamps
        fibdat(n).trialType(r).F415 = F415;
        fibdat(n).trialType(r).F465 = F465;
        fibdat(n).trialType(r).timeStamps = Values(r).timeStamps;
        fibdat(n).trialType(r).assetName = Values(r).assetName;
        fibdat(n).ts = ts2;

        % Subtract DC offset to get signals on top of one another
        meanSignal1 = meanSignal1 - dcSignal1;
        meanSignal2 = meanSignal2 - dcSignal2;        
        
        %Plot
        figure('WindowStyle', 'docked');
        subplot(2,1,1)
        plot(ts1, meanSignal1, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 3); hold on;
        plot(ts2, meanSignal2, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3);
        % Plot vertical line at solenoid onset, time = 0
        ylimit = get(gca, 'ylim');
        line([0 0], ylimit, 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)
        % Plot vertical line at tone onset, time = -1.5
        line([-1.5 -1.5], ylimit, 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)
        % Make a legend
        legend('415 nm','465 nm','Solenoid onset', 'AutoUpdate', 'off');
        % Create the standard error bands for the 415 signal
        XX = [ts1, fliplr(ts1)];
        YY = [meanSignal1 + stdSignal1, fliplr(meanSignal1 - stdSignal1)];
        % Plot filled standard error bands.
        h = fill(XX, YY, 'g');
        set(h, 'facealpha',.25,'edgecolor','none')
        % Repeat for 465
        XX = [ts2, fliplr(ts2)];
        YY = [meanSignal2 + stdSignal2, fliplr(meanSignal2 - stdSignal2)];
        h = fill(XX, YY, 'r');
        set(h, 'facealpha',.25,'edgecolor','none')
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('mV', 'FontSize', 12)
        title(sprintf([Cohort Mouse ' Day %d ' Values(r).assetName ' Response, %d Trials'], n, size(F465, 1)))      
        
        %--Subplot 2: average and plot the dF/F signals--
        fibdat(n).trialType(r).F465_dFF = F465_dFF;

        % Create mean signal, standard error of signal, and DC offset of dFF signal
        meanSignal3 = mean(F465_dFF);
        stdSignal3 = std(double(F465_dFF))/sqrt(size(F465_dFF,1));
        
        %Plot
        subplot(2,1,2)
        plot(ts2, meanSignal3, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
        % Plot vertical line at solenoid onset, time = 0
        ylimit = get(gca, 'ylim');
        line([0 0], ylimit, 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)
        % Plot vertical line at tone onset, time = -1.5
        line([-1.5 -1.5], ylimit, 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)
        % Make a legend
        legend('465 dF/F','Solenoid onset', 'AutoUpdate', 'off');
        % Create the standard error bands for the signal
        XX = [ts2, fliplr(ts2)];
        YY = [meanSignal3 + stdSignal3, fliplr(meanSignal3 - stdSignal3)];
        h = fill(XX, YY, 'r');
        set(h, 'facealpha',.25,'edgecolor','none')
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('dF/F', 'FontSize', 12)
        title(sprintf([Cohort Mouse ' Day %d ' Values(r).assetName ' Response, %d Trials'], n, size(F465_dFF, 1)))
    end
end

%save data
path = cd;
filename = [path filesep Cohort '_' Mouse '_allTrials_fibdat_dFF.mat'];
save(filename, 'fibdat');


end
