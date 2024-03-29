%Written in 2018 by SB and MRT, with edits 2019, 2020, 2021, and 2023
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%Instructions:
%   go to cohort data folder
%   Inputs:
%       string for cohort (ex, 'VK'
%       string for mouse = 'Mouse 1')
%       recalculate = 1 for yes, 2 for this day only, 0 for no
%   Outputs:
%       schultz: data structure containing behavior values
%       various plots of the analyzed data

function OUT = licksRelativeAnalysisDualCondExt_NI(cohort, mouse, recalculate)

currentFolder = pwd;

%filter value for later plotting
filtNum = 30;

if nargin<3 || isempty(recalculate)
    recalculate = 0;
end

%array of time stamps per trial we are interested in for anticipatory licking (from -3 to 3 with a step size dt)
dt=0.1;
t = -3: dt : 3;
%array of time stamps per trial we are interested in for detecting first lick post-reward (from -0 to 5 with a step size 0.01)
collectt = 0:0.01:5;
%time windows of interest
tSchultz = [-1.5 0]; %key analysis window for anticipatory licking detection
tMotiv = [0 1.5]; %analysis window for reward collection window detection


%if some days of behavior have already been analyzed, don't recalculate
fullDirPath = [currentFolder filesep mouse filesep '*analysis_tSchultz_' num2str(tSchultz(1)) '-' num2str(tSchultz(2)) 's.mat'];
checkDir = dir(fullDirPath);
if ~isempty(checkDir) && recalculate ~= 1
    load([checkDir.folder filesep checkDir.name]);
    %get the number of days of analysis
    dirUnsort = dir([currentFolder filesep mouse filesep 'Data/*.mat']);
    numFiles = size(dirUnsort, 1);
    daysAnalyzed = length(unique(schultz(1).sessionIndx)) - 1;
end
%for the first time running the analysis, or if the recalculate flag is set
%to 1
if recalculate == 1
    startFile = 1;
    %initialize the variables holding the analysis
    for k=1:2
        %k = 1 is the initial lower frequency tone - tone B, conditioning, starts as 0kHz
        %k = 2 is the initial higher frequency tone - tone A, extinction, 2.5kHz
        schultz(k).Tones = 0;
        schultz(k).performance = []; %Keeps track of licking amount over the tSchultz time window for each trial
        schultz(k).motivation = []; %Keeps track of licking amount over the tMotiv time window for each trial
        schultz(k).firstLick = []; %keeps track of time from cue initiation to first lick per trial
        schultz(k).sessionIndx = []; %keeps track of the trial at which a new day of behavior has started
        schultz(k).flip = 0; %keeps track of the two trials (one per k) at which the reward contingency flip occurs
        schultz(k).dayIndex = []; %keeps track of the day of behavior
        schultz(k).probeIndx = []; %keeps track of indices of the tone B probe trials
        schultz(k).anticipatorySpeed = []; %keeps track of rotary encoder rotations over the tSchultz time window
        schultz(k).collectionSpeed = []; %keeps track of rotary encoder rotations over the tMotiv time window
        schultz(k).fullSpeed = []; %keeps track of rotary encoder rotations over the full t time window
    end
    %if adding days of behavior to the analysis, mark how many days have been analyzed already
elseif recalculate == 2
    startFile = daysAnalyzed + 1;
end

%variable to hold speed, not trial dependent, across the whole behavior session
totalSpeed = [];

%if need to analyze any number of days
if recalculate ~= 0
    %load the trial data
    fullDirPath = [currentFolder filesep mouse filesep 'Data/*.mat'];
    dirUnsort = dir(fullDirPath);

    %make sure data is in the correct order by date
    for q = 1:size(dirUnsort, 1)
        monthday = NumbersInString(dirUnsort(q).name);
        dirUnsort(q).month = abs(monthday(3));
        dirUnsort(q).day = abs(monthday(4));
    end
    [B, I]=sort([dirUnsort.day]);
    %sort by day first
    dirUnsort = dirUnsort(I);
    %then by month
    [B, I]=sort([dirUnsort.month]);
    D = dirUnsort(I);

    numFiles = size(D, 1);

    %figure 1 is built as we go through
    figure(1); clf; %clear figure at beginning
    numBaselines = 0;

    %----------------GO THROUGH DAY BY DAY-------------------------------%
    for i = startFile:numFiles
        load([D(i).folder filesep D(i).name]);  %trialData variable saved from behavior code

        %keep track of session boundaries in the trials
        %note where the start of this day of trials is
        for k=1:length(schultz)
            schultz(k).sessionIndx(end+1) = length(schultz(k).performance);
        end

        dayType = trialData.dayType;
        toneDurationSec = 1.5;


        %--------------EXTRACT LICK, SUCROSE, AND TONE TIMES------------%
        %get indices of when a lick occurred
        lickTimesIdx = trialData.lickValues(2,:)>0;
        %get the actual time value matching the indices
        lickTimes = trialData.recordTimes(1, lickTimesIdx);



        %Get sucrose delivery times and times of silent/unrewarded trials
        sucroseStartTimes = zeros(3, size(trialData.trialInfo, 2));
        sucroseStartTimes(2, :) = trialData.trialInfo(2, :); %second row is the tone played

        %for rewarded trials, get the real time of solenoid opening from the NI card event.TimeStamps
        actualSucroseIndices = find(trialData.sucroseDeliveryValues > 1.5);
        actualSucrose = zeros(2, trialData.numRewardedTrials);
        %keep the first index of sucrose delivery
        actualSucrose(1, 1) = actualSucroseIndices(1);
        countThroughSucroseDeliv = 2;
        %grab first event.TimeStamp ms that the solenoid opens for every rewarded trial
        for k = 2:size(actualSucroseIndices, 2)
            if actualSucroseIndices(k) - actualSucroseIndices(k-1) > 1
                actualSucrose(1, countThroughSucroseDeliv) = actualSucroseIndices(k);
                countThroughSucroseDeliv = countThroughSucroseDeliv + 1;
            end
        end
        %translate indices into time stamps
        actualSucrose(2, :) = trialData.recordTimes(1, actualSucrose(1, :));

        %use the NI card downward blip at end of a stimulus session to get
        %the "sucrose delivery time" of unrewarded trials
        %(more accurate than the Matlab clock)
        blipIndices = find(trialData.sucroseDeliveryValues < -2);
        %get rid of any artifacts from the start of recording
        blipIndices = blipIndices(blipIndices > 100);
        actualBlip = zeros(2, trialData.numTrials);
        actualBlip(1, 1) = blipIndices(1);
        countthruBlip = 2;
        for k = 2:size(blipIndices, 2)
            if blipIndices(k) - blipIndices(k-1) > 1
                actualBlip(1, countthruBlip) = blipIndices(k);
                countthruBlip = countthruBlip + 1;
            end
        end
        %translate indices into time stamps
        actualBlip(2, :) = trialData.recordTimes(1, actualBlip(1, :));
        %subtract the two seconds of stimulus and the solenoid open time to
        %get the initiation time of no reward/reward omission
        actualBlip(3, :) = actualBlip(2, :) - 2 - (trialData.solenoidOpenDuration)/1000;

        %combine recorded actual and calculated non-actual solenoid opening times into one variable
        countThroughRewTrials = 1;
        for k = 1:size(trialData.trialInfo, 2)
            %for rewarded trials, use the actual time of solenoid opening
            if trialData.trialInfo(3, k) > 0
                sucroseStartTimes(1, k) = actualSucrose(2, countThroughRewTrials);
                countThroughRewTrials = countThroughRewTrials + 1;
                sucroseStartTimes(3, k) = 1; %note rewarded or not
                %unrewarded trials, use the estimate from the end of trial blips
            else
                sucroseStartTimes(1, k) = actualBlip(3, k);
            end
        end

        %get tone times - no reliable way to extract from NI trace, so calculate from solenoid time
        toneStartTimes = sucroseStartTimes(1, :) - trialData.rewardDelay - toneDurationSec;

        %Extract the tone information
        Tones = unique(trialData.trialInfo(2,:));
        if dayType == 1
            %during training, the tone that isn't zero is the conditioning tone
            baselineRewTone = max(Tones);
            schultz(2).Tones = baselineRewTone;
        elseif dayType > 1 && dayType < 3
            %on the flip day, we now have 3 tones in play: 0, prevRewTone, newRewTone
            baselineRewTone = schultz(2).Tones;
            newRewTone = Tones(find(Tones ~= baselineRewTone));
            newRewTone = newRewTone(find(newRewTone ~= 0));
            %since "unique" puts the tones in order, make sure the new rew tone is in the same place as zero was
            Tones = [newRewTone, baselineRewTone];
            %on the flip day, switch all the zeroes to the new rew tone (to make things easier)
            for k=1:size(trialData.trialInfo,2)
                if sucroseStartTimes(2, k) == 0
                    sucroseStartTimes(2, k) = newRewTone;
                end
            end
            schultz(1).Tones = newRewTone;

            %if a flip day, store the trial at which the flip occurs
            tmpflip(1) = 0;
            tmpflip(2) = 0;
            for k = 1:trialData.flipTrial - 1
                sIdx = find(Tones == sucroseStartTimes(2,k));
                if dayType == 2
                    schultz(sIdx).flip = schultz(sIdx).flip + 1;
                end
            end
            %now add to all the trials from the previous days
            for k=1:2
                if dayType == 2
                    schultz(k).flip = schultz(k).flip + length(schultz(k).performance);
                end
            end
        elseif dayType == 3
            %on testing days, ensure that the order of the tones is still correct (otherwise data will be flipped)
            baselineRewTone = schultz(2).Tones;
            newRewTone = Tones(find(Tones ~= baselineRewTone & Tones ~= 0));
            %since "unique" puts the tones in order, make sure the new rew tone is in the same place as zero was
            Tones = [newRewTone, baselineRewTone];
        end

        %------------GO THROUGH ALL TRIALS ON THIS DAY--------------------%

        %----------------------LICK CALCULATIONS---------------------------%
        L = zeros(length(Tones), length(t)); %to hold lick data
        
        %go through each trial
        for k=1:size(sucroseStartTimes, 2)
            SDT = sucroseStartTimes(1, k);  %this is the time, in sec, of the reward delivery around which we are analyzing licks
            RelLickTime = lickTimes - SDT; %time of licks relative to sucrose delivery
            I =  (RelLickTime>=t(1)) &  (RelLickTime<=t(end)); %licks in our t window of interest
            nT  = find(Tones == sucroseStartTimes(2,k)); %which tone

            %--All licks--
            %go through our "t" variable, and ask at each time if the IR beam was broken or a touch occurred
            Ltmp = t*0;
            for j=1:length(t)
                dtToNearestLick = min( abs(t(j) - RelLickTime(I)));
                if (dtToNearestLick <= dt/2)
                    L(nT,j) = L(nT,j) + 1;
                    Ltmp(j) = Ltmp(j) + 1;  %sum for just this one
                end
            end

            %--Time to first lick--
            %go through our "collectt" variable (only after reward delivery), to find the first lick after reward delivery
            Ltmp2 = collectt*0;
            firstlick = 0;
            for j=1:length(collectt)
                dtToNearestLick = min( abs(collectt(j) - RelLickTime(I)));
                if (dtToNearestLick <= 0.01/2)
                    Ltmp2(j) = Ltmp2(j) + 1;  %sum for just this one
                    if ~firstlick
                        schultz(nT).firstLick(end+1) = collectt(j); %the time of the first lick relative to sucrose delivery
                        firstlick = 1;
                    end
                end
            end

            %if no licks during the window, set time to first lick as 5.1s
            %(end of the time window we are checking)
            if firstlick == 0
                schultz(nT).firstLick(end+1) = 5.1;
            end

            %--Anticipatory licking calculation--
            %predetermined anticipatory licking window
            I =  (t>=tSchultz(1)) &  (t<=tSchultz(2));
            perf = sum(Ltmp(I))/sum(I); %for I being the anticipatory licking window, divide number of licks in that window by number of "checks" that is the window
            schultz(nT).performance(1, end+1) = perf; %row 1 is the licking performance
            schultz(nT).performance(2, end) = trialData.trialInfo(6, k); %row 2 is the recorded NI time when the tone started
            Iafter = (collectt>=tMotiv(1)) & (collectt<=tMotiv(2)); %assess motivation/reward collection
            schultz(nT).motivation(end+1) = sum(Ltmp2(Iafter))/sum(Iafter);%divide number of licks in the window by number of "checks" that is the window

            %if a probe trial, note the index
            if size(trialData.trialInfo, 1) > 7
                if trialData.trialInfo(8, k) > 0 && nT == 1 %if a probe trial
                    schultz(nT).probeIndx = [schultz(nT).probeIndx length(schultz(nT).performance)];
                    %nearby tone A trials for reference
                    schultz(2).probeIndx = [schultz(2).probeIndx length(schultz(2).performance)];
                end
            end

        end


        %----------------ROTATIONAL SPEED CALCULATIONS--------------------%
        %--Total speed across the session--
        %make an array with both time and the rpm
        rotationalSpeed = [trialData.recordTimes(1, :); trialData.rotationalSpeed];
        %save all rotational speed data for last three training days, flip,
        %and test days only (as these are the only days we care about)
        if i>(numFiles - 5)
            totalSpeed = [totalSpeed rotationalSpeed];
        end
        
        %--Speed during trial timepoints--
        %temp variables to hold the speed information across all trials
        speedSumTone1 = t*0;
        speedSumTone2 = t*0;
        speedSumTone1 = speedSumTone1(1:end-1);
        speedSumTone2 = speedSumTone2(1:end-1);

        %go through each trial
        for k = 1:size(trialData.trialInfo, 2)
            sucDelivTime = sucroseStartTimes(1, k); %this is the time, in sec, of the reward delivery around which we are analyzing RPM
            nT  = find(Tones == sucroseStartTimes(2,k)); %which tone
            relativeSpeedTime = [rotationalSpeed(1, :) - sucDelivTime; rotationalSpeed(2, :)]; %subtract to get the t window

            relevantSpeeds = [];
            relevantSpeeds(1, :) = relativeSpeedTime(1, (relativeSpeedTime(1, :)>=t(1)) & (relativeSpeedTime(1, :)<=t(end)));%get the relevant times
            relevantSpeeds(2, :) = relativeSpeedTime(2, (relativeSpeedTime(1, :)>=t(1)) & (relativeSpeedTime(1, :)<=t(end)));%get the rotational speed for the relevant times

            %temporary variable to hold the speeds per trial
            tempSpeedSumTone1 = t*0;
            tempSpeedSumTone2 = t*0;
            for j = 1:length(t)
                %get the speed at the time closest to j
                dtToNearestCheckTime = min( abs(t(j) - relevantSpeeds(1, :)));
                nearestCheckTime = [(t(j) + dtToNearestCheckTime) (t(j) - dtToNearestCheckTime)];
                %have a tolerance to deal with floating point minute differences
                tol = eps('single');
                if ismember(nearestCheckTime(1), relevantSpeeds(1, :))
                    speedAtNearestCheckTime = relevantSpeeds(2, abs(relevantSpeeds(1, :) - nearestCheckTime(1)) <= tol);
                else
                    speedAtNearestCheckTime = relevantSpeeds(2, abs(relevantSpeeds(1, :) - nearestCheckTime(2)) <= tol);
                end

                %add that speed to the temporary sum for j, rewarded or unrewarded
                if nT == 1
                    tempSpeedSumTone1(j) = tempSpeedSumTone1(j) + speedAtNearestCheckTime(1);
                else
                    tempSpeedSumTone2(j) = tempSpeedSumTone2(j) + speedAtNearestCheckTime(1);
                end
            end

            %get the derivative - want changes in edge counts from the NI card
            %tempSpeedSum is the number of edge counts total
            %we want the difference, ie the number of rotations that happened during this time period
            tempDiffSpeed1 = diff(tempSpeedSumTone1);
            tempDiffSpeed2 = diff(tempSpeedSumTone2);

            %add to the summed RPM across all trials
            speedSumTone1 = speedSumTone1 + tempDiffSpeed1;
            speedSumTone2 = speedSumTone2 + tempDiffSpeed2;

            %add tempDiffSpeed to our larger data structures
            %longitudinal analysis of anticipation window
            antWindow = (t>=tSchultz(1)) &  (t<=tSchultz(2));
            %longitudinal analysis of reward collection window
            rewCollectWindow = (t>=tMotiv(1)) &  (t<=tMotiv(2));
            if nT == 1 %tone 1 variables
                schultz(nT).anticipatorySpeed(end+1) = sum(tempDiffSpeed1(antWindow)); 
                schultz(nT).collectionSpeed(end+1) = sum(tempDiffSpeed1(rewCollectWindow));
                schultz(nT).fullSpeed = [schultz(nT).fullSpeed; tempDiffSpeed1];
            else %tone 2 variables
                schultz(nT).anticipatorySpeed(end+1) = sum(tempDiffSpeed2(antWindow)); 
                schultz(nT).collectionSpeed(end+1) = sum(tempDiffSpeed2(rewCollectWindow));
                schultz(nT).fullSpeed = [schultz(nT).fullSpeed; tempDiffSpeed2];
            end
        end


        %--------------PLOT LICKS IN -3 to 3S WINDOW PER BEHAVIOR SESSION----------------------%
        style = {'-' '--'}; %tone 1, then tone 2

        L=squeeze(L);
        if size(L,2)==1
            L = L';  %make sure time varies with columns, not rows
        end

        %color scheme for plotting the day by day data
        figure(1);
        %training sessions - greyscale, with lighter colors being more recent days
        if dayType == 1
            numBaselines = numBaselines + 1;
            cmap = gray;
            cmap(end-10:end,:)=[];
            imap = round(size(cmap,1)*(numBaselines)/numFiles);
            clr = cmap(imap,:);
        %flip session
        elseif dayType == 2
            clr = 'b';
        %testing sessions
        else
            switch (numFiles - i)
                case 0
                    %last run
                    clr = 'r';
                case 1
                    %penultimate
                    clr = 'm';
            end
        end

        %plot figure 1, the time-locked data, to visualize licking within
        %the t window
        for q=1:size(L,1)
            %plot the non-normalized data
            subplot(1,2,2); hold on;
            plot(t, L(q,:), ['o' style{q}], 'linewidth', 2, 'color', clr);
            grid on;
            set(gca,'fontsize',20);

            %plot the normalized data (normalize to number of trials)
            subplot(1,2,1); hold on;
            if q == 1
                plot(t,L(1,:)./trialData.numUnrewardedTrials, ['o' style{q}], 'linewidth', 2, 'color', clr); 
            else
                plot(t,L(2,:)./trialData.numRewardedTrials, ['o' style{q}], 'linewidth', 2, 'color', clr); 
            end
            grid on;
            set(gca,'fontsize',20);
        end
        drawnow; %update graphics
    end

    %get end of last day trial index
    for k=1:length(schultz)
        schultz(k).sessionIndx(end+1) = length(schultz(k).performance);
    end

    %add details to figure 1
    figure(1);
    for k=1:2
        subplot(1,2,k);
        ylim = get(gca,'ylim');
        %add vertical line to the graph at zero for reward delivery
        line([0 0], ylim, 'color', 'k', 'linewidth', 2); hold on;
        %add vertical line to graph at -1.5 for tone initiation
        line([-1.5 -1.5], ylim, 'color', 'k', 'linewidth', 2); hold on;
        hold off;
        xlabel('Time (seconds)');
        if k==1
            ylabel('Normalized Licks');
            title([mouse ': licks relative to sucrose delivery']);
        elseif k == 2
            ylabel('Licks');
            title([mouse ': licks relative to sucrose delivery']);
            legend({['Tone ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone ' num2str(schultz(2).Tones/1000) 'kHz']});
        end
    end

    %save tshultz window
    for k=1:length(schultz)
        schultz(k).tSchultz = tSchultz;
    end

    %Save day by day calculations to save computing time later
    filename = [currentFolder filesep mouse filesep cohort '_' mouse '_analysis_tSchultztmp_' num2str(tSchultz(1)) '-' num2str(tSchultz(2)) 's.mat'];
    save(filename, 'schultz');
    filename2 = [currentFolder filesep mouse filesep cohort '_' mouse '_analysis_runningtmp.mat'];
    save(filename2, 'totalSpeed');
end


%set styles for all following graphs
style1 = '-';  %licks
style2 = '--'; %motivation
%set colors for following graphs
clr = ['r', 'c']; %tone B is red, tone A is cyan, tone C(ref) is black

%------------------------FIGURE 2: ALL TRIALS PER MOUSE------------------%
figure(2);
clf; hold on;
lines = zeros(2, 4);
for k=1:length(schultz)
    Raw = schultz(k).performance(1, :);
    RawSpeed = schultz(k).anticipatorySpeed;
    SessIndx = schultz(k).sessionIndx;
    flip = schultz(k).flip(1);
    probeIndx = schultz(k).probeIndx;
    SessIndx(SessIndx==0)=[];
    motiv = schultz(k).motivation;

    %subplot for anticipatory licking
    subplot(2, 1, 1); hold on;
    yyaxis left
    lines(k, 1) = plot(Raw,style1,'color',clr(k),'linewidth',2);
    ylabel('Anticipatory licking (solid lines)');
    yyaxis right
    lines(k, 5) = plot(motiv, style2, 'color', clr(k), 'linewidth', 2);
    ylabel('Post-reward licking (dashed lines)');
    lines(k, 3)= plot(SessIndx, Raw(SessIndx),'o', 'markersize', 20, 'color', 'k','linewidth', 2);
    %only plot flipday marker if we have had a flip day
    if ~schultz(k).flip == 0
        lines(k, 4)= plot(flip, Raw(flip), 'x', 'markersize', 20, 'color', 'k', 'linewidth', 2);
    end
    plot(probeIndx, Raw(probeIndx), '*', 'markersize', 20, 'color', 'k', 'linewidth', 2);

    %subplot for speed
    subplot(2, 1, 2); hold on;
    lines(k, 2) = plot(RawSpeed, style1, 'color', clr(k), 'linewidth', 2);
    lines(k, 3)= plot(SessIndx, RawSpeed(SessIndx),'o', 'markersize', 20, 'color', 'k','linewidth', 2);
    %only plot flipday marker if we have had a flip day
    if ~schultz(k).flip == 0
        lines(k, 4)= plot(flip, RawSpeed(flip), 'x', 'markersize', 20, 'color', 'k', 'linewidth', 2);
    end
    plot(probeIndx, RawSpeed(probeIndx), '*', 'markersize', 20, 'color', 'k', 'linewidth', 2);

end

%label figure 2
figure(2);
for k = 1:2
    subplot(2, 1, k);
    set(0,'DefaultAxesFontSize',15)
    xlabel('Trial number');
    if k == 1
        title('All Trials: Licking Behavior');
        ylabel('Licking');
        if ~schultz(k).flip == 0
            legend([lines(1, 1) lines(2, 1) lines(1, 3) lines(1, 4)], {['Tone A'], ['Tone B'], 'Day change', 'Flip trial'});
        else
            legend([lines(1, 1) lines(2, 1) lines(1, 3)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz'], 'Day change'});
        end
    elseif k == 2
        title('All Trials: Running Behavior');
        ylabel('Average running speed');
        legend([lines(1, 2) lines(2, 2)], {['Tone ' num2str(schultz(1).Tones/1000) 'kHz speed'], ['Tone ' num2str(schultz(2).Tones/1000) 'kHz speed']});
    end
end

%----------FIGURE 3: ALL TRIALS BY DAY WITH BREAKS-----------------------%
%get the day index variable
for d = 1:numFiles
    for k = 1:length(schultz)
        schultz(k).sessionIndx = unique(schultz(k).sessionIndx); %remove doubles
        trialIdx = [schultz(k).sessionIndx(d), schultz(k).sessionIndx(d+1)];
        dt = 1/((trialIdx(2)-trialIdx(1)));
        daySplit = d-1:dt:(d-dt);
        for i = 1:length(daySplit)
            schultz(k).dayIndex(end+1) = daySplit(i);
        end
    end
end
%create an index for the trials based on where they occur in time
for d = 1:numFiles
   for k = 1:length(schultz)
      schultz(k).sessionIndx = unique(schultz(k).sessionIndx); %remove doubles
      trialIdx = [schultz(k).sessionIndx(d)+1, schultz(k).sessionIndx(d+1)];
      schultz(k).performance(3, trialIdx(1):trialIdx(2)) = (schultz(k).performance(2, trialIdx(1):trialIdx(2))/3600)+(d-1);
   end
end

%create the filter
filt = designfilt('lowpassfir', 'FilterOrder', filtNum, 'CutoffFrequency', .01);

figure(3); clf; hold on;
lines = zeros(2, 4);
%generate the filter with day breaks
for k = 1:2
    flip = schultz(k).flip(1);
    schultz(k).filtful = [];
    xlims = [schultz(k).sessionIndx flip];
    xlims = sort(xlims);
    Raw = schultz(k).performance(1, :);
    %code to insert a break in the filter between days and at the flip
    for i = 1:length(xlims)-1
        tmp = Raw(xlims(i)+1:xlims(i+1));
        firstval = Raw(xlims(i)+1);
        tmp = tmp - firstval;

        tmpfilt = filter(filt, tmp);
        tmpfilt = tmpfilt + firstval;
        schultz(k).filtful = [schultz(k).filtful tmpfilt];
    end
end

%plot
for k=1:length(schultz)
    timeIndex = schultz(k).performance(3, :);
    flip = schultz(k).flip(1);

    lines(k, 1) = plot(timeIndex, schultz(k).filtful,style1,'color',clr(k),'linewidth',2);
    if ~schultz(k).flip == 0
        lines(k, 4) = plot(timeIndex(flip), schultz(k).filtful(flip), 'x', 'markersize', 20, 'color', 'k', 'linewidth', 2);
    end
    startidxes = schultz(k).sessionIndx(2:end-1) + 1;
    plot(schultz(k).filtful(startidxes), '.', 'markersize', 10, 'color', 'k', 'linewidth', 2);
end

%label the graph
set(gca,'fontsize',30);
title([cohort ' ' mouse ' All Trials: Licking Behavior Time Locked']);
xlabel('Session');
ylabel('Anticipatory licking');
if ~schultz(1).flip == 0
    legend([lines(1, 1) lines(2, 1) lines(1, 4)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz'], 'Flip trial'}, 'Location', 'northwest');
else
    legend([lines(1, 1) lines(2, 1)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz']}, 'Location', 'northwest');
end
xticks(0:1:length(schultz(1).sessionIndx)-1);


%----------FIGURE 4: TIME TO FIRST LICK------------------------------------
figure(4); clf; hold on;
lines = zeros(2, 4);
%plot
for k=1:length(schultz)
    Raw = schultz(k).firstLick;
    DayIndx = schultz(k).dayIndex;
    flip = schultz(k).flip(1);
    performanceSmooth = filter(filt, Raw);
    if k == 1 && schultz(k).flip > 0 %only plot tone B trials after flip
        lines(k, 1) = plot(DayIndx(flip:end), performanceSmooth(flip:end), 'color', clr(k), 'linewidth', 2);
    else
        lines(k, 1) = plot(DayIndx, performanceSmooth,style1,'color',clr(k),'linewidth',2);
    end
    if ~schultz(k).flip == 0
        lines(k, 4) = plot(DayIndx(flip), performanceSmooth(flip), 'x', 'markersize', 20, 'color', 'k', 'linewidth', 2);
    end
end
%label figure
set(gca,'fontsize',30);
title([cohort ' ' mouse ' All Trials: Time to first lick']);
xlabel('Session');
ylabel('Time to first lick');
if ~schultz(1).flip == 0
    legend([lines(1, 1) lines(2, 1) lines(1, 4)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz'], 'Flip trial'}, 'Location', 'northwest');
else
    legend([lines(1, 1) lines(2, 1)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz']}, 'Location', 'northwest');
end
xticks(0:1:length(schultz(1).sessionIndx)-1);
ax = gca;
ax.YDir = 'reverse';

%---------------FIGURE 5: PLOT PROBE TRIALS------------------------------------
%only plot if there has actually been probe trials
if schultz(1).probeIndx > 0
    figure(5); clf; hold on;
    for k=1:2
        plot(schultz(k).performance(1, schultz(k).probeIndx), 'color', clr(k), 'linewidth', 2);
    end
    title([cohort ' ' mouse ' Probe trials']);
    legend('Probe', 'Nearby Tone A');
    xlabel('Trial');
    ylabel('Anticipatory licking');
    meanprobe = mean(schultz(1).performance(1, schultz(1).probeIndx));
    meannearby = mean(schultz(2).performance(1, schultz(2).probeIndx));
    schultz(2).sessionIndx = unique(schultz(2).sessionIndx);
    meanprobeday = mean(schultz(2).performance(1, (schultz(2).sessionIndx(10)+1):schultz(2).sessionIndx(11)));
    fprintf('Mean probe response: %0.4f\n', meanprobe); %info for exclusion criteria
    fprintf('Mean nearby response: %0.4f\n', meannearby); %info for exclusion criteria
    fprintf('Mean full probe day ant licking to tone A response: %0.4f\n', meanprobeday); %info for exclusion criteria
end

%--------------------FIGURE 6: RUNNING SPEED BY DAY---------------------
%create the filter
filt = designfilt('lowpassfir', 'FilterOrder', filtNum, 'CutoffFrequency', .01);
figure(6); clf; hold on;
lines = zeros(2, 4);
%plot
for k=1:length(schultz)
    timeIndex = schultz(k).performance(3, :);
    flip = schultz(k).flip(1);
    speedSmooth = filter(filt, schultz(k).anticipatorySpeed(1, :));
    collectSpeedSmooth =  filter(filt, schultz(k).collectionSpeed(1, :));

    lines(k, 1) = plot(timeIndex, speedSmooth, style1,'color',clr(k),'linewidth',2);
    plot(timeIndex, collectSpeedSmooth, style2, 'color', clr(k), 'linewidth', 2);
    if ~schultz(k).flip == 0
        lines(k, 4) = plot(timeIndex(flip), speedSmooth(flip), 'x', 'markersize', 20, 'color', 'k', 'linewidth', 2);
        plot(timeIndex(flip), collectSpeedSmooth(flip), 'x', 'markersize', 20, 'color', 'k', 'linewidth', 2);
    end
end
%label 
set(gca,'fontsize',30);
title([cohort ' ' mouse ' All Trials: Treadmill Speed By Session']);
xlabel('Session');
ylabel('Number of rotations per trial');
if ~schultz(1).flip == 0
    legend([lines(1, 1) lines(2, 1) lines(1, 4)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz'], 'Flip trial'}, 'Location', 'northwest');
else
    legend([lines(1, 1) lines(2, 1)], {['Tone: ' num2str(schultz(1).Tones/1000) 'kHz'], ['Tone: ' num2str(schultz(2).Tones/1000) 'kHz']}, 'Location', 'northwest');
end
xticks(0:1:length(schultz(1).sessionIndx)-1);

%Pass outputs
OUT.schultz = schultz;
end