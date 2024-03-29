%Written 2023 by SB
%Cleaned up and commented 3/2024 by SB
% MIT License
% Copyright (c) 2024 Michael R. Tadross

%function to combine spiking data from different cells

%Instructions:
%   run in a folder with spiking data extracted from plotSpikes from multiple cells
%   (i.e., data structures that are one channel (one single row) saved from the Templates_data structure output by plotSpikes)
%   input: condition = descriptive string for saving
%   output: allCells data structure with calculated spiking summary data for each cell as a new row

function spikingAcrossCells(condition)

currentFolder = pwd;

%load the data
folderDir = dir(currentFolder);
analysisFiles = [];
for q = 1:size(folderDir)
    if contains(folderDir(q).name, '_data')
        subDirPath = [currentFolder filesep '*_data.mat'];
        subFolder = dir(subDirPath);
        analysisFiles = [analysisFiles; subFolder];
        break;
    end
end

numFiles = size(analysisFiles, 1);
baselineWindow = [0 900]; %first 15 mins
infusionWindow = [5400 6300]; %15 mins 1 hour after the infusion

for n = 1:numFiles
    load([analysisFiles(n).folder filesep analysisFiles(n).name]);
    
    %extract spikeTimesSec and ISIs from the dat structure
    allCells(n).baselineNumSpikes = length(dat.spikeTimesSec(dat.spikeTimesSec<901));
    allCells(n).spikeTimesSec = dat.spikeTimesSec;
    allCells(n).ISIs = dat.ISIs;

    %initialize variables
    %baseline variables
    allCells(n).bFR = [];
    allCells(n).bNumBursts = [];
    allCells(n).bSpikesPerBurst = [];
    allCells(n).bSFB = [];
    allCells(n).bNumPauses = [];
    allCells(n).bLenPause = [];
    allCells(n).bPSI = [];
    allCells(n).btonicISI = [];
    %post-DART capture variables
    allCells(n).iFR = [];
    allCells(n).iNumBursts = [];
    allCells(n).iSpikesPerBurst = [];
    allCells(n).iSFB = [];
    allCells(n).iNumPauses = [];
    allCells(n).iLenPause = [];
    allCells(n).iPSI = [];
    allCells(n).itonicISI = [];
end

%go through each cell
for n = 1:numFiles
    %calculate FRs
    bcurrSpikeTimesSec = allCells(n).spikeTimesSec(allCells(n).spikeTimesSec > baselineWindow(1) & allCells(n).spikeTimesSec < baselineWindow(2));
    allCells(n).bFR = length(bcurrSpikeTimesSec)/900;
    icurrSpikeTimesSec = allCells(n).spikeTimesSec(allCells(n).spikeTimesSec > infusionWindow(1) & allCells(n).spikeTimesSec < infusionWindow(2));
    allCells(n).iFR = length(icurrSpikeTimesSec)/900;

    %-----BURST CALCULATIONS-------
    %Burst: 3-10 spikes of decreasing amplitude, started by an ISI of <80ms and ended by a >160ms ISI
    Burst_ISI_Thresh_Sec = [.08 0.16];
   
    %--Calculation for the baseline window
    baselineISI = diff(bcurrSpikeTimesSec); %Extract ISIs for this window
    baselineISI = rmmissing(baselineISI); %remove NaNs
    baselineBurstLen = []; %initialize
    TMP = 0;  %temporary variable to count the number of ISI in this burst
    for b=1:length(baselineISI)
        %if <80ms and burst has not started, begin the burst
        if TMP == 0
            if baselineISI(b) <= Burst_ISI_Thresh_Sec(1)
                TMP = TMP + 1;
            end
            %once started calculating burst, continue adding to the TMP
            %variable for each ISI < 160ms
        else
            if baselineISI(b)<= Burst_ISI_Thresh_Sec(2)
                TMP = TMP + 1;
            else %as soon as we have an ISI > 160ms, save burst length and reset
                baselineBurstLen(end+1) = TMP+1; %add 1 - because this is the count of ISIs, we need to consider an additional "ending" spike as within the burst
                TMP = 0;
            end
        end
    end
    baselineBurstLen = baselineBurstLen(baselineBurstLen>2); %must have 3 spikes or more to be considered a burst
    allCells(n).bNumBursts = length(baselineBurstLen);    %number of bursts
    allCells(n).bSpikesPerBurst = mean(baselineBurstLen);    %median spikes per burst
    allCells(n).bSFB = sum(baselineBurstLen) / (length(baselineISI)+1); %percentage of spikes fired in bursts out of all spikes

    %--Calculation for the post-DART window
    %15 min bin near the end of the recording
    infuseISI = diff(icurrSpikeTimesSec); %Extract ISIs for this window
    infuseISI = rmmissing(infuseISI); %remove nan
    infuseBurstLen = []; %initialize
    TMP = 0;  %temporary variable to count the number of ISI in this burst
    for b=1:length(infuseISI)
        if TMP == 0
        %if <80ms and burst has not started, begin the burst
            if infuseISI(b) <= Burst_ISI_Thresh_Sec(1) && TMP == 0
                TMP = TMP + 1;
            end
        else
            %we are already in a burst, so criteria is easier
            if infuseISI(b) <= Burst_ISI_Thresh_Sec(2)
                TMP = TMP + 1;
            else %as soon as we have an ISI > 160ms, save burst length and reset
                infuseBurstLen(end+1) = TMP+1; %add 1 - because this is the count of ISIs, we need to consider an additional "ending" spike as within the burst
                TMP = 0;
            end
        end
    end
    infuseBurstLen = infuseBurstLen(infuseBurstLen>2); %must have 3 spikes or more to be considered a burst
    allCells(n).iNumBursts = length(infuseBurstLen);    %number of bursts
    allCells(n).iSpikesPerBurst = mean(infuseBurstLen);    %median spikes per burst
    allCells(n).iSFB = sum(infuseBurstLen) / (length(infuseISI)+1); %percentage of spikes fired in bursts out of all spikes

    %----------Pausing analysis------------%
    %reduction in activity beyond average for this window. ISI > tonic ISI
    %calculate the baseline tonic ISI rate
    tonicRate = median(baselineISI);
    allCells(n).btonicISI = tonicRate;
    baseline_pause_thresh = tonicRate*2;
    %calculate the post-DART tonic ISI rate
    infusionRate = median(infuseISI);
    allCells(n).itonicISI = infusionRate;
    infusion_pause_thresh = infusionRate*2;

    baselinePauses = baselineISI(baselineISI > baseline_pause_thresh);
    transientPauses = infuseISI(infuseISI > infusion_pause_thresh);

    %number of pauses
    allCells(n).bNumPauses = length(baselinePauses);
    allCells(n).iNumPauses = length(transientPauses);
    %mean length of pauses
    allCells(n).bLenPause = mean(baselinePauses)/tonicRate;
    allCells(n).iLenPause = mean(transientPauses)/infusionRate;
    % the percent of ISI intervals that are a pause
    allCells(n).bPSI = length(baselinePauses)/length(baselineISI);
    allCells(n).iPSI = length(transientPauses)/length(infuseISI);
end

%save the grouped spiking analysis
filename = [currentFolder filesep condition '_grouped spiking analysis_n' int2str(numFiles) '.mat'];
save(filename, 'allCells');

end
