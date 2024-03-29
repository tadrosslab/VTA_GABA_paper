%Written 4/2022 by SB
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%function to get change in running at all times of a behavior session (not just when trials are happening)

%Instructions:
%   Copy all analysis_running data structures from one cohort into a folder
%   Go to that folder
%   Input: 'cohort', a string for later saving the calculations performed
%   Outputs: 
%       s: a data structure with across-mouse calculations, saved
%       plotting of RPM across mice

function runningAcrossMice(cohort)

currentFolder = pwd;

%load the trial data
folderDir = dir(currentFolder);
analysisFiles = [];
for q = 1:size(folderDir)
    if contains(folderDir(q).name, 'Mouse')
        subDirPath = [currentFolder filesep '*_analysis_running.mat'];
        subFolder = dir(subDirPath);
        analysisFiles = [analysisFiles; subFolder];
        break;
    end
end

numFiles = size(analysisFiles, 1);
filtNum = 5; %for filtfilt in later plotting

%initialize a large structure for holding data from all mice
for n = 1:numFiles
    allMice(n).totalRunning = [];
    allMice(n).dayIndex = [];
    allMice(n).RPM = [];
    allMice(n).RPMnorm = [];
    allMice(n).minIndex = [];
    allMice(n).resizeRPM = [];
    allMice(n).resizeRPMnorm = [];
end

%initialize vars for calculating data across mice
s.totRPM = [];
s.totRPMnorm = [];

%flip occurs ~1/4 of the way through the flip day, so plot as a 25/75 split
flipratio(1) = 25;
flipratio(2) = 100-flipratio(1);

%get files and copy data from each mouse into one structure
for i = 1:numFiles 
    load([analysisFiles(i).folder filesep analysisFiles(i).name]);
    
    allMice(i).totalRunning = totalSpeed(2, :);

    %find day indices based on recording times being 0 (i.e., the start of a recording)
    tmpDayIndex = find(totalSpeed(1, :) == 0);
    %first two indices are 0 per session, so don't double count - only extract the first index per session
    allMice(i).dayIndex = [1 tmpDayIndex(2) tmpDayIndex(4) tmpDayIndex(6) tmpDayIndex(8) tmpDayIndex(end)];

    %get rotations per minute (roughly 60000 NI card "checks")
    for n=2:length(allMice(i).dayIndex)
        NICount = allMice(i).dayIndex(n-1);
        while NICount < (allMice(i).dayIndex(n)-60000)
            currCPR = sum(diff(totalSpeed(2, NICount:(NICount+59999))));
            currRPM = currCPR/100; %Our rotary encoder has 100 cycles per revolution, so divide by 100 to get the number of revolutions
            allMice(i).RPM = [allMice(i).RPM currRPM];
            NICount = NICount + 60000;
        end
        %record the last minute of this session
        allMice(i).minIndex = [allMice(i).minIndex length(allMice(i).RPM)];
    end
    
    %get mean of pre-DART days for normalization
    normind = [1 allMice(i).minIndex(1) allMice(i).minIndex(2) allMice(i).minIndex(3) allMice(i).minIndex(4) allMice(i).minIndex(5)];    
    meanRPM = mean(allMice(i).RPM(normind(1):normind(4)));
    allMice(i).RPMnorm = allMice(i).RPM./meanRPM;
   
    %just in case, to present bugs, resize to 60 mins if necessary    
    resizeIdx = [60 60 60 60 60];
    for n=1:5
        allMice(i).resizeRPM = [allMice(i).resizeRPM imresize(allMice(i).RPM(normind(n):normind(n+1)), [1 resizeIdx(n)], 'nearest')];  
        allMice(i).resizeRPMnorm = [allMice(i).resizeRPMnorm imresize(allMice(i).RPMnorm(normind(n):normind(n+1)), [1 resizeIdx(n)], 'nearest')];
    end
    
    %add raw and normalized data to the s variable as a new row per mouse
    s.totRPM = [s.totRPM; allMice(i).resizeRPM];
    s.totRPMnorm = [s.totRPMnorm; allMice(i).resizeRPMnorm];
end

%save the s variable with the resized and combined data
filename = [currentFolder filesep cohort '_grouped running analysis_n' int2str(numFiles) '.mat'];
save(filename, 's');

%---------PLOTTING---------------
td = 1/60; %60 minutes per session
x = -3:td:1.99;
grey = [.5 .5 .5];
lines = zeros(1, 3);
filt = designfilt('lowpassfir', 'FilterOrder', filtNum, 'CutoffFrequency', .01);

%---------FIGURE 1:RPM------------%
s.m = mean(s.totRPM); %mean
SEM = std(s.totRPM)/sqrt(numFiles); %SEM
s.m = [s.m; s.m+SEM; s.m-SEM]; %combine into one figure

figure(1); clf; hold on;
for i = 1:3
    meansmooth(i, :) = filter(filt, s.m(i, :)); %filter the data
end

%plot the data
h2 = patch([x fliplr(x)], [meansmooth(2, :) fliplr(meansmooth(3, :))], 'r', 'edgecolor', 'none');
set(h2, 'facealpha', .3);
plot(x, meansmooth(1, :), 'color', 'r', 'linewidth', 2);
%add details to the graph
ylimit = get(gca, 'ylim');
line([0 0], ylimit, 'color', 'k');
set(gca,'fontsize',20);
title([cohort ' running, n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Rotary encoder RPM');
xticks(-1:1:2);
xlim([-3 2]);
xlimit = get(gca, 'xlim');

%---------FIGURE 2:RPM normalized------------%
s.m = mean(s.totRPMnorm); %mean
SEM = std(s.totRPMnorm)/sqrt(numFiles); %SEM
s.m = [s.m; s.m+SEM; s.m-SEM]; %combine for simplicity into one variable

figure(2); clf; hold on;
for i = 1:3
    meansmooth(i, :) = filter(filt, s.m(i, :)); %filter the data
end

%plot the data
h2 = patch([x fliplr(x)], [meansmooth(2, :) fliplr(meansmooth(3, :))], 'r', 'edgecolor', 'none');
set(h2, 'facealpha', .3);
plot(x, meansmooth(1, :), 'color', 'r', 'linewidth', 2);
%add details to the graph
ylimit = get(gca, 'ylim');
line([0 0], ylimit, 'color', 'k');
set(gca,'fontsize',20);
title([cohort ' running, n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Rotary encoder RPM normalized');
xticks(-1:1:2);
xlim([-3 2])
xlimit = get(gca, 'xlim');

%---------FIGURE 3:RPM normalized with day breaks for filtering------------%
xlims = [0 60 120 180 240 300]; %the day boundaries
figure(3); clf; hold on;
filtfull = [];
%apply the filter for each behavior session independently
for k = 1:length(xlims)-1
    tmp = s.m(:, xlims(k)+1:xlims(k+1));
    firstval = s.m(:, xlims(k)+1); 
    for i = 1:3
        tmpsmooth(i, :) = filter(filt, tmp(i, :)); 
    end
    filtfull = [filtfull tmpsmooth];
end
%plot the data
h2 = patch([x fliplr(x)], [filtfull(2, :) fliplr(filtfull(3, :))], 'r', 'edgecolor', 'none');
set(h2, 'facealpha', .3);
plot(x, filtfull(1, :), 'color', 'r', 'linewidth', 2);
%add graph details
ylimit = get(gca, 'ylim');
line([0 0], ylimit, 'color', 'k');
set(gca,'fontsize',20);
title([cohort ' running, n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Rotary encoder RPM normalized');
xticks(-1:1:2);
xlim([-3 2])
xlimit = get(gca, 'xlim');

end