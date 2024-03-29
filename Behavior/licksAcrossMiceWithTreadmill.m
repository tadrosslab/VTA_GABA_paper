%Written 12/12/2018 by SB, with edits 2019, 2020, and 2021
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%function to load the schultz data structure per mouse and combine data
%across a cohort

%Instructions:
%   Copy all analysis_tSchultz data structures from one cohort into a folder
%   Go to that folder
%   Input: 'cohort', a string for later saving the calculations performed
%   Outputs: 
%       p: a data structure with across-mouse calculations, saved
%       plotting of various features across mice

function licksAcrossMiceWithTreadmill(cohort)

currentFolder = pwd;

%load all the analysis files
folderDir = dir(currentFolder);
analysisFiles = [];
for q = 1:size(folderDir)
    if contains(folderDir(q).name, 'Mouse')
        subDirPath = [currentFolder filesep '*_analysis_tSchultz*.mat'];
        subFolder = dir(subDirPath);
        analysisFiles = [analysisFiles; subFolder];
        break;
    end
end

numFiles = size(analysisFiles, 1);
filtNum = 7; %for filtfilt in later plotting
probes = [];

%initialize a large structure for holding data from all mice
for n = 1:numFiles
    for k = 1:2
        allMice(n).schultz(k).performance = [];
        allMice(n).schultz(k).sessionIndx = [];
        allMice(n).schultz(k).flip = [];
        allMice(n).schultz(k).motivation = [];
        allMice(n).schultz(k).dayIndex = []; 
        allMice(n).schultz(k).firstLick = [];
        allMice(n).schultz(k).anticipatorySpeed = [];
        allMice(n).schultz(k).collectionSpeed = [];
        allMice(n).schultz(k).anticipatorySpeedNorm = [];
        allMice(n).schultz(k).collectionSpeedNorm = [];
        allMice(n).normPerformance(k).performance = [];
        allMice(n).normPerformance(k).norm = [];
        allMice(n).normPerformance(k).FL = [];
        allMice(n).normPerformance(k).motiv = [];
        allMice(n).normPerformance(k).startperformance = [];
        allMice(n).normPerformance(k).aSpeed = [];
        allMice(n).normPerformance(k).cSpeed = [];
        allMice(n).normPerformance(k).aSpeedNorm = [];
        allMice(n).normPerformance(k).cSpeedNorm = [];
    end
end

%initialize vars for calculating data across mice
for k = 1:2
    p(k).performance = [];
    p(k).norm = [];
    p(k).firstlick = [];
    p(k).m = [];
    p(k).SEM = [];
    p(k).flipratio = [];
    p(k).startperformance = [];
    p(k).probes = [];
    p(k).aSpeed = [];
    p(k).cSpeed = [];
    p(k).aSpeedNorm = [];
    p(k).cSpeedNorm = [];
end

%flip occurs ~1/4 of the way through the flip day, so plot as a 25/75 split
flipratio(1) = 25;
flipratio(2) = 100-flipratio(1);

%get files and perform analysis/resizing
for i = 1:numFiles 
    load([analysisFiles(i).folder filesep analysisFiles(i).name]);

    for k = 1:2     
        allMice(i).schultz(k).performance = schultz(k).performance(1, :);
        allMice(i).schultz(k).sessionIndx = unique(schultz(k).sessionIndx); %double check removed doubles
        allMice(i).schultz(k).flip = schultz(k).flip(1);
        allMice(i).schultz(k).motivation = schultz(k).motivation;
        allMice(i).schultz(k).dayIndex = schultz(k).dayIndex; 
        allMice(i).schultz(k).firstLick = schultz(k).firstLick;
        allMice(i).schultz(k).probeIndx = unique(schultz(k).probeIndx);
        allMice(i).schultz(k).anticipatorySpeed = schultz(k).anticipatorySpeed;
        allMice(i).schultz(k).collectionSpeed = schultz(k).collectionSpeed;
    end
    
    %get mean of final training day for normalization
    normind = [];
    %only get mean of tone A anticipatory licking, not the background
    %licking
    for k = 2
        %get trial indices of training 10 tone A
        for n = 1:length(allMice(i).schultz(k).sessionIndx)
            if allMice(i).schultz(k).sessionIndx(n) > allMice(i).schultz(k).flip(1)
                normind = [(allMice(i).schultz(k).sessionIndx(n-2)+1) allMice(i).schultz(k).sessionIndx(n-1)]; %norm to training 10 only
                break
            end
        end
        %calculate baseline means for normalization
        m = mean(allMice(i).schultz(k).performance(1, normind(1):normind(2)));
        maSpeed = mean(allMice(i).schultz(k).anticipatorySpeed(1, normind(1):normind(2)));
        mcSpeed = mean(allMice(i).schultz(k).collectionSpeed(1, normind(1):normind(2)));
        fprintf('mean is: %0.4f for mouse %s\n', m, analysisFiles(i).name);
    end

    %probes are only in tone B
    %get the mean probe response of each mouse
    for k = 1
        probes(i) = mean(allMice(i).schultz(1).performance(allMice(i).schultz(1).probeIndx));
    end

    %normalize and resize trials (different numbers across mice) into 100
    %"trials" per behavior session
    for k = 1:2
        %normalize to the calculated means
        allMice(i).schultz(k).norm = allMice(i).schultz(k).performance(1,:)./m;
        allMice(i).schultz(k).anticipatorySpeedNorm = allMice(i).schultz(k).anticipatorySpeed(1, :)./maSpeed;
        allMice(i).schultz(k).collectionSpeedNorm = allMice(i).schultz(k).collectionSpeed(1, :)./mcSpeed;
        
        %--Final training and testing--
        %get indices of start of flip day, start of flip, and 2 days before and after
        %give 100 "slots" per day and .25/.75 each for pre and post flip
        for n = 1:length(allMice(i).schultz(k).sessionIndx)
            if allMice(i).schultz(k).sessionIndx(n) > allMice(i).schultz(k).flip(1)
                days = [allMice(i).schultz(k).sessionIndx(n-3)
                    allMice(i).schultz(k).sessionIndx(n-2)
                   allMice(i).schultz(k).sessionIndx(n-1)
                   schultz(k).flip(1)
                   allMice(i).schultz(k).sessionIndx(n)
                   allMice(i).schultz(k).sessionIndx(n+1)];
               days(:, 2) = [100 100 flipratio(1) flipratio(2) 100 100];
               break
            end
        end

        %resize all variables to our preallocated slots, creating one large
        %row vector
        for td = 1:length(days)-1
            %performance
            allMice(i).normPerformance(k).performance = [allMice(i).normPerformance(k).performance imresize(allMice(i).schultz(k).performance(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];
           
            %normalized performance
            allMice(i).normPerformance(k).norm = [allMice(i).normPerformance(k).norm imresize(allMice(i).schultz(k).norm(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];
           
            %time to first lick
            allMice(i).normPerformance(k).FL = [allMice(i).normPerformance(k).FL imresize(allMice(i).schultz(k).firstLick(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];
                        
            %Anticipatory speed
            allMice(i).normPerformance(k).aSpeed = [allMice(i).normPerformance(k).aSpeed imresize(allMice(i).schultz(k).anticipatorySpeed(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];  
            
            %Collection speed
            allMice(i).normPerformance(k).cSpeed = [allMice(i).normPerformance(k).cSpeed imresize(allMice(i).schultz(k).collectionSpeed(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];  
            
            %normalized Anticipatory speed
            allMice(i).normPerformance(k).aSpeedNorm = [allMice(i).normPerformance(k).aSpeedNorm imresize(allMice(i).schultz(k).anticipatorySpeedNorm(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];
            
            %normalized Collection speed
            allMice(i).normPerformance(k).cSpeedNorm = [allMice(i).normPerformance(k).cSpeedNorm imresize(allMice(i).schultz(k).collectionSpeedNorm(1, days(td, 1):days(td+1, 1)), [1 days(td, 2)], 'nearest')];  
        end
        
        %--Initial Learning--
        %get indices of all training days to look at anticipatory licking
        %across training
        %give 100 "slots" per day
        startdays = [allMice(i).schultz(k).sessionIndx(1)
            allMice(i).schultz(k).sessionIndx(2)
            allMice(i).schultz(k).sessionIndx(3)
            allMice(i).schultz(k).sessionIndx(4)
            allMice(i).schultz(k).sessionIndx(5)
            allMice(i).schultz(k).sessionIndx(6)
            allMice(i).schultz(k).sessionIndx(7)
            allMice(i).schultz(k).sessionIndx(8)
            allMice(i).schultz(k).sessionIndx(9)
            allMice(i).schultz(k).sessionIndx(10)
            allMice(i).schultz(k).sessionIndx(11)];
        startdays(:, 2) = [100 100 100 100 100 100 100 100 100 100 100]; %variable to look at beginning of performance
        startdays(1, 1) = 1;
        
        %resize all variables to our preallocated slots
        for sd = 1:length(startdays)-1
            %give 100 "slots" per day 
            %starting performance, no normalization wanted
            allMice(i).normPerformance(k).startperformance = [allMice(i).normPerformance(k).startperformance imresize(allMice(i).schultz(k).performance(1, startdays(sd, 1):startdays(sd+1, 1)), [1 startdays(sd, 2)], 'nearest')];
        end
        
        %add data per mouse as a new row into our p variable
        p(k).performance = [p(k).performance; allMice(i).normPerformance(k).performance];
        p(k).norm = [p(k).norm; allMice(i).normPerformance(k).norm];
        p(k).firstlick = [p(k).firstlick; allMice(i).normPerformance(k).FL];
        p(k).startperformance = [p(k).startperformance; allMice(i).normPerformance(k).startperformance];        
        p(k).aSpeed = [p(k).aSpeed; allMice(i).normPerformance(k).aSpeed];
        p(k).cSpeed = [p(k).cSpeed; allMice(i).normPerformance(k).cSpeed];
        p(k).aSpeedNorm = [p(k).aSpeedNorm; allMice(i).normPerformance(k).aSpeedNorm];
        p(k).cSpeedNorm = [p(k).cSpeedNorm; allMice(i).normPerformance(k).cSpeedNorm];
        
        %get the probe values from our normalized schultz.norm anticipatory licking variable
        try
            p(k).probes = [p(k).probes; allMice(i).schultz(k).norm(allMice(i).schultz(k).probeIndx)];
        catch
            fprintf('no probe data!');
        end
    end
 
end

%save p variable
filename = [currentFolder filesep cohort '_grouped analysis_n' int2str(numFiles) '.mat'];
save(filename, 'p');

%--------------PLOT DATA---------------------
%initialize various variables for the graphs
x = -2:.01:1.99;
startidx = 20;
flipidx = 200 + flipratio(1) + 1;
grey = [.5 .5 .5];
lines = zeros(1, 3);
filt = designfilt('lowpassfir', 'FilterOrder', filtNum, 'CutoffFrequency', .01);

%---------FIGURE 1: Plot non-normalized anticipatory licking------------%
for k = 1:2
    p(k).m = mean(p(k).performance); %mean
    SEM = std(p(k).performance)/sqrt(numFiles); %SEM
    p(k).m = [p(k).m; p(k).m+SEM; p(k).m-SEM]; %combine into one variable
end

figure(1); clf; hold on;
for k = 1:2   
    for i = 1:3
        meansmooth(i, :) = filter(filt, p(k).m(i, :)); %filter the mean and SEM bars
    end

    if k == 1
        %plot the background licking (pre-flip)
        ph1 = patch([x(startidx:flipidx) fliplr(x(startidx:flipidx))], [meansmooth(2, startidx:flipidx) fliplr(meansmooth(3, startidx:flipidx))], grey, 'edgecolor', 'none');
        set(ph1, 'facealpha', .3);
        %plot the conditioning ant licking
        h1 = patch([x(flipidx:end) fliplr(x(flipidx:end))], [meansmooth(2, flipidx:end) fliplr(meansmooth(3, flipidx:end))], 'r', 'edgecolor', 'none');
        set(h1, 'facealpha', .3);
        lines(1, 1)=plot(x(startidx:flipidx), meansmooth(1, startidx:flipidx), 'color', grey, 'linewidth', 2);
        lines(1, 2)=plot(x(flipidx:end), meansmooth(1, flipidx:end), 'color', 'r', 'linewidth', 2);
    else
        %plot the extinction ant licking
        h2 = patch([x(startidx:end) fliplr(x(startidx:end))], [meansmooth(2, startidx:end) fliplr(meansmooth(3, startidx:end))], 'c', 'edgecolor', 'none');
        set(h2, 'facealpha', .3);
        lines(1, 3)=plot(x(startidx:end), meansmooth(1, startidx:end), 'color', 'c', 'linewidth', 2);
    end
end
%label graph
ylimit = get(gca, 'ylim');
line([flipratio(1)/100 flipratio(1)/100], ylimit, 'color', 'k');
set(gca,'fontsize',20);
title([cohort ', n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Anticipatory licking');
xticks(-2:1:3);
xlim([-1 2])
ylim([0 1]);
xlimit = get(gca, 'xlim');
plot(xlimit, [.25 .25], ':', 'color', 'k')
plot(xlimit, [.5 .5], ':', 'color', 'k');
plot(xlimit, [.75 .75], ':', 'color', 'k');
legend([lines(1, 3) lines(1, 2) lines(1, 1)], {'Tone A', 'Tone B', 'Background'}, 'Location', 'northwest');

%---------FIGURE 2: plot time to first lick------------%
for k = 1:2
    p(k).firstlick(p(k).firstlick(:, :) > 5.5) = 5.5; %collapse all long lick delays into one value
    p(k).m = mean(p(k).firstlick); %mean
    SEM = std(p(k).firstlick)/sqrt(numFiles); %SEM
    p(k).m = [p(k).m; p(k).m+SEM; p(k).m-SEM]; %combine into one variable
end

figure(2); clf; hold on;
for k = 1:2
    %generate filtered mean and SEM data, with a filter break between each
    %behavioral session
    for i = 1:3
        lmeansmooth(i, 1:100) = filtfilt(filt, p(k).m(i, 1:100));
        lmeansmooth(i, 101:200) = filtfilt(filt, p(k).m(i, 101:200));
        lmeansmooth(i, 201:300) = filtfilt(filt, p(k).m(i, 201:300));
        lmeansmooth(i, 301:400) = filtfilt(filt, p(k).m(i, 301:400));
    end

    %plot the conditioning time to first lick lines
    %no point plotting the baseline for this graph
    if k == 1
        lh1 = patch([x(flipidx:300) fliplr(x(flipidx:300))], [lmeansmooth(2, flipidx:300) fliplr(lmeansmooth(3, flipidx:300))], 'r', 'edgecolor', 'none');
        set(lh1, 'facealpha', .3);
        lh1 = patch([x(301:end) fliplr(x(301:end))], [lmeansmooth(2, 301:end) fliplr(lmeansmooth(3, 301:end))], 'r', 'edgecolor', 'none');
        set(lh1, 'facealpha', .3);
        plot(x(flipidx:300), lmeansmooth(1, flipidx:300), 'color', 'r', 'linewidth', 2);      
        plot(x(301:end), lmeansmooth(1, 301:end), 'color', 'r', 'linewidth', 2);
    %plot all the extinction time to first lick lines
    else
        lh2 = patch([x(startidx:100) fliplr(x(startidx:100))], [lmeansmooth(2, startidx:100) fliplr(lmeansmooth(3, startidx:100))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(startidx:100), lmeansmooth(1, startidx:100), 'color', 'c', 'linewidth', 2);
        
        lh2 = patch([x(101:200) fliplr(x(101:200))], [lmeansmooth(2, 101:200) fliplr(lmeansmooth(3, 101:200))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(101:200), lmeansmooth(1, 101:200), 'color', 'c', 'linewidth', 2);
        
        lh2 = patch([x(201:300) fliplr(x(201:300))], [lmeansmooth(2, 201:300) fliplr(lmeansmooth(3, 201:300))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(201:300), lmeansmooth(1, 201:300), 'color', 'c', 'linewidth', 2);
        
        lh2 = patch([x(301:end) fliplr(x(301:end))], [lmeansmooth(2, 301:end) fliplr(lmeansmooth(3, 301:end))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(301:end), lmeansmooth(1, 301:end), 'color', 'c', 'linewidth', 2);
    end
end

%label figure
line([flipratio(1)/100 flipratio(1)/100], [-2 5], 'color', 'k');
set(gca,'fontsize',50);
title([cohort ', n=' int2str(numFiles)]);
ylabel('Time to first lick (s)');
xlabel('Session');
xticks(-2:1:3);
xlim([-1 2]);
ylim([0 5]);
set(gca, 'YDir', 'reverse');

%---------FIGURE 3: Plot normalized anticipatory licking with filter breaks------------%
for k = 1:2
    p(k).m = mean(p(k).norm); %mean
    SEM = std(p(k).norm)/sqrt(numFiles); %SEM
    p(k).m = [p(k).m; p(k).m+SEM; p(k).m-SEM]; %combine variables
end

figure(3); clf; hold on;
for k = 1:2
    %generate filtered mean and SEM data, with a filter break between each
    %behavioral session
    for i = 1:3
        anmeansmooth(i, 1:100) = filtfilt(filt, p(k).m(i, 1:100));
        anmeansmooth(i, 101:200) = filtfilt(filt, p(k).m(i, 101:200));   
        anmeansmooth(i, 201:300) = filtfilt(filt, p(k).m(i, 201:300));
        anmeansmooth(i, 301:400) = filtfilt(filt, p(k).m(i, 301:400));                  
    end

    if k == 1
        %plot the background licking (pre-flip)
        anph1 = patch([x(startidx:100) fliplr(x(startidx:100))], [anmeansmooth(2, startidx:100) fliplr(anmeansmooth(3, startidx:100))], grey, 'edgecolor', 'none');
        set(anph1, 'facealpha', .3);
        anph1 = patch([x(101:200) fliplr(x(101:200))], [anmeansmooth(2, 101:200) fliplr(anmeansmooth(3, 101:200))], grey, 'edgecolor', 'none');
        set(anph1, 'facealpha', .3);
        anph1 = patch([x(201:flipidx) fliplr(x(201:flipidx))], [anmeansmooth(2, 201:flipidx) fliplr(anmeansmooth(3, 201:flipidx))], grey, 'edgecolor', 'none');
        set(anph1, 'facealpha', .3);
        %plot the conditioning anticipatory licking (post-flip)
        anh1 = patch([x(flipidx:300) fliplr(x(flipidx:300))], [anmeansmooth(2, flipidx:300) fliplr(anmeansmooth(3, flipidx:300))], 'r', 'edgecolor', 'none');
        set(anh1, 'facealpha', .3);
        anh1 = patch([x(301:end) fliplr(x(301:end))], [anmeansmooth(2, 301:end) fliplr(anmeansmooth(3, 301:end))], 'r', 'edgecolor', 'none');
        set(anh1, 'facealpha', .3);
        lines(1, 1)=plot(x(startidx:100), anmeansmooth(1, startidx:100), 'color', grey, 'linewidth', 3);
        lines(1, 1)=plot(x(101:200), anmeansmooth(1, 101:200), 'color', grey, 'linewidth', 3);
        lines(1, 1)=plot(x(201:flipidx), anmeansmooth(1, 201:flipidx), 'color', grey, 'linewidth', 3);
        lines(1, 2)=plot(x(flipidx:300), anmeansmooth(1, flipidx:300), 'color', 'r', 'linewidth', 3);
        lines(1, 2)=plot(x(301:end), anmeansmooth(1, 301:end), 'color', 'r', 'linewidth', 3);
    else
        %plot the training ant licking
        lh2 = patch([x(startidx:100) fliplr(x(startidx:100))], [anmeansmooth(2, startidx:100) fliplr(anmeansmooth(3, startidx:100))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(startidx:100), anmeansmooth(1, startidx:100), 'color', 'c', 'linewidth', 2);
        
        lh2 = patch([x(101:200) fliplr(x(101:200))], [anmeansmooth(2, 101:200) fliplr(anmeansmooth(3, 101:200))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(101:200), anmeansmooth(1, 101:200), 'color', 'c', 'linewidth', 2);
        
        %plot the extinction ant licking
        lh2 = patch([x(201:300) fliplr(x(201:300))], [anmeansmooth(2, 201:300) fliplr(anmeansmooth(3, 201:300))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(201:300), anmeansmooth(1, 201:300), 'color', 'c', 'linewidth', 2);
        
        lh2 = patch([x(301:end) fliplr(x(301:end))], [anmeansmooth(2, 301:end) fliplr(anmeansmooth(3, 301:end))], 'c', 'edgecolor', 'none');
        set(lh2, 'facealpha', .3);
        plot(x(301:end), anmeansmooth(1, 301:end), 'color', 'c', 'linewidth', 2);
    end
end

%label graph
line([flipratio(1)/100 flipratio(1)/100], [0 2], 'color', 'k');
set(gca,'fontsize',20);
title([cohort ' with norm, n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Anticipatory licking (normalized)');
xticks(-2:1:3);
xlim([-1 2]);
xlimit = get(gca, 'xlim');
plot(xlimit, [.5 .5], ':', 'color', 'k');
plot(xlimit, [1 1], ':', 'color', 'k');
legend([lines(1, 3) lines(1, 2) lines(1, 1)], {'Tone A', 'Tone B', 'Background'}, 'Location', 'northwest');
ylim([0 1.5]);


%---------FIGURE 4: Plot learning rates during training------------%
for k = 1:2
    p(k).m = mean(p(k).startperformance); %mean
    SEM = std(p(k).startperformance)/sqrt(numFiles); %SEM
    p(k).m = [p(k).m; p(k).m+SEM; p(k).m-SEM]; %combine variables
end
x = 0:.01:9.99;
xbins = [0 100 200 300 400 500 600 700 800 900 1000];
figure(4); clf; hold on;
for k = 1:2
    %generate filtered mean and SEM data, with a filter break between each
    %behavioral session
    for i = 1:3
        for n=1:10
            startmeansmooth(i, xbins(n)+1:xbins(n+1)) = filtfilt(filt, p(k).m(i, xbins(n)+1:xbins(n+1)));
        end
    end
    if k == 1
        %plot the background licking (pre-flip)
        for n = 1:10
            sph1 = patch([x(xbins(n)+1:xbins(n+1)) fliplr(x(xbins(n)+1:xbins(n+1)))], [startmeansmooth(2, xbins(n)+1:xbins(n+1)) fliplr(startmeansmooth(3, xbins(n)+1:xbins(n+1)))], [.6 .6 .6], 'edgecolor', 'none');
            set(sph1, 'facealpha', .3);
            lines(1, 1)=plot(x(xbins(n)+1:xbins(n+1)), startmeansmooth(1, xbins(n)+1:xbins(n+1)), 'color', [.6 .6 .6], 'linewidth', 3);
        end
    else
        %plot the training tone A licking
        for n=1:10
        sh2 = patch([x(xbins(n)+1:xbins(n+1)) fliplr(x(xbins(n)+1:xbins(n+1)))], [startmeansmooth(2, xbins(n)+1:xbins(n+1)) fliplr(startmeansmooth(3, xbins(n)+1:xbins(n+1)))], 'c', 'edgecolor', 'none');
        set(sh2, 'facealpha', .3);
        lines(1, 3)=plot(x(xbins(n)+1:xbins(n+1)), startmeansmooth(1, xbins(n)+1:xbins(n+1)), 'color', 'c', 'linewidth', 3);
        end
    end
end

%label plot
set(gca,'fontsize',20);
title([cohort ', n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Anticipatory licking');
xticks(0:1:10);
xlim([0 10]);
xlimit = get(gca, 'xlim');
legend([lines(1, 3) lines(1, 1)], {'Tone A', 'Background'}, 'Location', 'northwest');
ylim([0 1]);
grid on;

%-----------------FIGURE 5: Plot running speed during anticipatory window--------------------------
x = -2:.01:1.99;
for k = 1:2
    p(k).m = mean(p(k).aSpeed); %mean
    SEM = std(p(k).aSpeed)/sqrt(numFiles); %SEM
    p(k).m = [p(k).m; p(k).m+SEM; p(k).m-SEM]; %combine variables
end

figure(5); clf; hold on;
for k = 1:2   
    %generate filtered mean and SEM data
    for i = 1:3
        meansmooth(i, :) = filter(filt, p(k).m(i, :));
    end

    if k == 1
        %plot pre-flip background running
        ph1 = patch([x(startidx:flipidx) fliplr(x(startidx:flipidx))], [meansmooth(2, startidx:flipidx) fliplr(meansmooth(3, startidx:flipidx))], grey, 'edgecolor', 'none');
        %plot post-flip conditioning running
        h1 = patch([x(flipidx:end) fliplr(x(flipidx:end))], [meansmooth(2, flipidx:end) fliplr(meansmooth(3, flipidx:end))], 'r', 'edgecolor', 'none');
        set(h1, 'facealpha', .3);
        set(ph1, 'facealpha', .3);
        lines(1, 1)=plot(x(startidx:flipidx), meansmooth(1, startidx:flipidx), 'color', grey, 'linewidth', 2);
        lines(1, 2)=plot(x(flipidx:end), meansmooth(1, flipidx:end), 'color', 'r', 'linewidth', 2);
    else
        %plot training/extinction running
        h2 = patch([x(startidx:end) fliplr(x(startidx:end))], [meansmooth(2, startidx:end) fliplr(meansmooth(3, startidx:end))], 'c', 'edgecolor', 'none');
        set(h2, 'facealpha', .3);
        lines(1, 3)=plot(x(startidx:end), meansmooth(1, startidx:end), 'color', 'c', 'linewidth', 2);
    end
    
end
%label plot
ylimit = get(gca, 'ylim');
line([flipratio(1)/100 flipratio(1)/100], ylimit, 'color', 'k');
set(gca,'fontsize',20);
title([cohort ' Running, n=' int2str(numFiles)]);
xlabel('Session');
ylabel('Treadmill rotations per trial');
xticks(-2:1:3);
xlim([-2 2])
legend([lines(1, 3) lines(1, 2) lines(1, 1)], {'Tone A', 'Tone B', 'Background'}, 'Location', 'northwest');
end
