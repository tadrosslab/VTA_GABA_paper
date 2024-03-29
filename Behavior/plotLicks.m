%Written 11/2023 by SB with input from MRT
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%Function to plot licks as ticks across rewarded trials

%Instructions:
%   Inputes:
%       load trial data of interest into the workspace, then run this code
%       trialRange: a two number array of the trials you are interested in
%           (ex, [1 10])
%   output: plot of licks

function plotLicks(trialData, trialRange)

lickTimesIdx = trialData.lickValues(2,:)>0;
%get the actual time value matching the indices
lickTimes = trialData.recordTimes(1, lickTimesIdx);

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
%translate trial indices into time stamps
actualSucrose(2, :) = trialData.recordTimes(1, actualSucrose(1, :));

%array of time stamps we are interested (from -3 to 3 around reward delivery with a step size dt)
dt=0.01;
t = -1.5: dt : 1.5;

Ltmp = zeros(length(actualSucrose), length(t));
for k=1:size(actualSucrose, 2)
    SDT = actualSucrose(2, k);  %this is the time, in sec, of the reward delivery around which we are analyzing licks
    RelLickTime = lickTimes - SDT; %time of licks relative to reward delivery
    
    I =  (RelLickTime>=t(1)) &  (RelLickTime<=t(end));
    for j=1:length(t)
        dtToNearestLick = min( abs(t(j) - RelLickTime(I)));
        if (dtToNearestLick <= dt/2)
            Ltmp(k, j) = Ltmp(k, j) + 1;  %sum for just this one
        end
    end
end

%make each trial a new row
Ltmp(Ltmp==0) = NaN;
for i = 1:size(Ltmp,1)
    Ltmp(i, :) = Ltmp(i, :) * i;
end
figure(1); clf;
plot(t, Ltmp(trialRange(1):trialRange(2), :), 'color', 'k', 'LineWidth', 10);
xlim([-3 3]);
ylim([trialRange]);
ax = gca;
ax.YDir = 'reverse';
end