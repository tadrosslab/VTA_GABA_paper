%Written by SB

%Code to run a solenoid to (1) run liquid through tubing or (2) check
%volume delivered by the solenoid

function calibrateSolenoid_NI

sol = daq.createSession('ni');
sol.Rate = 5000;

addAnalogOutputChannel(sol, 'Dev1', 'ao1', 'Voltage');   %the solenoid channel

offDuration = .5; %wait 500ms between each click
offDurScans = offDuration*5000;

%consistently runs until the user enters a value of '0' for the open
%duration query
while 1 
    stim = [];
    openDuration = input('Set open duration in ms as: '); %user types input into command window
    numClicks = input('How many clicks: '); %user types input into command window
    if openDuration == 0
        break
    end
    openDurScans = openDuration*5;

    %prep the NI stimulation data
    stim = zeros((offDurScans+openDurScans), 1); %one long column
    stim(1:openDurScans) = 5; %openDurScans set to high
    stim(openDurScans + 1:end) = 0; %offDurScans set to low
    stimFull = repmat(stim, numClicks, 1); %repeat the stimulus numClick times

    queueOutputData(sol, stimFull);
    startForeground(sol); %run the solenoid
end

stop(sol);
end