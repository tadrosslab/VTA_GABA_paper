%Written by SB and MRT

%Function for running the conditional learning lick assay box with the treadmill setup and two tones
%Input: 
%   1) solOpenDuration: ms the solenoid should be open for

%Output:
%   1) trialData: structure containing the mouse information and all the time stamps
%       *output vars in order of saving*
%       a)trialData.cohort: cohort name
%       b)trialData.mouse: mouse number
%       c)trialData.cage: cage number
%       d)trialData.group: control, cannula saline, cannula DART, etc.
%       e)trialData.trialType: training (pre flip) and testing (post flip)
%       f)trialData.recordTimes: scalar vector of the record times; 
%       g)trialData.toneValues: Array of time stamps indicating every time a tone was presented with rows 1) frequency of tone. Data points should line up with trialData.recordTimes
%       h)trialData.sucroseDeliveryValues: Scalar vector of the voltage output to the solenoid. Data points should line up with trialData.recordTimes
%       i)trialData.lickValues: Array with rows: 1) IR voltage 2) binary yes/no (1/0) as to whether that measurement counted as a lick while the code was running. Data points should line up with trialData.recordTimes
%       j)trialData.videoTimes: Scalar vector of time stamps indicating every time video data was retrieved. Time stamps are negative if there was an error in video retrieval.
%       k)trialData.rotationalSpeed: Scalar vector with rows (1) Edgecount from rotary encoder. Data points should line up with trialData.recordTimes
%       l)trialData.trialInfo: overall reference of trials. Array with rows (1) trial number (2) tone frequency (3) 1/0 rewarded or not (4) interval 
%               (5) time started from matlab clock (6) time stim output
%               starts (aka when tone plays/would play) from matlab clock
%               (7) time solenoid closed from NI clock (8) 1/0 probe or not
%       m)trialData.flipTrial: the trial that the flip occurred on. 0 if not a flip day
%       n)trialData.lickThreshold: threshold in volts above which is considered a lick
%       o)trialData.baselineThreshold: IR threshold baseline
%       p)trialData.numTrials: total number of trials run
%       q)trialData.numRewardedTrials: total number of rewarded trials run
%       r)trialData.numUnrewardedTrials: total number of unrewarded trials run. Includes probe trials
%       s)trialData.numProbeTrials: total number of probe trials run
%       t)trialData.interval: delay interval
%       u)trialData.solenoidOpenDuration: solenoid opening in ms
%       v)trialData.rewardDelay: delay between tone and reward (or start of next trial)
%       w)trialData.toneDurationSec: length of the tone
%       x)trialData.rewardedSoundHz: rewarded tone frequency
%       y)trialData.unrewardedSoundHz: unrewarded tone frequency
%       z)trialData.dayType: 1 for training, 2 for flip, 3 for testing
%       aa)trialData.dateTime: string with the date and start time 
%       bb)trialData.box: 1 or 2 for the behavior box

function lickVid_dualConditionExtinction_NI(solOpenDuration)

%make trialData and videoData global variables so they are not lost even if the code crashes
clear global trialData
global trialData videoData NumVideoFrames NI_INPUT
videoData = zeros(300, 500, 100000, 'uint8'); %initialize memory
NumVideoFrames = 0;

%__________________________________________________________________________
%-----------------USER INPUT----------------------
%run time in num trials
numTrials = 600;

%run time in num seconds
numSeconds = 3600;

toneA = 2500;
toneB = 11000;

%really 5-15, when 2 seconds at end of previous trial is factored in
intervalMinMaxSec = [3 13];

%trial on which to flip the tones on a relevant flip day
numFlipTrial = 61;

%__________________________________________________________________________
%------------------DO NOT CHANGE----------------------
%set default values for the function input parameters if they are empty
if nargin<1 || isempty(solOpenDuration)
    %To modulate the amount of fluid delivery, set solenoid open delay in ms
    solOpenDuration = 60;
    fprintf('Default solenoid open duration is %2f ms\n', solOpenDuration);
end

%--------------tone--------------------
%ask through a dialog box which tone is rewarded
rewardedSound = questdlg('Which tone (Hz) is rewarded?', 'Rewarded Tone', string(toneA), string(toneB), string(toneA));
%handle response
switch rewardedSound
    case string(toneA)
        rewardedTone = toneA;
        unrewardedTone = toneB;
    case string(toneB)
        rewardedTone = toneB;
        unrewardedTone = toneA;
end

%tell what day it is through a dialog box
whatDay = questdlg('What type of day is it?', 'Set Day Type', 'Train', 'Flip', 'Test', 'Train');
%handle answer
switch whatDay
    case 'Train'
        whatTrainDay = questdlg('What type of training day is it?', 'Set Day Type', 'No probe', 'Probe', 'No probe');
        switch whatTrainDay
            case 'No probe'
                dayType = 1;
                isTesting = 0;
            case 'Probe'
                dayType = 1.1;
                isTesting = 0;
        end
    case 'Flip'
        dayType = 2;
        isTesting = 0;
    case 'Test'
        dayType = 3;
        isTesting = 1;
end

%-----------initialize vars-------------------
%initialize our globale CACHE variable (to allow master loop to interact with NI card)
NI_INPUT = [];  %empty to start (waiting for card to generate data)

%convert solenoid open time to seconds
solOpenDurationSec = solOpenDuration/1000;

%Hard Code parameters
toneDurationSec = 1.50;
rewardDelaySec = 0;

%directory to write out data
outputDirectory = ['C:' filesep 'Users' filesep 'tadrosslab' filesep 'Dropbox (TadrossLab)' filesep 'BIG_DATA' filesep 'Projects' filesep 'Schultz' filesep 'Data' filesep];

%user input mouse information
userInput = inputdlg({'Cohort', 'Mouse', 'Cage', 'Experimental Group', 'Day type?'}, 'User Input', [1 40]);
trialData.cohort = char(userInput(1, :));
trialData.mouse = char(userInput(2, :));
trialData.cage = char(userInput(3, :));
trialData.group = char(userInput(4, :));
trialData.trialType = char(userInput(5, :));

%initialize overlarge variables to hold the data
%scalar vector of the record times; needs to be a double
trialData.recordTimes = zeros(2, 5000000, 'double');
%Array of time stamps indicating every time a tone was presented with rows 1) frequency of tone. Data points should line up with trialData.recordTimes
trialData.toneValues = zeros(1, 5000000, 'double'); 
%Scalar vector of the voltage output to the solenoid. Data points should line up with trialData.recordTimes
trialData.sucroseDeliveryValues = zeros(1, 5000000, 'double'); 
%Array with rows: 1) IR voltage 2) binary yes/no (1/0) as to whether that measurement counted as a lick while the code was running. Data points should line up with trialData.recordTimes
trialData.lickValues = zeros(2, 5000000, 'double');
%Scalar vector of time stamps indicating every time video data was retrieved. Time stamps are negative if there was an error in video retrieval.
trialData.videoTimes = zeros(1, size(videoData,3));
%Scalar vector with rows (1) Edgecount from rotary encoder. Data points should line up with trialData.
trialData.rotationalSpeed = zeros(1, 5000000, 'double'); 

%overall reference of trials. Array with rows (1) trial number (2) tone frequency (3) 1/0 rewarded or not (4) interval 
%               (5) time started from matlab clock (6) time stim output
%               starts (aka when tone plays/would play) from matlab clock
%               (7) time solenoid closed from NI clock (8) 1/0 probe or not
trialData.trialInfo = [];

%use as a count for how much of the preallocated arrays have been filled, and then to trim at the end before saving
inputCount = 1;

trialNum = 1; %start with first trial
lastLickTime = 0;
lastTrialTime = 0;
numRewardedTrials = 0;
numUnrewardedTrials = 0;
numProbeTrials = 0;
trialData.flipTrial = 0;
FlagStimIsNotOver = 0;
freqToneHz = 0;
saveTime = 420;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------NI card setup----------------------
%create the record session
RecordSession = daq.createSession('ni');
RecordSession.Rate = 1000; %lower resolution for recording than for stimulating
RecordSession.IsContinuous = true;
%set up our inputs
%Change 'Dev1' as necessary for your NI Card's name in the computer
ch1 = addAnalogInputChannel(RecordSession, 'Dev1', 'ai0', 'Voltage'); %lick port detection
ch2 = addCounterInputChannel(RecordSession, 'Dev1', 0, 'EdgeCount');  %rotary encoder
ch3 = addAnalogInputChannel(RecordSession, 'Dev1', 'ai3', 'Voltage'); %tone output measured
ch4 = addAnalogInputChannel(RecordSession, 'Dev1', 'ai2', 'Voltage'); %solenoid output measured
ch5 = addAnalogInputChannel(RecordSession, 'Dev1', 'ai9', 'Voltage'); %read input flag for end of trial

%set channel ranges
ch1.Range = [-5, 5];
ch3.Range = [-5, 5];
ch4.Range = [-5, 5];
ch5.Range = [-5, 5];

%set up listeners
listenerDA = addlistener(RecordSession, 'DataAvailable', @newDataAvailable);
listenerEO = addlistener(RecordSession, 'ErrorOccurred', @newError);

%create the output session
StimSession = daq.createSession('ni');
StimSession.Rate = 100000; %scans per second (100 kHz)
StimSession.IsContinuous = false;
%set up our four outputs
addAnalogOutputChannel(StimSession, 'Dev1', 'ao0', 'Voltage');       %tone
addAnalogOutputChannel(StimSession, 'Dev1', 'ao1', 'Voltage');       %solenoid
addDigitalChannel(StimSession, 'Dev1', 'port0/line2', 'OutputOnly'); %+5V for IR light
addDigitalChannel(StimSession, 'Dev1', 'port0/line0', 'OutputOnly'); %flag for end of trial

%1 second of baseline output data
baselineOutput = zeros(StimSession.Rate,4);
baselineOutput(:, 3) = 1; %set Ir emit to high

%send out one scan background output to set pins to correct settings
outputSingleScan(StimSession, [0 0 1 0]);

%--------------ir setup----------------------
%determine IR baseline for this mouse/lickport setup
%needs to happen before background starts
confirmBaseline = 0;
while ~confirmBaseline
    baselineArray = zeros(1, 500);
    for i = 1:500
       %get 10 readings and determine the minimum
       baselineVoltVal = RecordSession.inputSingleScan;
       baselineArray(1, i) = baselineVoltVal(1);
    end
    baselineVolt = mean(baselineArray(1, 400:500));
    %confirm with user that minimum is okay - if no (ie if mouse licks constantly), recheck
    question = ['Is ' num2str(baselineVolt) ' V an acceptable baseline?'];
    answer = questdlg(question);
    if strcmpi(answer, 'Yes')
        confirmBaseline = 1;
    end
end
outputSingleScan(StimSession, [0 0 1 0]);

%set the threshold as 1.0V over baseline
thresholdLickVal = baselineVolt + 1.0;
trialData.lickThreshold = thresholdLickVal;
trialData.baselineThreshold = baselineVolt;

%--------------------video setup------------------
%also create video capture and VideoLogger objects
imaqreset; %reset to make sure gigE camera is available
videoCapture = videoinput('gige', 'acA1300-60gm');
videoCapture.FramesPerTrigger = 1; %default was 10, but I don't think it matters if you use getsnapshot, so probably no benefit
s = videoCapture.source; %get handle to the camera settings
s.PacketSize = 9014;  %***this is likely the biggest improvement -- default was 1000 or so; so hopefully 9x less likely to have a packet dropped
s.GainRaw = 2;  %this make the video brighter, but shouldn't affect performance
s.AutoFunctionAOIHeight = size(videoData,1);  %here was my attempt to crop on camera, but doesn't seem to have an effect.
s.AutoFunctionAOIWidth = size(videoData,2);
s.CenterX = 'True';
s.CenterY = 'True';
triggerconfig(videoCapture, 'manual');

start(videoCapture);
set(0, 'DefaultFigureWindowStyle', 'docked');
close all; figure(1); clf; colormap(gray); %open the figure
%set(gcf,'CurrentCharacter', ''); %make sure that no keypress is lingering from last time

fprintf('Press any key when ready to start trials\n');
hImage = [];
while isempty(get(gcf,'CurrentCharacter'))
    try
        data = getsnapshot(videoCapture);
        %crop
        rMid = round(size(data,1)/2); cMid = round(size(data,2)/2);
        delR = size(videoData,1)/2;   delC = size(videoData,2)/2;
        data = data(rMid-delR+1:rMid+delR, cMid-delC+1:cMid+delC); %crop
        if isempty(hImage)
            %first time drawing an image
            hImage = image(data); 
            axis image; axis off;
        else
            %faster version - replace image each time
            set(hImage, 'CData', data);
        end
        %drawnow;
    catch
        disp('getsnapshot fail');
        disp(lasterr);
    end
end

%next two lines should be juxtaposed to make tic/toc time approximatley equal to NI time
startBackground(RecordSession);
expClock = tic;  

%Experiment starts now
startTime = clock; 

%__________________________________________________________________________
%--------------------MAIN LOOP-----------------------
while trialNum <= numTrials && toc(expClock) <= numSeconds
    %------------intermediate save code-----
    if toc(expClock) > saveTime
        dbstop if error
        save([outputDirectory filesep trialData.cohort filesep 'Mouse ' trialData.mouse filesep trialData.cohort '_' trialData.mouse '_' datestr(startTime, 'yyyy-mm-dd-HHMM') '_trialData.mat'], 'trialData');
        saveTime = saveTime + 420;
        dbclear all
    end
        
    trialStart = toc(expClock);
    randIntervalSec = randi(intervalMinMaxSec);
    fprintf('Trial number is %3.0f with an interval of %2.2f seconds\n', trialNum, randIntervalSec);
    %while it has been less than randIntervalSec since the last lick
    while ((toc(expClock) - max(lastLickTime, lastTrialTime)) < randIntervalSec)   &&   (toc(expClock) <= numSeconds)
        checkLickSpeed();
    end
    
    %------------flip day code--------------
    %for the first set number of trials, continue with the previous paradigm
    %on the flip trial on a flip day, switch the tones
    if dayType == 2 && trialNum == numFlipTrial
        trialData.flipTrial = numFlipTrial;
        %switch the values for rewardedTone and unrewardedTone
        temp1 = rewardedTone;
        temp2 = unrewardedTone;
        unrewardedTone = temp1;
        rewardedTone = temp2;
        %note that we are now testing
        isTesting = 1;
    end
    %----------probe code-----------------
    probeReady = 0;
    wasProbe = 0;
    if dayType > 1 && dayType < 3
        %have happen every 10 mins on probe training day and flip day
        if toc(expClock) > (numProbeTrials*600) + 600 && numProbeTrials < 6
            if isTesting ~= 1 %no probe trials while testing
                probeReady = 1;
            end
        end
    end
   
    %--------tone/reward code-------------
    if toc(expClock) <= numSeconds        
        x = rand; %random number between 0 and 1 to determine which tone gets played
        %50% probability of rewarded trial; rewarded tone always plays
        if x < 0.5
            bGiveReward = 1;
            freqToneHz = rewardedTone;
            numRewardedTrials = numRewardedTrials + 1;
            disp("Current trial is rewarded");
        %50% probability of unrewarded trial
        else
            bGiveReward = 0;
            numUnrewardedTrials = numUnrewardedTrials + 1;
            %testing, play the unrewarded sound
            if isTesting || probeReady
                freqToneHz = unrewardedTone;
                if probeReady
                    numProbeTrials = numProbeTrials + 1;
                    fprintf("Current trial is probe trial number: %1.0f\n", numProbeTrials);
                    %mark if this trial was a probe trial
                    wasProbe = 1;
                end
            %training, no unrewarded sound
            else
                freqToneHz = 0;
            end
            disp("Current trial is unrewarded");
        end
        
        %generate NI card stimulus
        stimDurSec  = toneDurationSec + rewardDelaySec + solOpenDurationSec; 
        stimDurSamp = (stimDurSec * StimSession.Rate) + 1; 
        tStimSec = ((0:stimDurSamp-1)/StimSession.Rate)';  %time vector in seconds
        STIM = zeros(stimDurSamp, 4);
        STIM(:,3) = 1; %always keep this at +5V
        if bGiveReward
            %indices when solenoid should be open
            nDispense = find(tStimSec>=(toneDurationSec+rewardDelaySec) & tStimSec<(toneDurationSec + rewardDelaySec + solOpenDurationSec));
            %set these indices to +5V
            STIM(nDispense, 2) = 5;
        end
        if freqToneHz>0
            %indices when tone should be playing
            nToneLast = max(find(tStimSec<toneDurationSec));
            while rem(2 * freqToneHz * tStimSec(nToneLast),1) > 1e-6
                nToneLast = nToneLast + 1;
            end
            amplitude = .05;
            STIM(1:nToneLast, 1) = amplitude * sin(2 * pi * freqToneHz * tStimSec(1:nToneLast));  
        end
        %add two second of baseline stim so that any fluctuations when
        %background output ends are not when we are measuring licking
        %baselineOutput is already 1s long
        STIM = [STIM; baselineOutput; baselineOutput];
        STIM(end-1000:end, 4) = 1; %flag the end of the stimulus so the input knows the trial is over

        StimSession.queueOutputData(STIM);
        toneStartTime = toc(expClock);
        startBackground(StimSession);
        FlagStimIsNotOver = 1;
        %need to wait until sucrose is completely delivered before moving on
        while FlagStimIsNotOver
            checkLickSpeed();
        end
        
        %to ensure no attempt to overwrite while background output still happening
        %pause until background session is completely done to prevent
        %overwriting errors - and record the amount of time paused
        tryCount = 0;
        errCount = 0;
        amtPaused = 0;
        while tryCount == errCount
            try
                outputSingleScan(StimSession, [0 0 1 0]);
            catch
                warning('Attempted to overwrite a background session. Pausing for another .1 seconds.');
                errCount = errCount + 1;
                pauseClock = tic;
                while toc(pauseClock) < .1
                    checkLickSpeed();
                end
                amtPaused = amtPaused + .1;
            end
            tryCount = tryCount + 1;
        end
        lastSolTime = lastTrialTime - 2 - amtPaused;
        %update trialInfo for this trial
        trialData.trialInfo = [trialData.trialInfo [trialNum; freqToneHz; bGiveReward; randIntervalSec; trialStart; toneStartTime; lastSolTime; wasProbe]];
    %if the interval is not reached before time runs out, subtract this trial from the total final count
    else
        trialNum = trialNum - 1;
    end

    %if fluid is not licked, don't deliver more
    %only relevant on trials where reward was delivered
    %modify last trial time to reflect addition to output and pausing to ensure no overwriting
    %want mouse to have licked since the fluid was delivered   
    if bGiveReward
        while (lastLickTime < lastSolTime)  &&  (toc(expClock) <= numSeconds)
            checkLickSpeed();
        end
    end
    %start a new trial
    trialNum = trialNum + 1;
end

%end of for loop, check licks for a few seconds to record from the last trial
pauseTime = tic();
while toc(pauseTime) < 5
    checkLickSpeed();
end

dbstop if error
%display total trials actually run
trialNum = trialNum - 1;
fprintf("Total trials run: %3.0f\n", trialNum);
fprintf("Total rewarded trials: %3.0f\n", numRewardedTrials);
fprintf("Total unrewarded trials: %3.0f\n", numUnrewardedTrials);
fprintf("Total probe trials: %3.0f\n", numProbeTrials);

%save misc variables
trialData.numTrials = trialNum;
trialData.numRewardedTrials = numRewardedTrials;
trialData.numUnrewardedTrials = numUnrewardedTrials;
trialData.numProbeTrials = numProbeTrials;
trialData.interval = intervalMinMaxSec;
trialData.solenoidOpenDuration = solOpenDuration;
trialData.rewardDelay = rewardDelaySec;
trialData.toneDurationSec = toneDurationSec;
trialData.rewardedSoundHz = rewardedTone;
trialData.unrewardedSoundHz = unrewardedTone;
trialData.dayType = dayType;
trialData.dateTime = datestr(startTime, 'yyyy-mm-dd-HHMM');
trialData.box = 1;

%trim the scalar vectors
trialData.recordTimes = trialData.recordTimes(:,1:inputCount);
trialData.lickValues = trialData.lickValues(:,1:inputCount);
trialData.rotationalSpeed = trialData.rotationalSpeed(1,1:inputCount);
trialData.toneValues = trialData.toneValues(1, 1:inputCount);
trialData.sucroseDeliveryValues = trialData.sucroseDeliveryValues(1, 1:inputCount);

%trim the videoTimes
trialData.videoTimes = trialData.videoTimes(1:NumVideoFrames);

%save the data
save([outputDirectory filesep trialData.cohort filesep 'Mouse ' trialData.mouse filesep trialData.cohort '_' trialData.mouse '_' datestr(startTime, 'yyyy-mm-dd-HHMM') '_trialData.mat'], 'trialData');

%end NI card session
StimSession.IsDone;
RecordSession.IsDone;
delete(listenerDA);
delete(listenerEO);

dbclear all
%save video
try
    delete(videoCapture); %close to make sure gigE camera doesn't complain
catch
    disp('problem deleting videoCapture');
end
videoFilename = [outputDirectory filesep trialData.cohort filesep 'Mouse ' trialData.mouse filesep trialData.cohort '_' trialData.mouse '_' datestr(startTime, 'yyyy-mm-dd-HHMM') '_video.avi'];
if ~isempty(dir(videoFilename))
    confirmDelete = input('video file already exists; Delete? y/n: ', 's');
    if confirmDelete ~= 'n'
        delete(videoFilename);
        disp('deleted');
    end
end
videoLogger = VideoWriter(videoFilename, 'Motion JPEG AVI');
open(videoLogger);
%GUI waitbar to visualize how much progress has been made in saving the video
h = waitbar(0,['Saving Video (' num2str(NumVideoFrames) ' frames)']);
for k = 1:NumVideoFrames
    waitbar(k/NumVideoFrames)
    writeVideo(videoLogger, videoData(:,:,k));
end
close(h)
close(videoLogger);


%__________________________________________________________________________
%nested function
%----------------------checkLickSpeed--------------------------
%function to 1) check licks 2) get rotary encoder speed 3) grab frame from camera
    function checkLickSpeed
        %if there is new input data to grab
        %only add input data if the input data has been updated
        if ~isempty(NI_INPUT)
            TMPBIG = NI_INPUT;
            NI_INPUT = [];
            
            TMP = TMPBIG;
            
            %invert all variables to have time vary with column, not row
            %copy cache output into variables
            preTime = (TMP(:,1))';
            lickVoltVal =  (TMP(:,2))';
            rpm = (TMP(:,3))';
        
            %read back the outputs
            soundStim = (TMP(:,4))';
            solStim = (TMP(:,5))';
           
            numSamples = length(preTime);
            endOfSamples = inputCount+numSamples-1;
            
            %the time that the cache is retrieved; every 100ms
            postTime = repmat(toc(expClock), 1, size(TMP,1)); %fake it -- same value, but same # of rows as STIM        

            wasLick = lickVoltVal > thresholdLickVal;  %this will do vector boolean
            %update arrays
            trialData.recordTimes(:,inputCount:endOfSamples) = [preTime; postTime];
            trialData.lickValues(:,inputCount:endOfSamples) = [lickVoltVal; wasLick];
            trialData.rotationalSpeed(1,inputCount:endOfSamples) = rpm;
            trialData.toneValues(1,inputCount:endOfSamples) = soundStim;
            trialData.sucroseDeliveryValues(1,inputCount:endOfSamples) = solStim;
            if any(wasLick) %indicates a touch of the lick port occured
                fprintf('lick occured at time %4.2f and was %4.2f volts\n', preTime(find(wasLick,1,'last')), lickVoltVal(find(wasLick,1,'last')));
                lastLickTime = preTime(find(wasLick,1,'last')); %these are all the same value (see "fake it" above)
            end
            
            %update inputCount with the length of the new chunk of data
            inputCount = inputCount + numSamples;
            
            %flag and time for the end of the previous trial, whether
            %rewarded or not
            if FlagStimIsNotOver
                wasTrialOver = (TMPBIG(:, 6))' > 3;
                %check the last chunk of data for completion of stim
                if any(wasTrialOver)
                    preTimeBIG = (TMPBIG(:,1))';
                    FlagStimIsNotOver = 0;  %now it's NOT not over!!! 
                    lastTrialTime =  preTimeBIG(find(wasTrialOver,1,'last'));
                end
            end
            
            global DEBUG;
            DEBUG = [DEBUG; toc(expClock) preTime(1)];
        end
        
        %run video all the time
        tVideoNow = toc(expClock);
        %don't grab a frame too often;  max is every 0.1 sec
        if (NumVideoFrames==0) || (abs(trialData.videoTimes(NumVideoFrames)) < tVideoNow-0.1 )
            %take a new video frame
            NumVideoFrames = NumVideoFrames + 1;
            try
                data = getsnapshot(videoCapture);
                %crop
                rMid = round(size(data,1)/2); cMid = round(size(data,2)/2);
                delR = size(videoData,1)/2;   delC = size(videoData,2)/2;
                data = data(rMid-delR+1:rMid+delR, cMid-delC+1:cMid+delC); %crop
                set(hImage, 'CData', data);%faster version                
                
                trialData.videoTimes(NumVideoFrames) = tVideoNow;
                videoData(:,:,NumVideoFrames) = data; %save
                drawnow;
            catch
                trialData.videoTimes(NumVideoFrames) = -tVideoNow;
                disp('getsnapshot fail');
                disp(lasterr);
                %check number of negative video times
                %if over some amouint of time, restart the camera (takes about 1-2 seconds)
                m = NumVideoFrames;
                while trialData.videoTimes(m) < 0
                    m = m-1;
                end
                %get time between current time and last succesful video grab
                timeOfCameraFailure = tVideoNow - abs(trialData.videoTimes(m));
                %if more than two seconds has passed
                if abs(timeOfCameraFailure) > 2
                    disp('Resetting videoCapture');
                    %reset the camera
                    imaqreset; %reset to make sure gigE camera is available
                    videoCapture = videoinput('gige');
                    triggerconfig(videoCapture, 'manual');
                    start(videoCapture);
                end
            end
        end
     end

%__________________________________________________________________________
%nested function
%--------------functions for the NI card listeners----------------------------
    function newDataAvailable(src,event)      
        %global NI_INPUT
        %push data into our global read cache (the master loop is waiting and checking this periodically)
        NI_INPUT = [NI_INPUT; [event.TimeStamps event.Data]];
    end

    function newError(src,event)
        %output an error
        disp(getReport(event.Error));       
    end
end