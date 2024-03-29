%Written 4/2022 by SB and MRT
%Cleaned up and commented 3/2024 by SB
% MIT License
% Copyright (c) 2024 Michael R. Tadross

%function to plot cumulative spikes and firing rate of templates extracted from spyking circus

%Instructions:
%   run in the folder with the result-v1 sorted data extracted from Spyking Circus
%   input: chWithCell = array with the channel of interest per template, in order of the templates
%   output: Templates_data structure
%           saves some calculated features and also the full spikeTimesSec variable for future analysis

function plotSpikes(chWithCell)

set_Nt = 3; %modify this to fit the Nt used while sorting

D = dir('*result-v1*');
fName = D(1).name;

%display the particular template
tmpD = dir('*templates-v1*');
tmpName = tmpD(1).name;

set(0,'DefaultFigureWindowStyle','docked');
figure(1); clf;
figure(2); clf;
figure(3); clf;
LegTxt = {''};

%---Extract the templates---
%this matrix has a size that is twice the number of templates 2k
%Only the first k elements are the real templates.
templates_size = double(h5read(tmpName, '/temp_shape'));
N_e = templates_size(1); %number of electrodes
N_t = templates_size(2); %temporal width of templates
N_temps = templates_size(3)/2;

temp_x = double(h5read(tmpName, '/temp_x') + 1);
temp_y = double(h5read(tmpName, '/temp_y') + 1);
temp_z = double(h5read(tmpName, '/temp_data'));
templates = sparse(temp_x, temp_y, temp_z, N_e*N_t, templates_size(3)); 

temps_data = [];

%---Initialize the data structure---
for t = 1:N_temps
    Templates_data(t).channel = chWithCell(t);
    Templates_data(t).spikeTimesSec = [];
    Templates_data(t).ISIs = [];
    Templates_data(t).firingRates = [];
    Templates_data(t).meanISIs = [];
    Templates_data(t).tempWidth = [];
    Templates_data(t).tempMin = [];
end

for k = 1:N_temps
    spikes = double(h5read(fName, ['/spiketimes/temp_' num2str(k-1)])); %times at which spikes are occurring)
    amps = double(h5read(fName, ['/amplitudes/temp_' num2str(k-1)]));

    amps=amps(1,:);
    %to calculate the real times at which spikes are occurring, divide
    %spikes by the recording sampling rate in Hz
    SpikeTimeSec = spikes/30000; 
    
    %code to eliminate duplicate times
    for q=2:length(SpikeTimeSec)
        if SpikeTimeSec(q)<=SpikeTimeSec(q-1)
            SpikeTimeSec(q)=SpikeTimeSec(q-1)+1e-5;
        end
    end
    Templates_data(k).spikeTimesSec = SpikeTimeSec;
    CumSpikes    = 1:length(spikes); %number of cumulative spikes
    
    %resample SpikeTimeSec into firing rates
    ResampTimeSec  = 0:1:max(SpikeTimeSec);
    ResampCumSpike = interp1(SpikeTimeSec, CumSpikes, ResampTimeSec);
    ResampAmp      = interp1(SpikeTimeSec, amps,  ResampTimeSec);
   
    %rate estimates
    FiringRate    = diff(ResampCumSpike)./diff(ResampTimeSec); %this variable length is the number of seconds in the recording
    RateTimeSec   =  (ResampTimeSec(1:end-1) + ResampTimeSec(2:end))/2;
    
    SM = 10; %smooth window (in seconds)
    LegTxt{k} = ['channel ' num2str(chWithCell(k))];
    
    %Plot to visualize amplitude, cumulative spikes, and firing rate over time per cell
    figure(1); subplot(3, 1, 1); 
    plot(ResampTimeSec/60, smooth(ResampAmp,SM), '-'); grid on; hold on; ylim([0 2]);
    xlabel('Time (minutes)');
    ylabel('Spike amplitude');
    legend(LegTxt);
    
    figure(1);subplot(3, 1, 2);
    plot(ResampTimeSec/60, ResampCumSpike, '-'); grid on; hold on;
    xlabel('Time (minutes)');
    ylabel('Cumulative spikes');
    legend(LegTxt);
    
    figure(1);subplot(3, 1, 3);
    plot(RateTimeSec/60, smooth(FiringRate,SM), '-', 'linewidth', 2); grid on; hold on;
    xlabel('Time (minutes)');
    ylabel('Firing Rate (Hz)');
    legend(LegTxt);
    

    %----Firing Rate-------%
    overallFR = mean(FiringRate, 'omitnan');
    baselineFR = mean(FiringRate(1:900), 'omitnan'); %first 15 min, or 900s, is the baseline window
    baselinestd = std(FiringRate(1:900), 'omitnan');
    Templates_data(k).firingRates = [overallFR baselineFR];
    
    
    %----Template width and plotting--------%
    figure(k+1); clf; subplot(2, 1, 1);
    temp_i = full(reshape(templates(:, k), N_t, N_e));
    tempWithCell = temp_i(:, chWithCell(k));
    x=0:(set_Nt/N_t):set_Nt;
    x = x(2:N_t+1);
    plot(x, tempWithCell); hold on;
    title(['Channel ' num2str(chWithCell(k)) ' template']);
    xlabel('Time (ms)');
    ylabel('Amplitude (mV)');
   
    %calculate template width: duration of temp_i from spike initiation > .25 for the first time to the maximal negative component (min)
    peak = min(tempWithCell);
    temp_duration = ((find(tempWithCell==peak)-find(tempWithCell > .25,1))/61)*N_t; %this is in ms    
    temps_data = [temps_data; peak temp_duration];
    Templates_data(k).tempWidth = temp_duration;
    Templates_data(k).tempMin = peak;
    %sanity check option: plot raw data and see if these ms match up with the spike width in real time

    
    %---------ISI calculations-------------%
    %calculate ISI across whole session
    isis = diff(SpikeTimeSec);
    overallisi = mean(isis, 'omitnan');
    Templates_data(k).ISIs = isis;
    
    %15 mins (900s) baseline prior to infusion
    baselineISI = diff(SpikeTimeSec(SpikeTimeSec<900));
    baselineISI = rmmissing(baselineISI);
    baselinemean = mean(baselineISI);
    
    %15 min bin near the end of the recording  
    infusionISI = diff(SpikeTimeSec(SpikeTimeSec>5400 & SpikeTimeSec<6300));
    infusionISI = rmmissing(infusionISI); 
    infusionmean = mean(infusionISI);   
    
    Templates_data(k).meanISIs = [overallisi baselinemean infusionmean];
  
    % histogram of log of ISIs before and after infusion 
    % log is used to pull out smaller ISI intervals
    % use to visualize any changes in activity between baseline and infusion 
    bins = 0:.0025:1.5;
    figure(k+1); subplot(2, 1, 2);
    h = histogram(baselineISI, bins); hold on;
    histogram(infusionISI, bins); hold on;
    title(['Channel ' num2str(chWithCell(k)) ' interspike intervals']);
    xlabel('ISIs (s)');
    ylabel('count');    
    
end

%plot template width versus minimum peak per cell to help with decision of whether a cell
%is dopamine or not
figure(k+2); clf;
scatter(temps_data(:, 2), temps_data(:, 1), 'filled', 'o');
xlabel('Template width');
ylabel('Minimum peak');

%save Templates_data
path = cd;
filename = [path '/Templates_data'];
save(filename, 'Templates_data');

end
