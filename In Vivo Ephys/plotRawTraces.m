%Written 05/2023 by SB
%Cleaned up and commented 3/2024 by SB
% MIT License
% Copyright (c) 2024 Michael R. Tadross

%function to plot raw traces from open ephys unsorted continuous data

%Instructions:
%   run in folder with raw continuous.dat data
%   input: integer value of the electrode waveform you want to plot
%   output: a plot

function plotRawTraces(electrode)
numelectrodes = 1;

currentfolder = pwd;
datafile = [currentfolder filesep 'continuous.dat'];

FID_shank1_read = fopen(datafile);
fseek(FID_shank1_read, (electrode-1)*2, 'bof');

precision = '1*int16';
skip = (16-numelectrodes)*2;

shank1_data = fread(FID_shank1_read, [1 inf], precision, skip);
plot(shank1_data);
end