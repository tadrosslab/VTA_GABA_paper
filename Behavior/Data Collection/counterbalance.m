%Written by SB
%function to give the order to run the mice for each day of habituation and behavior
%saves a matrix in the cohort folder

function counterbalance(cohort, numMice)

outputDirectory = cd;
numDays = 14;
mouseOrder = [];

for i = 1:numDays
    mouseOrder = [mouseOrder; randperm(numMice)];
end

save([outputDirectory filesep cohort '_counterbalance_order'], 'mouseOrder');
end