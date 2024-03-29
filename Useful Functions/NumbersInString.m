%Written by MRT

% MIT License
% Copyright (c) 2024 Michael R. Tadross

function [Numbers, Starts, Stops] = NumbersInString(Str)
%function [Numbers, Starts, Stops] = NumbersInString(Str)
%
%finds a list of numbers embedded in a string
%returns a vector of the numbers = Numbers (1 x N, where N = # of numbers found)
%also returns the indices of these in the string = Starts, Stops (1 x N).

%default returns
Numbers = [];
Starts = [];
Stops = [];

k=1;
while k<=length(Str)
    %find the first digit of a number
    while k<length(Str) && ~IsReallyReal(str2double(Str(k))) 
        k = k+1;
    end

    %we are done if we haven't found a number
    if ~IsReallyReal(str2double(Str(k)))
        return;
    end
    
    %now find the span of characters that are numeric
    j=k;
    while j<length(Str) && IsReallyReal(str2double(Str(k:j+1)))
        j=j+1;
    end

    %make sure there wasn't a +/- or decimal point at the start of the number
    while k>1 && IsReallyReal(str2double(Str(k-1:j)))
        k=k-1;
    end
    
    %if the number ends in a decimal point, then really it ended before the decimal
    %and it's actually just a perioed, as in the filename  e2.mrt
    if Str(j)=='.'
        j=j-1;
    end
    
    %add to the list
    Numbers(end+1)      = str2double(Str(k:j));
    Starts(end+1) = k;
    Stops(end+1) =  j;

    k=j+1;
end


function bIs = IsReallyReal(Num)

bIs = isreal(Num) && ~isnan(Num);

