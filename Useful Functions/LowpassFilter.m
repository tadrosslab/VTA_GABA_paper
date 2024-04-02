%Written by MRT

% MIT License
% Copyright (c) 2024 Michael R. Tadross

function [lowydata] = LowpassFilter(xdata, ydata ,fcutoff, bPreZeroEnds)
%function [lowydata] = LowpassFilter(xdata, ydata ,fcutoff, bPreZeroEnds)
%
%set fcutoff==0 for no filtering
%
%bPreZeroEnds is an optional parameter, default value is 0
%             if set to 1, then ydata is shifted so that ydata(1) and
%             ydata(end) = 0 prior to filtering, then unshifted

lowydata = ydata;  %allocate memory for return
if fcutoff==0
    return;
end
try
    test= bPreZeroEnds;
catch
    bPreZeroEnds = 0; %default
end

%zeropad the beginning of the trace
Npad = round(length(xdata)/10);
ydatapad = zeros(1,Npad);
dT = (xdata(2)-xdata(1));		% sampling interval in sec
fmax = 1/dT;			  	% this is the maximum freq (hz)

for k=1:size(ydata,1)
    
    %MRT - prior to filtering, set ends of trace to 0 by subtracting line
    %this avoids discontinuities at the ends
    if bPreZeroEnds
        NumSampEnds = round(2*fmax/fcutoff);
        NumSampEnds = max(NumSampEnds, 1);  %at least 1
        NumSampEnds = min(NumSampEnds, floor(length(xdata)/2)); %no more than half the data
        P = polyfit(xdata([1:NumSampEnds end-NumSampEnds+1:end]), ydata(k,[1:NumSampEnds end-NumSampEnds+1:end]), 1);
        ydata(k,:) = ydata(k,:) - polyval(P, xdata);
    end

    %Start PGP code
    ydata1row = [ydatapad,ydata(k,:)];

    %calculate inputs to GAUSS filter
    N = length(xdata) + Npad;		% number of points
    N2=2^ceil(log(N)/log(2));

    % pgp's GAUSS filter.
    wc = 2*pi*fcutoff;
    w=0:2*pi/(N2*dT):2*pi*(N2-1)/(N2*dT);				% radians/sec 0:2*pi/dT
    w(N2/2+2:N2)=-2*pi/(2*dT)+2*pi/(N2*dT):2*pi/(N2*dT):-2*pi/(N2*dT);	% Adjust so +/- 2*pi/2T
    lfilt=GAUSS(w,wc,dT);

    %n = round(N2*fcutoff/fmax)			% end of square pulse
    %lfilt = zeros(1,N2);
    %lfilt(1:n) = ones(1,n);
    %lfilt(N2-n+1:N2) = ones(1,n);		% this is the crude filter

    %ydata = get(kids(j),'ydata');
    ydata1row = real(ifft(fft(ydata1row,N2).*lfilt));
    lowydata(k,:) = ydata1row(1+Npad:length(xdata)+Npad);
    
    %MRT - undo the zeroing
    if bPreZeroEnds
        lowydata(k,:) = lowydata(k,:) + polyval(P, xdata);
    end
end

function H = GAUSS(w,wc,dT)
% GAUSSIAN Filter Components with corner (3dB) frequency at wc.

sigma = 0.1325*(2*pi/wc);

W=(-sigma*sigma/2)*w.^2;
H=exp(W);

% Gaussian Filters have no delay

