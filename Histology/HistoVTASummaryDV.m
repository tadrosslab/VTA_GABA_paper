%Written 2018 by MRT, edited 2021 by SB
%Commented and cleaned up 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%Function to summarize quantification performed by HistoOfflineAnalysis
%Instructions:
%   Run in a folder with HistoROI_Analysis.mat
%   Output: result of HistoOfflineAnalysis; summarized quantification as variables in the workspace
%       SOI: Values in summed ROIs
%       SOIDV: Values in summed ROI rows (3 rows: dorsal, medial, and ventral)
%       In both:
%           H1 and H2 first column = 'virus' from HistoOfflineAnalysis (ch 2)
%           H1 and H2 second column = 'dye' from HistoOfflineAnalysis (ch 1)

function [SOIDV, SOI] = HistoVTASummaryDV(OUTPUT)

load('HistoROI_Analysis.mat');
%extract background fluorescence, averaged over every section to get a
%whole brain value
[ByMouse.BKG, ~, NumBlobs] = BkgFromHist(sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3), OUTPUT.binIntensity, 1);
if NumBlobs~=1
    warning('BkgFromHist NumBlobs should be 1 for this Fix hack to work');
end
NSec = length(OUTPUT.Bkg);

%Zero-Out any data where Red < 1950;
VirusThresh = 1950
for h=1:2
    for k=1:NSec
        for r=1:size(OUTPUT.Hemisphere(h).Block,1)
            for c=1:size(OUTPUT.Hemisphere(h).Block,2)
                OUTPUT.Hemisphere(h).Block(r,c,k).count2D(OUTPUT.binIntensity<VirusThresh,:) = 0;
            end
        end
    end
end

%get the sums of each hemisphere
[ByMouse.H1, ByMouse.Pix1]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
[ByMouse.H2, ByMouse.Pix2]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);

NUMROW = size(OUTPUT.Hemisphere(1).Block,1); %number of dorsal-->ventral subdividions

BySec.H1  = zeros(NSec, 2);
BySec.H2  = zeros(NSec, 2);
BySec.Pix1  = zeros(NSec, 1);
BySec.Pix2  = zeros(NSec, 1);

BySecDV.H1  = zeros(NUMROW, 2, NSec);
BySecDV.H2  = zeros(NUMROW, 2, NSec);

%go through each section
for k=1:NSec
    BySec.BKG = ByMouse.BKG; %copy for completeness
    %row, column, section
    [BySec.H1(k,:), BySec.Pix1(k,:)]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    [BySec.H2(k,:), BySec.Pix2(k,:)]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    for nR=1:NUMROW
        %column, section
        BySecDV.H1(nR,:,k)  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(nR,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
        BySecDV.H2(nR,:,k)  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(nR,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    end
end

%find the strongest 15-sections (using maximum dye capture, which is column 2)
DyePerSection = BySec.H1(:,2) + BySec.H2(:,2);
WIDTH = 15; %number of sections for VTA quantification
for k=1:NSec-WIDTH+1
    TMP(k) = sum(DyePerSection(k:k+WIDTH-1));
end
[~,kBest] = max(TMP);

%Sections Of Interest - summed over all rows
SOI.BKG = ByMouse.BKG; %copy for completeness
SOI.SecIndx = kBest:kBest+WIDTH-1;
SOI.H1   = sum(BySec.H1(SOI.SecIndx,:), 1); %summed across sections SOI.SecIndx, column 1 is virus, column 2 is dye
SOI.H2   = sum(BySec.H2(SOI.SecIndx,:), 1); %summed across sections SOI.SecIndx, column 1 is virus, column 2 is dye

%Sections of Interest - summed over each row in a dorsal-ventral gradient
SOIDV.BKG = ByMouse.BKG; %copy for completeness
SOIDV.SecIndx = kBest:kBest+WIDTH-1;
SOIDV.H1   = sum(BySecDV.H1(:,:,SOI.SecIndx), 3); %summed across sections SOI.SecIndx and row, column 1 is virus, column 2 is dye
SOIDV.H2   = sum(BySecDV.H2(:,:,SOI.SecIndx), 3); %summed across sections SOI.SecIndx and row, column 1 is virus, column 2 is dye

set(0,'DefaultFigureWindowStyle','docked')

end


%subfunctions to run HistoVTASummaryDV
function [BKG, PixelIdxList, NumBlobs] = BkgFromHist(Hist2D, binIntensity, bDebug)
%function [BKG, PixelIdxList, NumBlobs] = BkgFromHist(Hist2D, binIntensity, bDebug)
%
%Calculates tissue fluroescence from Hist2D
%
%INPUTS:
%   Hist2D: 2D histogram; rows=virus; cols=dye; value=frequency.  Note: assume square (same #rows and columns)
%   binIntensity:  intensity value corresponding to each bin (assume identical for virus and dye)
%   bDebug: (optional)  set=1 to show segmentation (=0 by default, if not provided)
%
%OUTPUTS
%   BKG = [virusBKG  dyeBKG]

try
    bDebug; %see if user specified
catch
    bDebug = 0; %set flag to show the segmentation
end
binIntensity = binIntensity(:); %make this a column vector (needed later)

%Segmentation:
ThreshImage = Hist2D/max(max(Hist2D));  %normalize 0 to 1
ThreshImage = im2bw(ThreshImage, 0.1);  %thresholded version (binary)
Blobs  = regionprops(ThreshImage, {'PixelIdxList' 'FilledArea' 'Centroid' }); %segment the blobs
NumBlobs = length(Blobs);
if NumBlobs~=2
    %we expect two blobs --  glass and tissue
    if length(Blobs)==1
        warning('Only one dominant Background intensity ... expected only two (glass & tissue) ... be careful.');
    else
        warning('More than two dominant Background intensities ... expected only two (glass & tissue) ... be careful.');
    end
    bDebug = 1;
end
if bDebug
    figure(1);
    imagesc(log10(Hist2D)); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('background');
    figure(2);
    imagesc(ThreshImage); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('background');
end
%use binCentroids to figure out which blob is tissue (higher value centroid) and which is glass (lower value centroid)
binCentroids = mean(vertcat(Blobs.Centroid),2);
[~,I] = max(binCentroids);

%next, calculate the average intensities for each bin in the 2D-histrogram. For this, we need to know
%how often this bin occurrs (frequency), and its corresponding intensity of virus (row),and dye (col)
PixelIdxList = Blobs(I).PixelIdxList;  %bins in 2D histogram corresponding to tissue
Frequency = Hist2D(PixelIdxList);      %frequency of each bin
[row,col] = ind2sub(size(Hist2D), PixelIdxList);  %intensity of virus(row) and dye (col) for each bin

%calculate means, using histogram weighted averaging
MeanRowVirus = sum(Frequency.*binIntensity(row))/sum(Frequency);
MeanColDye = sum(Frequency.*binIntensity(col))/sum(Frequency);

%output
BKG = [MeanRowVirus  MeanColDye];
end


function [SUM, NumPix] = SumFromHist(Hist2D, binIntensity, BKG)
%function [SUM, NumPix] = SumFromHist(Hist2D, binIntensity, BKG)
%
%Calculates tissue fluroescence from Hist2D
%
%INPUTS:
%   Hist2D: 2D histogram; rows=virus; cols=dye; value=frequency.  Note: assume square (same #rows and columns)
%   binIntensity:  intensity value corresponding to each bin (assume identical for virus and dye)
%   BKG = [virusBKG  dyeBKG]  (see BkgFromHist)
%
%OUTPUTS
%   SUM = [virusSUM dyeSUM]

%this takes care of the background subtraction
binIntensity_Virus = binIntensity - BKG(1);
binIntensity_Dye   = binIntensity - BKG(2);

%this gets the marginal distributions
Hist1D_Virus = sum(Hist2D, 2);  %Virus=dimension 1; sum over the other dimension (2)
Hist1D_Dye   = sum(Hist2D, 1);  %Dye=dimension 2; sum over the other dimension (1)

%calculate means, using histogram weighted averaging
virusSUM = sum(Hist1D_Virus(:).*binIntensity_Virus(:));
dyeSUM = sum(Hist1D_Dye(:).*binIntensity_Dye(:));
NumPix = sum(Hist1D_Dye(:));

SUM = [virusSUM dyeSUM];

end