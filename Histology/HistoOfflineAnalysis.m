%Written 2018 by MRT
%Cleaned up and commented 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross

%Loads images and ROIs from HistoROI.mat and quantifies fluorescence of
%specified channels in binned rows and columns per ROI
%Instructions:
%   Run in a folder with ROIs drawn in a saved HistoROIs.mat and the corresponding images
%   Output: HistoROI_Analysis.mat

function HistoOfflineAnalysis

%set parameters
Mag = 1;  %=1 is the full spatial resolution;

%channels - change as necessary
virus = 2; %dTomato channel = viral expression
dye   = 1; %Alexa647 channel = dye capture

%spatially divide into NxM grid
NROW = 3;  %dorsal to ventral
NCOL = 4;  %medial to lateral

SAVEFILE = [pwd '\HistoROI_Analysis.mat'];

if ~isfile('HistoROI.mat')
    disp('No ROI file')
end

%load ROIs
HISTO = load('HistoROI.mat');
HISTO = HISTO(1).HISTO;
HISTO.RootDir = pwd;

TotalSection = 0;
for k=1:length(HISTO.SubDir)
    TotalSection = TotalSection + length(HISTO.SubDir(k).RoiCell);
end

OUTPUT = [];

%Initialize Histogram
OUTPUT.ch = [virus dye];
gamma=10; E = 65536*((0:500)/500).^(gamma); E=unique(round(E));  %generates ~300 bins, with high resolution for low intensity (gamma = 10)
OUTPUT.binEdges = E;  %intensity resolution
OUTPUT.binIntensity = 0.5*OUTPUT.binEdges(1:end-1) +  0.5*OUTPUT.binEdges(2:end);  %intensity of each bin.
OUTPUT.NROW = NROW;
OUTPUT.NCOL = NCOL;


%dock new figure windows
set(0,'DefaultFigureWindowStyle','docked')

nSection = 0;
%loop through each brain section
for k=1:length(HISTO.SubDir)
    SectionROI = length(HISTO.SubDir(k).RoiCell);  %number of ROI pairs for this section (can be 0, or 1 pair;  unlikely to be more than 1 pair)
    if SectionROI>0
        %load the raw data at full resolotion
        [imageData] = bfopenIX83([HISTO.RootDir '\' HISTO.SubDir(k).Name], Mag);
        try
            if HISTO.bVS200
                imageData = rot90(imageData, -1); %counterclockwize by 90
            end
        catch
        end
        nCol = size(imageData,2);
        
        %show the image (DAPI channel, which is always last), just for purposes of visual reassurance
        close all;
        imagesc(imageData(:,:,end));
        set(gca,'ColorMap', gray);
        set(gca,'ydir','reverse'); axis image; hold on;
        title(['Sample #'  num2str(k)  '  of  '  num2str(length(HISTO.SubDir))]);
        drawnow;
        
        BW_BKG = [];
    end
    
    %loop through the ROI pairs in this section
    for j=1:SectionROI
        nSection = nSection + 1;
        disp(['Analyzing Image ' num2str(k) ' of ' num2str(length(HISTO.SubDir)) '.  Section#' num2str(nSection) ' of ' num2str(TotalSection)]);
        
        Fields = {'Nuc'  'Cyto'}; %These names are a historical artifact.  Here, LEFT  = Nuc;  RIGHT = Cyto
        for nHemisphere=1:length(Fields)
            ROI = HISTO.SubDir(k).RoiCell(j).(Fields{nHemisphere});  %first LEFT, then RIGHT
            ROI.Pos = ROI.Pos/(2^Mag);  %convert to current magnification
            
            %draw the ROI and get the ROI mask (pixels inside)
            tmpH = feval(ROI.Type, gca, ROI.Pos);
            drawnow;
            BW = createMask(tmpH);
            drawnow;
            if isempty(BW_BKG)
                BW_BKG = (~BW);  %pixels not in this ROI
            else
                BW_BKG = (~BW) & BW_BKG;  %pixels not in any ROI so far
            end
            
            %section up the ROI into NROW x NCOL blocks
            TMP = find(sum(BW,1));
            ColEdges = round(min(TMP):((1+max(TMP)-min(TMP))/NCOL):1+max(TMP));
            TMP = find(sum(BW,2));
            RowEdges = round(min(TMP):((1+max(TMP)-min(TMP))/NROW):1+max(TMP));
            
            %determine if we are on the LHS or RHS of image
            %show symbol on dorso-medial corner.
            if ColEdges(1)-1 > nCol-ColEdges(end)
                bFlipCol = 0;  %we are on the right side, so naturally going from medial to lateral
                plot(ColEdges(1), RowEdges(1), 'ow', 'linewidth', 2);
            else
                bFlipCol = 1;  %we are on the left side, so need to flip columns to go from medial to lateral
                plot(ColEdges(end), RowEdges(1), 'ow', 'linewidth', 2);
            end
            
            %plot grid
            for nC = 1:length(ColEdges)
                plot(ColEdges([nC nC]), RowEdges([1 end]), ':w');
            end
            for nR = 1:length(RowEdges)
                plot(ColEdges([1 end]), RowEdges([nR nR]), ':w');
            end
            
            %Get the fancy new curvy grid - fits the ROIs better than a square grid
            drawnow;
            ROIset = ROI_Subdivide(ROI, NROW, NCOL);
            if bFlipCol
                %left of image, need to change columns order to go from medial-->lateral
                ROIset = ROIset(:, end:-1:1);
            end
            
            %go through all rows and all columns
            for nR = 1:length(RowEdges)-1
                for nC = 1:length(ColEdges)-1
                    
                    tmpH = feval(ROIset(nR,nC).Type, gca, ROIset(nR,nC).Pos);
                    setColor(tmpH,'r')
                    drawnow;
                    BW2 = createMask(tmpH);
                    drawnow;
                    
                    %pixels in curvy grid (also make sure they are in the original ROI)
                    I = find( BW & BW2  );
                    
                    Pixels = zeros(length(I), 2);
                    for q = 1:length(OUTPUT.ch) %for the channels we care about
                        tmp = imageData(:,:,OUTPUT.ch(q)); %full image, single color
                        Pixels(:,q) = tmp(I);
                    end
                    
                    OUTPUT.Hemisphere(nHemisphere).Block(nR, nC, nSection).ROI = ROIset(nR,nC);
                    OUTPUT.Hemisphere(nHemisphere).Block(nR, nC, nSection).count2D = histcounts2(Pixels(:,1), Pixels(:,2), OUTPUT.binEdges, OUTPUT.binEdges);
                end
            end
        end  %finished one Section/Hemisphere
        
        %BKG histogram for this section:
        I = find(BW_BKG);
        
        Pixels = zeros(length(I), 2);
        for q = 1:length(OUTPUT.ch) %for the channels we care about
            tmp = imageData(:,:,OUTPUT.ch(q)); %full image, single color
            Pixels(:,q) = tmp(I);
        end
        OUTPUT.Bkg(nSection).count2D = histcounts2(Pixels(:,1), Pixels(:,2), OUTPUT.binEdges, OUTPUT.binEdges);
        
    end
end
%save
save(SAVEFILE, 'OUTPUT', '-v7.3');

%%%nested function to subdivide the drawn ROIs into NROWs and NCOLs%%%
    function ROIset = ROI_Subdivide(ROI, NROW, NCOL)
        %function ROIset = ROIsubdivide(ROI, NROW, NCOL)
        %
        %ROI = an ROIpoly -- really should be polygonal (but code can be easily adapted
        %NROW = #rows
        %NCOL = #columns
        
        if ~strcmpi(ROI.Type, 'impoly')
            error('ROIsubdivide is not yet compatible with ROI Type other than impoly');
        end
        
        %convert to an xy list that is resampled
        NSAMP = 4000*NROW*NCOL;  %start with very dense oversampling -- this will be downsampled later
        xy = curvspace(ROI.Pos([1:end 1],:), NSAMP+1); %resample
        [x, y] = poly2cw(xy(:,1),xy(:,2));
        xy = [x y]; %make sure this is counter-clockwise (using inverted y-axis convention)
        xy(end,:) = [];  %remove repeated point
        
        %find the best corners
        %left top
        %left bottom
        %right bottm
        %right top
        OuterCorners = [ ...
            min(xy(:,1))  min(xy(:,2)); ...
            min(xy(:,1))  max(xy(:,2)); ...
            max(xy(:,1))  max(xy(:,2)); ...
            max(xy(:,1))  min(xy(:,2)); ...
            ];
        %find the closest point along path to each of the four corners
        for c=1:4
            xyC = repmat(OuterCorners(c,:),size(xy,1),1); %desired corner
            DIST = (xy - xyC).^2; %[(x-xC).^2  (y-YC).^2]
            DIST = sqrt(sum(DIST, 2));
            [~,kBest(c)] = min(DIST); %index of closest point
        end
        %circ permute to first kBest
        xyTMP = xy(mod((1:NSAMP)+kBest(1)-2,NSAMP)+1, :);
        kBest = mod(kBest-kBest(1),NSAMP)+1;
        kBest(5) = size(xyTMP,1)+1;  %append last point
        
        %resample each segment so that top/bottom/left/right all have the same number of samples
        xy = [];
        for c=1:4
            xySEG = curvspace(xyTMP(kBest(c):kBest(c+1)-1,:), NSAMP/4);
            xy = [xy; xySEG];
        end
        
        %get and plot the four corners
        Corners = xy(1:NSAMP/4:NSAMP, :);
        plot(Corners(:,1), Corners(:,2), 'ow', 'linewidth', 5);
        
        %next step is extracting the vertical & horizontal outlines
        xyV1 = xy(1+0*NSAMP/4:1+1*NSAMP/4, :);
        xyV2 = xy(1+2*NSAMP/4:1+3*NSAMP/4, :); xyV2=xyV2(end:-1:1,:); %reverse order
        
        xyH2 = xy(1+1*NSAMP/4:1+2*NSAMP/4, :);
        xyH1 = xy([1+3*NSAMP/4:4*NSAMP/4 1], :); xyH1=xyH1(end:-1:1,:); %reverse order
        
        %calculating weighted averages of vertical outlines (e.g., 25%L/75%R etc.)
        for r=0:NROW
            w1 = (NROW-r)/NROW;
            w2 =        r/NROW;
            TMP = w1*xyH1 + w2*xyH2;  %weighted sum of boundary lines
            
            %stretch/rotate/transform, etc.
            Ends0 = TMP([1 end], :);
            Ends1 = [xyV1(1+r*NSAMP/4/NROW,:); xyV2(1+r*NSAMP/4/NROW,:)];
            [~,~,P] = procrustes(Ends1,Ends0,'reflection',false);
            TMP = P.b*TMP*P.T + repmat(P.c(1,:),size(TMP,1),1);
            
            %store & plot each horizontal line
            xyH{r+1} = TMP;
            plot(TMP(:,1), TMP(:,2), 'c', 'linewidth', 2);
        end
        
        for c=0:NCOL
            w1 = (NCOL-c)/NCOL;
            w2 =        c/NCOL;
            TMP = w1*xyV1 + w2*xyV2;  %weighted sum of boundary lines
            
            %stretch/rotate/transform, etc.
            Ends0 = TMP([1 end], :);
            Ends1 = [xyH1(1+c*NSAMP/4/NCOL,:); xyH2(1+c*NSAMP/4/NCOL,:)];
            [~,~,P] = procrustes(Ends1,Ends0,'reflection',false);
            TMP = P.b*TMP*P.T + repmat(P.c(1,:),size(TMP,1),1);
            
            %store & plot each vertical line
            xyV{c+1} = TMP;
            plot(TMP(:,1), TMP(:,2), 'c', 'linewidth', 2);
        end
        
        %determine intersecitons of the curvy lines
        for r=0:NROW
            tmpH = xyH{r+1};  %this is the horizontal line corresponding to ROW #r
            for c=0:NCOL
                tmpV = xyV{c+1};  %this is the vertical line corresponding to COL #c
                
                if r==0 || c==0 || r==NROW || c==NCOL
                    %boundary -- these are always perfect
                    VxH{r+1,c+1} = 1 + r*NSAMP/4/NROW; %index within vertical line near intersection
                    HxV{r+1,c+1} = 1 + c*NSAMP/4/NCOL; %index within horizontal line near intersection
                    xy0{r+1,c+1} = tmpV(VxH{r+1,c+1},:); %= tmpH(HxV,:); %always the same
                else
                    %use fancy function from Douglas Schwarz
                    [x0,y0,VxH{r+1,c+1},HxV{r+1,c+1}] = intersections(tmpV(:,1),tmpV(:,2),tmpH(:,1),tmpH(:,2));
                    xy0{r+1,c+1} = [x0 y0];
                end
                %plot the intersection
                plot(xy0{r+1,c+1}(1), xy0{r+1,c+1}(2), 'oc', 'linewidth', 2);
                
                if r>0 && c>0
                    ROIset(r,c) = ROI;
                    
                    %assemble ROI for line segments (this was tricky, but I got it).
                    ROIset(r,c).Pos = vertcat( ...
                        xy0{r,c}, ...
                        xyV{c}(ceil(VxH{r,c}):floor(VxH{r+1,c}),:), ...
                        xy0{r+1,c}, ...
                        xyH{r+1}(ceil(HxV{r+1,c}):floor(HxV{r+1,c+1}),:), ...
                        xy0{r+1,c+1}, ...
                        xyV{c+1}(floor(VxH{r+1,c+1}):-1:ceil(VxH{r,c+1}),:), ...
                        xy0{r,c+1}, ...
                        xyH{r}(floor(HxV{r,c+1}):-1:ceil(HxV{r,c}),:) ...
                        );
                    %eliminate redundant points
                    ROIset(r,c).Pos = unique(ROIset(r,c).Pos,'rows','stable');
                    
                    %simplify
                    ROIset(r,c).Pos = DouglasPeucker(ROIset(r,c).Pos, 0.001, 0);
                end
            end
        end
        
    end
end


