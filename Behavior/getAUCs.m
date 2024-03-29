%Written 3/2024 by SB

% MIT License
% Copyright (c) 2024 Michael R. Tadross


%load _analysis_tSchultz into the workspace - this is the p variable

%flip happens at 225, in our resized data structures made in
%licksAcrossMiceWithTreadmill
postFlipIdx = [226 400];

%get area under curve of the normalized data p(k).norm
%first row is tone B conditioning
conditioningAUC = trapz(p(1).norm(:, postFlipIdx(1):postFlipIdx(2)), 2);
%second row is tone A extinction
extinctionAUC = trapz(p(2).norm(:, postFlipIdx(1):postFlipIdx(2)), 2);

%as there are 100 "x values" per day, divide by 100 to get the AUC unit in
%behavioral sessions
conditioningAUC = conditioningAUC/100;
extinctionAUC = extinctionAUC/100;

