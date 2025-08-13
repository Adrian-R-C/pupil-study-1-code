function [outputSignal, outlierIndexes, XMEDIAN, XSIGMA] = nHampel(inputSignal,adjacentSampleHamp, NSIGMA, nTimesHampel)
% nHampel detects and interpolates outliers in time series
% data using Hampel filtering. This function iterates multiple time, as specified by the user. This function identifies
% rapid signal changes that exceed statistical thresholds, marks surrounding
% points as outliers, and replaces them with smoothed values. 
%
% See Leys, Christophe, et al. « Detecting Outliers: Do Not Use Standard
% Deviation around the Mean, Use Absolute Deviation around the Median ».
% Journal of Experimental Social Psychology, vol. 49, nᵒ 4, juillet 2013,
% p. 764‑66. DOI.org (Crossref),
% https://doi.org/10.1016/j.jesp.2013.03.013.
% 
% [Y1, Y2, Y3, Y4] = nHampel(X1, X2, X3, X4) applies
% Hampel outlier detection to the input signal, then
% interpolates detected outlier regions using linear interpolation. 
% 
% outputs:
% Y1 : cleaned signal with outliers replaced by smoothed values
% Y2 : binary vector marking outlier positions (1=outlier, 0=normal)
% Y3 : local medians across all time windows
% Y4 : standard deviations across all time windows
% 
% inputs:
% X1 : original data to be processed
% X2 : number of adjacent samples used for Hampel filter window
% X3 : statistical threshold multiplier (points > mean + X3 * std are outliers)
% X4 : number of iterations for Hampel filtering
% 
% See also hampel
% 
% Code created on July 11, 2023 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University.

outputSignal=inputSignal;
outlierIndexes=zeros(size(outputSignal));

if nTimesHampel > 0
    for cpt=1:nTimesHampel
        [outputSignal, outlierIndexesBuffer, XMEDIAN, XSIGMA]=hampel(outputSignal,adjacentSampleHamp,NSIGMA);
        outlierIndexes=or(outlierIndexesBuffer,outlierIndexes);
    end
end

end

