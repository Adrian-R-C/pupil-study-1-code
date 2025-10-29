function [outputSignal, outlierIndexes] = derivativeHampelInterp(inputSignal, adjacentSampleHamp, NSIGMA, nTimesHampel, adjacentSampleInterp)
% derivativeHampelInterp detects and interpolates outliers in time series
% data using derivative-based Hampel filtering. This function identifies
% rapid signal changes that exceed statistical thresholds, marks surrounding
% points as outliers, and replaces them with interpolated values. The
% function computes the absolute derivative of the input signal, applies 
% n-times Hampel filtering to detect outliers, marks surrounding regions for
% interpolation, and uses linear interpolation to replace outlier values.
% Edge cases are handled by adjusting interpolation windows at signal boundaries.
%
% See Leys, Christophe, et al. « Detecting Outliers: Do Not Use Standard
% Deviation around the Mean, Use Absolute Deviation around the Median ».
% Journal of Experimental Social Psychology, vol. 49, nᵒ 4, juillet 2013,
% p. 764‑66. DOI.org (Crossref),
% https://doi.org/10.1016/j.jesp.2013.03.013.

%   [Y1, Y2] = derivativeHampelInterp(X1, X2, X3, X4, X5) applies
%   Hampel outlier detection to the derivative of the input signal, then
%   interpolates detected outlier regions using linear interpolation. 
%   
%   outputs:
%   Y1 : cleaned signal with outliers replaced by interpolated values
%   Y2 : binary vector marking outlier positions (1=outlier, 0=normal)
%
%   inputs:
%   X1 : original data to be processed
%   X2 : number of adjacent samples used for Hampel filter window
%   X3 : statistical threshold multiplier (points > mean + X3 * std are outliers)
%   X4 : number of iterations for Hampel filtering
%   X5 : number of adjacent samples around each outlier to interpolate
%
%   Code created on July 11, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

outputSignal=inputSignal;
outlierIndexes=zeros(size(outputSignal));

deriv=abs(diff(inputSignal));
interpWindow=2*adjacentSampleInterp+1;
xAxis=1:length(inputSignal);

%Hampel outliers detector
[derivHampel outlierDerivIndexesBinary XMEDIAN, XISGMA]=nHampel(deriv,adjacentSampleHamp, NSIGMA, nTimesHampel);

outlierDerivIndexes=find(outlierDerivIndexesBinary);

for cpt=1:length(outlierDerivIndexes)
        
        if (outlierDerivIndexes(cpt) - adjacentSampleInterp ) <= 1 %i.e. Left edge
            outputSignal(2:outlierDerivIndexes(cpt)+adjacentSampleInterp)=NaN;
            outlierIndexes(2:outlierDerivIndexes(cpt)+adjacentSampleInterp)=1;
            
        elseif ( outlierDerivIndexes(cpt) + adjacentSampleInterp) >= length(inputSignal)  %i.e. Right edge
            outputSignal(outlierDerivIndexes(cpt)-adjacentSampleInterp:end-1)=NaN;
            outlierIndexes(outlierDerivIndexes(cpt)-adjacentSampleInterp:end-1)=1;
      
        else
            outputSignal(outlierDerivIndexes(cpt)-adjacentSampleInterp:outlierDerivIndexes(cpt)+adjacentSampleInterp)=NaN;
            outlierIndexes(outlierDerivIndexes(cpt)-adjacentSampleInterp:outlierDerivIndexes(cpt)+adjacentSampleInterp)=1;

        end
end


indexToRemove = isnan(outputSignal);

outputSignal(indexToRemove) = [];
sum(isnan(outputSignal));

xAxisToRemove=xAxis;
xAxisToRemove(indexToRemove) = [];

outputSignal=interp1(xAxisToRemove, outputSignal, xAxis, 'linear');
sum(isnan(outputSignal));


end

