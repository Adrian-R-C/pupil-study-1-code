function [outputSignal, indexModifiedPoints] = processingMethod(inputSignal,structContainingFilteringParameters, diameterOrSurface, detrending)
% processingMethod processes raw pupillary data by artefact correction, blink detection, filtering, and detrending.
%
%   [Y1, Y2] = processingMethod(X1, X2, X3, X4)
%   applies a multi-step processing pipeline to clean raw pupil diameter or pupil surface data.
%   The method detects and corrects artefacts such as blinks, interpolates missing segments,
%   low-pass filters and detrends the signal, and removes residual outliers using the Hampel filter.
%
%   Steps:
%     0. (optional) initial Hampel filtering to detect missing data.
%     1. blink detection via:
%          - thresholding the first derivative of the signal
%          - thresholding the signal itself
%        blink segments are interpolated linearly.
%     2. low-pass filtering at 4 Hz, with optional detrending.
%     3. final Hampel filtering to remove remaining outliers.
%
%   outputs:
%     Y1 : cleaned signal after interpolation, filtering, detrending, and outlier removal
%     Y2 : struct with logical indices marking modified points
%         .hampel             : Points removed by Hampel filtering
%         .interpolation      : Points replaced by interpolation (blink correction)
%         .total              : Points modified by either Hampel or interpolation
%
%   inputs:
%     X1 : raw pupillary signal (vector)
%     X2 : struct with Hampel filtering parameters:
%         .adjacentSampleHamp1     : Number of adjacent samples for Hampel filter window
%         .nSigma1                 : Threshold multiplier
%         .nTimesHampel1           : Number of Hampel filtering iterations
%     X3 : string, 'diameter' or 'surface' (sets blink detection thresholds)
%     X4 : string, 'detrend' to remove linear trend, anything else to skip detrending
%
%   Code created on July 11, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University

outputSignal=inputSignal;
indexModifiedPoints.hampel=zeros(size(inputSignal));
indexModifiedPoints.interpolation=zeros(size(inputSignal));
indexModifiedPoints.total=zeros(size(inputSignal));


% Step 0 (optional): missing data detection with hampel filtering
outputSignal=hampel(outputSignal, 5);



% Step 1: blink detection by thresholding the 1st order derivative
if diameterOrSurface == "diameter"
    threshold=40;
else
    threshold=60000;    
end

firstOrderDerivative = diff(outputSignal);
absFirstOrderDerivativeBuffer = abs(firstOrderDerivative);
absFirstOrderDerivative=[absFirstOrderDerivativeBuffer(1) ; absFirstOrderDerivativeBuffer];

indexModifiedPoints.interpolation=absFirstOrderDerivative > threshold;
% Step 1bis: blink detection by thresholding the signal itself
if diameterOrSurface == "diameter"
    threshold2=100;
else
    threshold2=25000;    
end

indexModifiedPoints.interpolation=or(indexModifiedPoints.interpolation, abs(detrend(outputSignal)) > threshold2);


indexModifiedPoints.interpolation = addFewMillisecondsBeforeAndAfterInterpol(indexModifiedPoints.interpolation, 0.050, 300);

if sum(indexModifiedPoints.interpolation)
                    
        bufferX=1:length(inputSignal);
        bufferXtruncated=bufferX;
        bufferYtruncated=inputSignal;

        bufferXtruncated(logical(indexModifiedPoints.interpolation))=NaN;                    
        bufferYtruncated(logical(indexModifiedPoints.interpolation))=NaN;

        
        bufferXtruncated(1)=bufferX(1);
        bufferXtruncated(end)=bufferX(end);
        bufferYtruncated(1)=inputSignal(1);
        bufferYtruncated(end)=inputSignal(end);
        
        bufferXtruncated(isnan(bufferXtruncated))=[];                    
        bufferYtruncated(isnan(bufferYtruncated))=[];                    

        
        outputSignal=interp1(bufferXtruncated,bufferYtruncated,bufferX, 'linear');
%         outputSignal=bufferYtruncated;
end

% Step 2: denoising by lowpassing the blink-corrected signal and then detrending
% it if wanted
outputSignal=lowpass(outputSignal,4,300); 

if detrending=="detrend"
    outputSignal=detrend(outputSignal);
end

%%% Step 3: removing the last outliers with hampel procedure
        % adjacentSampleHamp1=30; %37*2 +1 = 75pts = 1/4 de s
        % nSigma1=2;
        % nTimesHampel1=1;
%We extract the filtering parameters
adjacentSampleHamp1=structContainingFilteringParameters.adjacentSampleHamp1;
nSigma1=structContainingFilteringParameters.nSigma1;
nTimesHampel1=structContainingFilteringParameters.nTimesHampel1;

[outputSignal, indexModifiedPoints.hampel, XMEDIAN, XISGMA] = nHampel( outputSignal , adjacentSampleHamp1 , nSigma1 , nTimesHampel1);
 if size(indexModifiedPoints.hampel, 2) > 1
        indexModifiedPoints.hampel=indexModifiedPoints.hampel';
 end

% indexModifiedPoints.hampel=indexModifiedPoints.hampel';
indexModifiedPoints.total=or(indexModifiedPoints.hampel,indexModifiedPoints.interpolation);
% 
% % Step 4: we eventually replace the interpolated points by NaN
% outputSignal(logical(indexModifiedPoints.interpolation))=NaN;

%END
if size(outputSignal, 2) > 1
    outputSignal=outputSignal';
end
end

