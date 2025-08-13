function [outputSignal, indexModifiedPoints] = noProcessing(inputSignal,structContainingFilteringParameters)
% noProcessing does... well nothing actually, but it's necessary, I swear !
%
%   Code created on March 26, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

outputSignal=inputSignal;
indexModifiedPoints.hampel=zeros(size(inputSignal));
indexModifiedPoints.interpolation=zeros(size(inputSignal));
indexModifiedPoints.total=zeros(size(inputSignal));

end

