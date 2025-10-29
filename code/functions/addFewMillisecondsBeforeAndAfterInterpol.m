function [indexModifiedPointsByInterpolation] = addFewMillisecondsBeforeAndAfterInterpol(indexModifiedPointsByInterpolation, duration, samplingFreq)
% addFewMillisecondsBeforeAndAfterInterpol extends blink detection periods
% by marking additional time points as TRUE before and after each detected
% blink, enabling removal of a specified time window surrounding blink
% events.
%
%   outputs:
%   Y = addFewMillisecondsBeforeAndAfterInterpol(X1, X2, X3) dilates the
%   TRUE regions in binary vector X1 by duration X2 for a signal with
%   sampling rate X3.
%
%   inputs:
%   X1 : binary vector with 1 indicating detected blink locations
%   X2 : extension duration (in seconds)
%   X3 : signal sampling frequency (Hz)
%   Y  : extended blink detection
%   Example:
%       result = addFewMillisecondsBeforeAndAfterInterpol([0, 0, 1, 0, 0], 1, 1)
%       % Returns: [0, 1, 1, 1, 0]
%
%   Code created on February 05, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

   howManyOnesToAdd = round(duration * samplingFreq);

    if any(indexModifiedPointsByInterpolation ~= 0 & indexModifiedPointsByInterpolation ~= 1)
        error('Needs to contain 0 and 1 only FFS');
    end

    bufferVector = indexModifiedPointsByInterpolation;

    indexes = find(bufferVector == 1);

    for cpt = 1:length(indexes)
        deb = max(1, indexes(cpt) - howManyOnesToAdd);
        fin = min(length(indexModifiedPointsByInterpolation), indexes(cpt) + howManyOnesToAdd);
        bufferVector(deb:fin) = 1;
    end

    indexModifiedPointsByInterpolation=logical(bufferVector);
end

