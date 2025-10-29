function [rawTime, rawSignalGaze, rawSignalDiameter, rawLabel] = csvSpliter(csvPath)
% csvSpliter reads CSV data and separates it into structured output vectors
% containing timepoints, gaze direction, pupil diameter, and event labels
% from eye-tracking recordings sampled at 300Hz.
%
%   [Y1, Y2, Y3, Y4] = csvSpliter(X) loads and parses CSV file data, extracting temporal information, gaze
%   coordinates, pupil measurements, and experimental event markers.
%   
%   outputs:
%   Y1 : time vector (sampled at 300Hz)
%   Y2 : gaze direction data matrix containing left eye X/Y and
%   right eye X/Y coordinates [OGx, OGy, ODx, ODy] 
%   Y3 : pupil diameter matrix with left and right eye measurements [OGx,
%   OGy, ODx, ODy]   
%   Y4 : event labels vector (F=Fixation, I=Image, P=Pause, R=Response)
%
%   input:
%   X : file path to the CSV data file
%
%   Code created on July 10, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.


opts = detectImportOptions(csvPath);
opts.VariableTypes{10}='string';

preview(csvPath,opts);

csvRaw=readtable(csvPath, opts);

rawTime=table2array(csvRaw(:,1));

rawLXdeg=table2array(csvRaw(:,2));
rawLYdeg=table2array(csvRaw(:,3));
rawRXdeg=table2array(csvRaw(:,4));
rawRYdeg=table2array(csvRaw(:,5));
rawSignalGaze=[rawLXdeg  rawLYdeg  rawRXdeg  rawRYdeg];

rawLDXpix=table2array(csvRaw(:,6));
rawLDYpix=table2array(csvRaw(:,7));
rawRDXpix=table2array(csvRaw(:,8));
rawRDYpix=table2array(csvRaw(:,9));

%
rawSurfaceLeft=pi*(rawLDXpix/2).*(rawLDYpix/2); %surface ellipse = pi * r1 * r2 = (pi/4) * d1 * d2
rawSurfaceRight=pi*(rawRDXpix/2).*(rawRDYpix/2);
% rawSignalDiameter=[rawLDXpix  rawLDYpix  rawRDXpix  rawRDYpix rawSurfaceLeft rawSurfaceRight];


rawDxMean=(rawRDXpix+rawLDXpix)/2; 
rawDyMean=(rawRDYpix+rawLDYpix)/2;
rawSignalDiameter=[rawLDXpix  rawLDYpix  rawRDXpix  rawRDYpix rawDxMean rawDyMean rawSurfaceLeft rawSurfaceRight];


rawLabel=table2array(csvRaw(:,10));
rawLabel(ismissing(rawLabel)) = "";

end


