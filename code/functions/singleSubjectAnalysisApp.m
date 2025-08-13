function [trials timeBlock, rawSignalDiameterBlock, rawSignalBlockIndexes, rawTime, rawSignalGaze, rawSignalDiameter, rawLabel] = singleSubjectAnalysisApp(participantPath, app)
% singleSubjectAnalysisApp splits raw data into 72 trials.
%
%   [Y1, Y2, Y3, Y4, Y5, Y6, Y7] = responseTimeTrials(X1, X2) 
%
%   outputs:
%     Y1 : structure containing all trials
%     Y2 : time vector split into four blocks
%     Y3 : data vector split into four blocks
%     Y4 : event onsets and offsets for each block
%     Y5 : time vector
%     Y6 : gaze raw data 
%     Y7 : pupil raw data
%     Y8 : event onsets and offset
%
%   inputs:
%     X1 : path to the CSV file obtained during the categorization
%     task containing the data
%     X2 : MATLAB app
%
%   Code created on July 17, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University
cd(participantPath)
csvPath=[participantPath filesep dir('Export*').name];

participantPath

% On met le path du csv dans la fonction csvSpliter et on récupère nos 10
% vecteurs (le vecteur 10 de Labels est de type cell)

[rawTime, rawSignalGaze, rawSignalDiameter, rawLabel] = csvSpliter(csvPath);

rawTime=rawTime/1000; %ms en s

% On regarde l'ensemble des labels (fixation, stimulus, réponse, pause) 

eventTypeAndOnset=eventSorter(rawLabel);


%On sépare les signaux dans les 4 blocs (pour le temps, et LDx, LDy, RDx,
%RDy, Moyenne X, Moyenne Y, surface OG, surface OD on les sépare en 4 blocs);
[timeBlock, rawSignalDiameterBlock, rawSignalBlockIndexes, rawTime, rawSignalDiameter] = fourBlocks(rawTime, rawSignalDiameter, eventTypeAndOnset, app) ;

% 
% % Filtrage des signaux d'intéret (LDx,LDy,RDx,RDy, Surface OG, Surface OD)
% for block=1:4
%     
%     bufferRawBlock=rawSignalDiameterBlock{block};
%     [nPoints nSignalTypes]=size(bufferRawBlock);
%     for signalOfInterest=1:nSignalTypes
%         
%         bufferRawSignalType=bufferRawBlock(:,signalOfInterest);
%         
%         [cleanSignalDiameterBlock{block}(:,signalOfInterest) indexModifiedPointsBlock{block}(:,signalOfInterest) ]=processingHampel(bufferRawSignalType);
%     
%     end
%     
% end

% Concatenation des 4 blocs et reconstruction des signaux entiers 
% cleanSignalDiameter=rawSignalDiameter;
% indexModifiedPoints=zeros(size(rawSignalDiameter));
% 
% cleanSignalDiameter( rawSignalBlockIndexes(1):rawSignalBlockIndexes(2) , : ) = cleanSignalDiameterBlock{1};
% cleanSignalDiameter( rawSignalBlockIndexes(3):rawSignalBlockIndexes(4) , : ) = cleanSignalDiameterBlock{2};
% cleanSignalDiameter( rawSignalBlockIndexes(5):rawSignalBlockIndexes(6) , : ) = cleanSignalDiameterBlock{3};
% cleanSignalDiameter( rawSignalBlockIndexes(7):rawSignalBlockIndexes(8) , : ) = cleanSignalDiameterBlock{4};

% indexModifiedPoints( rawSignalBlockIndexes(1):rawSignalBlockIndexes(2) , : ) = indexModifiedPointsBlock{1};
% indexModifiedPoints( rawSignalBlockIndexes(3):rawSignalBlockIndexes(4) , : ) = indexModifiedPointsBlock{2};
% indexModifiedPoints( rawSignalBlockIndexes(5):rawSignalBlockIndexes(6) , : ) = indexModifiedPointsBlock{3};
% indexModifiedPoints( rawSignalBlockIndexes(7):rawSignalBlockIndexes(8) , : ) = indexModifiedPointsBlock{4};

% On coupe en 72 essais ;
trialDuration=app.trialDuration;
trials = splitTrialsApp(rawSignalDiameter,eventTypeAndOnset,rawTime,participantPath, trialDuration);


end



