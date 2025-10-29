function [timeBlock, rawSignalBlock, rawSignalBlockIndexes, rawTime, rawSignalDiameter] = fourBlocks(rawTime, rawSignalDiameter, eventTypeAndOnset, app)
% fourBlocks segments continuous eye-tracking recordings into experimental
% blocks, removing unwanted data during participant
% rest periods between blocks.
%
% [Y1, Y2, Y3, Y4, Y5] = fourBlocks(X1, X2, X3, X4) partitions raw
% eye-tracking data into four experimental blocks based on event markers,
% handling timing irregularities and signal discontinuities that may occur
% when breaks are initiated prematurely.
% 
% outputs:
% Y1 : cell array containing time vectors for each of the four blocks
% Y2 : cell array containing pupil signal data for each block
% Y3 : cell array with event onsets and offsets for each block
% Y4 : extended time vector with additional points to handle early breaks
% Y5 : extended signal vector with replicated end points for timing correction
% 
% inputs:
% X1 : original time vector from eye-tracking recording
% X2 : raw pupil diameter signal
% X3 : event triggers and experimental labels
% X4 : MATLAB application for data transfer purpose
%
% 1 trial = 1 Fixation-Image-Response cycle
% 4 blocks of 18 F-I-R cycles divided by 3 pauses
% The signals : OGx, OGy, ODx, ODy MeanX MeanY SurfG SurfD 
% (OG = left eye, OD = right eye)
%
% Code created on July 10, 2023 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University

indexImage=find([eventTypeAndOnset.eventType] == 'I');
indexResponse=find([eventTypeAndOnset.eventType] == 'R');
indexFixation=find([eventTypeAndOnset.eventType] == 'F');
indexPause=find([eventTypeAndOnset.eventType] == 'P');

% Because the trials 18, 36, 48 and 72 are not well designed (pause coming
% too soon), we need to extend the pause onset by a few seconds (10s ?)
secondToAddToLastBlock=10;
pauseExtentedOffsetForBlock4=round(app.samplingFrequency*secondToAddToLastBlock);
% pauseExtentedOffsetForBlock4=0;
additionalVector=rawTime(end):1/app.samplingFrequency:(rawTime(end)+secondToAddToLastBlock);
rawTime=[rawTime(1:end-1);additionalVector' ];
pointsToAdd=length(rawTime)-length(rawSignalDiameter);
rawSignalDiameter=[rawSignalDiameter; repmat(rawSignalDiameter(end,:), pointsToAdd,1)];
% %
indexFixationBlocksBegin=[0 1 2 3]*(72/4)+1;

onsetBlock1=eventTypeAndOnset(indexFixation(indexFixationBlocksBegin(1))).eventOnset ;
offsetBlock1=eventTypeAndOnset(indexPause(1)).eventOnset;

onsetBlock2=eventTypeAndOnset(indexFixation(indexFixationBlocksBegin(2))).eventOnset ;
offsetBlock2=eventTypeAndOnset(indexPause(2)).eventOnset;

onsetBlock3=eventTypeAndOnset(indexFixation(indexFixationBlocksBegin(3))).eventOnset ;
offsetBlock3=eventTypeAndOnset(indexPause(3)).eventOnset;

onsetBlock4=eventTypeAndOnset(indexFixation(indexFixationBlocksBegin(4))).eventOnset ;
offsetBlock4=eventTypeAndOnset(indexPause(4)).eventOnset + pauseExtentedOffsetForBlock4;

rawSignalBlock1=rawSignalDiameter(onsetBlock1:offsetBlock1, :);
rawSignalBlock2=rawSignalDiameter(onsetBlock2:offsetBlock2, :);
rawSignalBlock3=rawSignalDiameter(onsetBlock3:offsetBlock3, :);
rawSignalBlock4=rawSignalDiameter(onsetBlock4:offsetBlock4, :);

timeBlock1=rawTime(onsetBlock1:offsetBlock1);
timeBlock2=rawTime(onsetBlock2:offsetBlock2);
timeBlock3=rawTime(onsetBlock3:offsetBlock3);
timeBlock4=rawTime(onsetBlock4:offsetBlock4);

rawSignalBlock{1} = [rawSignalBlock1 ];
rawSignalBlock{2} = [rawSignalBlock2 ];
rawSignalBlock{3} = [rawSignalBlock3 ];
rawSignalBlock{4} = [rawSignalBlock4 ];

timeBlock{1} = timeBlock1 ;
timeBlock{2} = timeBlock2 ;
timeBlock{3} = timeBlock3 ;
timeBlock{4} = timeBlock4 ;

rawSignalBlockIndexes = [ onsetBlock1 ;
                          offsetBlock1 ;
                          onsetBlock2 ;
                          offsetBlock2
                          onsetBlock3 ;
                          offsetBlock3
                          onsetBlock4 ;
                          offsetBlock4 ;];
end

