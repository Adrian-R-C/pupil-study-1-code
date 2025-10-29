function [trials] = splitTrialsApp(rawSignalDiameter,eventTypeAndOnset,rawTime,participantPath,trialDuration)
% splitTrialsApp splits raw pupillary signal into individual trials.
%
%   Y = splitTrialsApp(X1, X2, X3, X4, X5) takes a raw pupillary signal
%   along with event information and metadata, and splits it into 
%   individual trials with associated timing, pupil diameter data,
%   categorization responses, baseline values, recognition outcomes,
%   stimulus characteristics, interpolation information.
%
%   output:
%     Y : structure array where each element contains:
%         - onsetBaseline : sample index of baseline start
%         - onset         : sample index of trial start
%         - offset        : sample index of trial end
%         - time          : time vector for trial duration
%         - rawSignalDiameter : raw pupil diameter samples for the trial
%         - cleanSignalDiameter : processed pupil diameter samples
%         - noise         : noise estimate (raw - clean)
%         - indexModifiedPointsByInterpolation : binary mask for interpolated points
%         - modifiedPointsPercentageByInterpolation : percentage of interpolated points
%         - indexModifiedPointsByManualInterpolation : binary mask for manually interpolated points
%         - modifiedPointsPercentageByManualInterpolation : percentage of manual interpolations
%         - indexModifiedPointsByHampel : binary mask for Hampel-filtered points
%         - modifiedPointsPercentageByHampel : percentage of Hampel outliers
%         - indexModifiedPoints : binary mask of all modified points
%         - modifiedPointsPercentage : percentage of all modified points
%         - cleanSignalBaseline : baseline period data
%         - baselineMeanValue : mean pupil diameter during baseline
%         - undetrendedBaseline : baseline value before detrending
%         - cleanSignalDiameterMinusBaseline : baseline-corrected trial signal
%         - eventFull      : full event string from input
%         - responseTime   : categorization response time (ms)
%         - categorizationResponse : participant’s categorization choice
%         - foreground     : foreground category code (e.g., 'ani', 'obj')
%         - categorizationResult : binary flag for correct categorization
%         - background     : background category code (e.g., 'nat', 'art')
%         - congruence     : trial congruence flag ('C' or 'I')
%         - imageID        : stimulus image identifier
%         - set            : stimulus set identifier
%         - snr            : signal-to-noise ratio (if available)
%         - isRecognized   : binary flag for recognition in subsequent test
%         - trustLevel     : recognition confidence score
%
%   inputs:
%     X1 : raw pupil diameter signal
%     X2 : structure array with fields:
%           - eventType : event type code ('I' for image, etc.)
%           - eventOnset : onset sample index of the event
%           - event      : full event string
%     X3 : time vector 
%     X4 : path to participant’s data directory 
%     X5 : trial duration in seconds
%
%   Code created on July 10, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University

eventArrayFull={eventTypeAndOnset.event};

eventArray=[eventTypeAndOnset.eventType];
nEvent=sum([eventArray] ~= 0);

eventArrayI=[eventTypeAndOnset.eventType] == 'I';
indexImage=find(eventArrayI);
nTrial=sum(eventArrayI);

onsetArray=[eventTypeAndOnset.eventOnset];
onsetsTrials=onsetArray(indexImage);
nSecondsAfter=trialDuration;
nSecondsBefore=0.2;
offsetsTrials=onsetsTrials+nSecondsAfter*300-1 ; % Les essais durent 6s
onsetsBaseline=onsetsTrials-nSecondsBefore*300+1 ; % Baseline = 200ms avant


txtV1=[dir('*gph').name(1:end-3) 'txt'];
txtV2=[dir('*Reconnaissance*').name];

txtV1Path=[participantPath filesep txtV1];
txtV2Path=[participantPath filesep txtV2];

responses = howIsItCategorized(txtV1Path);
responseTime = responseTimeTrials(txtV1Path);

nSignaltypes=8;

for cptTrial=1:nTrial
    cptTrial;
    trials(cptTrial).onsetBaseline=onsetsBaseline(cptTrial);
    trials(cptTrial).onset=onsetsTrials(cptTrial);
    trials(cptTrial).offset=offsetsTrials(cptTrial);
    
    trials(cptTrial).time=rawTime( onsetsTrials(cptTrial):offsetsTrials(cptTrial), : );

    trials(cptTrial).rawSignalDiameter=rawSignalDiameter( onsetsTrials(cptTrial):offsetsTrials(cptTrial), : );
%     trials(cptTrial).cleanSignalDiameter=cleanSignalDiameter( onsetsTrials(cptTrial):offsetsTrials(cptTrial), : );
    trials(cptTrial).cleanSignalDiameter=NaN(size(trials(cptTrial).rawSignalDiameter));

%     trials(cptTrial).noise=trials(cptTrial).rawSignalDiameter-trials(cptTrial).cleanSignalDiameter;
%     trials(cptTrial).indexModifiedPoints=indexModifiedPoints( onsetsTrials(cptTrial):offsetsTrials(cptTrial), : );
    trials(cptTrial).noise=NaN(size(trials(cptTrial).rawSignalDiameter));
    
    trials(cptTrial).indexModifiedPointsByInterpolation=zeros(size(trials(cptTrial).rawSignalDiameter));
    trials(cptTrial).modifiedPointsPercentageByInterpolation=sum( trials(cptTrial).indexModifiedPointsByInterpolation ) / length( trials(cptTrial).indexModifiedPointsByInterpolation );

    trials(cptTrial).indexModifiedPointsByManualInterpolation=zeros(size(trials(cptTrial).rawSignalDiameter));
    trials(cptTrial).modifiedPointsPercentageByManualInterpolation=sum( trials(cptTrial).indexModifiedPointsByManualInterpolation ) / length( trials(cptTrial).indexModifiedPointsByManualInterpolation );

    trials(cptTrial).indexModifiedPointsByHampel=zeros(size(trials(cptTrial).rawSignalDiameter));
    trials(cptTrial).modifiedPointsPercentageByHampel=sum( trials(cptTrial).indexModifiedPointsByHampel ) / length( trials(cptTrial).indexModifiedPointsByHampel );

    trials(cptTrial).indexModifiedPoints=zeros(size(trials(cptTrial).rawSignalDiameter));
    trials(cptTrial).modifiedPointsPercentage=sum( trials(cptTrial).indexModifiedPoints ) / length( trials(cptTrial).indexModifiedPoints );
        
%     trials(cptTrial).cleanSignalBaseline=cleanSignalDiameter( onsetsBaseline(cptTrial):onsetsTrials(cptTrial), : );
    trials(cptTrial).cleanSignalBaseline=NaN(length(onsetsBaseline(cptTrial):onsetsTrials(cptTrial)), nSignaltypes);
    
    trials(cptTrial).baselineMeanValue=mean(trials(cptTrial).cleanSignalBaseline) ;
    trials(cptTrial).undetrendedBaseline=trials(cptTrial).baselineMeanValue;

%     trials(cptTrial).cleanSignalDiameterMinusBaseline=trials(cptTrial).cleanSignalDiameter- trials(cptTrial).baselineMeanValue;
    trials(cptTrial).cleanSignalDiameterMinusBaseline=NaN(size(trials(cptTrial).rawSignalDiameter));

    trials(cptTrial).eventFull=eventArrayFull{indexImage(cptTrial)};
    trials(cptTrial).responseTime=responseTime(cptTrial);
    trials(cptTrial).categorizationResponse=responses(cptTrial);

    trials(cptTrial).foreground= string(trials(cptTrial).eventFull(7:9));
    trials(cptTrial).categorizationResult=isWellCategorized(trials(cptTrial).categorizationResponse,trials(cptTrial).foreground);

    trials(cptTrial).background= string(trials(cptTrial).eventFull(11:13));
    
    trials(cptTrial).congruence= trials(cptTrial).eventFull(20);
    
    trials(cptTrial).imageID= trials(cptTrial).eventFull(22:23);
    trials(cptTrial).set= trials(cptTrial).eventFull(18);

%     trials(cptTrial).snr=snr( trials(cptTrial).rawSignalDiameter , trials(cptTrial).noise );
    trials(cptTrial).snr=[];
   
    
    [trials(cptTrial).isRecognized trials(cptTrial).trustLevel]=isRecognized( trials(cptTrial).eventFull(7:23) , txtV2Path );


    %interpolation rate

end

end