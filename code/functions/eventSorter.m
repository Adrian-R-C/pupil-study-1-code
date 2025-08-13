function [eventTypeAndOnset] = eventSorter(rawLabel)
% eventSorter extracts and organizes experimental events from raw label data
% to facilitate analysis of pupillary responses time-locked to specific events.
% This function retrieves events (fixation crosses, image presentation, pauses
% between blocks, and participant responses) and structures them with their
% temporal onsets for subsequent pupillometric analysis.
%
%   Y = eventSorter(X) parses raw event labels and
%   creates a structured array containing event information and timing data
%   for each detected experimental event.
%   
%   Event codes:
%   F = Fixation cross presentation
%   I = Image/stimulus presentation  
%   P = Pause between experimental blocks
%   R = Participant response
%
%   Code created on July 10, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

nEvents=length( rawLabel ) ;    %Number of events (stimuli onsets, fixation, etc..)

cpt=0;
for sampleIndex=1:nEvents
    
    %Type of event
    %F=Fixation
    %I=Image
    %P=Pause
    %R=Response 
    
    if length( rawLabel{ sampleIndex , 1} ) ~= 0       %i.e. There is an event
        cpt=cpt+1;

        event=rawLabel{ sampleIndex , 1};                %We look at the type of event
        eventType=event(1); %From cell to characters

        eventTypeAndOnset( cpt ).eventType=eventType; %Full name of event
        eventTypeAndOnset( cpt ).event=event; %First letter of event
        eventTypeAndOnset( cpt ).eventOnset=sampleIndex; %Index of event
    end
end
end
