function [folderParticipant] = subjectCodeAndFolder( folderParticipantPath )
% subjectCodeAndFolder gets participants study code to subsequently
% identify them as 'young' or 'old' subject
%
%   Y = subjectCodeAndFolder(X) 
%
%   output:
%     Y : participant study code (XXjXX format for younger participants, XXaXX
%     for older participants)
%
%   input:
%     X : participant path
%
%   Code created on July 10, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University

folderParticipantBuffer = dir(folderParticipantPath);

nSubject=length(folderParticipantBuffer)-2;

for noSubject=1:nSubject
    if  folderParticipantBuffer(noSubject+2).isdir == 1 
        folderParticipant(noSubject).subjectCode = folderParticipantBuffer(noSubject+2).name; %+2 car on prend en compte . et ..
        folderParticipant(noSubject).subjectCodeFolder = folderParticipantBuffer(noSubject+2).folder;
    end
end

end

