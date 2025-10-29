function responseTime = responseTimeTrials(txtPath)
% responseTimeTrials extracts response times for each trial in the categorization task.
%
%   Y = responseTimeTrials(X) reads the .txt behavioral file produced during the 
%   categorization task and returns the response time (RT) for each trial.
%
%   output:
%     Y : Vector of response times (in milliseconds) for all trials in the task.
%
%   input:
%     X : Path to the .txt file obtained during the categorization task.
%
%   Code created on July 13, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University


txt=importdata(txtPath);

cpt=1;
for line=1:length(txt)
    
    %buffer=split(string(txt{line}), ';')
    buffer=txt{line};
    
    if strfind(buffer, 'N/A')
        responseTime(cpt)=NaN;
        cpt=cpt+1;
        
    elseif  strfind(txt{line}, 'ms') %Y a un temps de réponse écrit
        responseTime_=split(string(buffer), ';');
        responseTime_=convertStringsToChars(responseTime_(3));
        responseTime(cpt)=str2num(responseTime_(1:end-3));
        cpt=cpt+1;
    end
    
end
        
end