function [recognizedOrNot trustLevel]= isRecognized( eventFull , txtV2Path )
% isRecognized determines recognition status and confidence level for a given
% stimulus based on post-encoding memory test performance data.
%
% [Y1, Y2] = isRecognized(X1, X2) reads
% recognition test results to identify whether a specific stimulus was
% subsequently recognized and retrieves the associated confidence rating.
% 
% outputs:
% Y1 : recognition status ('R' = recognized, 'N' = not recognized)
% Y2 : confidence level/trust rating associated with the recognition
% response (R, K, or G)
% 
% inputs:
% X1 : stimulus identifier string to be matched against recognition data
% X2 : file path to the recognition test results file
%
% Code created on July 13, 2023 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University

txt=importdata(txtV2Path);

inputStim=eventFull;

for line=5:length(txt)
    
    buffer=txt{line};
    buffer2=split(string(buffer), ';');
        
    if length(buffer2)==7
        stim=convertStringsToChars(buffer2(1));
   
    
        newStim(1:3)=stim(6:8);
        newStim(4)='_';
        newStim(5:7)=stim(10:12);
        newStim(8)='_';
        newStim(9:12)=stim(1:4);
        newStim(13)='_';
        newStim(14:17)=stim(14:17);
        
        string(newStim);
        
        if string(newStim)==string(inputStim) && buffer2(2) =="Oui"
            recognizedOrNot='R';
            trustLevel=buffer2(5);
        elseif string(newStim)==string(inputStim) && buffer2(2)=="Non"
            recognizedOrNot='N';
            trustLevel=buffer2(5);
        end
    
    end

end

end

