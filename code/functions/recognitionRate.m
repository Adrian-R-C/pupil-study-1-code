function recognitionRates=recognitionRate(learnedSetV1, txtV2Path)
% recognitionRate computes recognition rates for multiple experimental conditions.
%
%   Y = recognitionRate(X1, X2) calculates the proportion of correct responses 
%   during a recognition test, both overall and separately for several stimulus 
%   conditions.
%
%   output:
%     Y : 1x7 vector containing:
%         Y(1) - Overall recognition rate (all stimuli)
%         Y(2) - Recognition rate for animal foregrounds
%         Y(3) - Recognition rate for object foregrounds
%         Y(4) - Recognition rate for natural backgrounds
%         Y(5) - Recognition rate for artificial backgrounds
%         Y(6) - Recognition rate for congruent scenes
%         Y(7) - Recognition rate for incongruent scenes
%
%   inputs:
%     X1 : Learned set identifier (e.g., 'set1' or 'set2')
%     X2 : Path to the .txt file obtained during recognition testing
%     containing behavioral data
%
%   Code created on July 11, 2023 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University


txt=importdata(txtV2Path);

goodAnswer=0;
goodAnswerAnimal=0;
goodAnswerObject=0;
goodAnswerNatural=0;
goodAnswerArtificial=0;
goodAnswerCongruent=0;
goodAnswerIncongruent=0;

for line=5:length(txt)
    
    buffer=txt{line};
    buffer2=split(string(buffer), ';');
        
    if length(buffer2)==7
        stim=convertStringsToChars(buffer2(1));
   
        presentedSetV2=stim(1:4);
        stimulusConditionForeground=string(stim(6:8));
        stimulusConditionBackground=string(stim(10:12));
        stimulusConditionCongruence=string(stim(14));
        
        if (string(learnedSetV1)==string(presentedSetV2) && buffer2(2) =="Oui") || (string(learnedSetV1)~=string(presentedSetV2) && buffer2(2) =="Non")
            goodAnswer=goodAnswer+1;
            if stimulusConditionForeground=="ani"
                goodAnswerAnimal=goodAnswerAnimal+1;
            else
                goodAnswerObject=goodAnswerObject+1;
            end
            
            if stimulusConditionBackground=="nat"
                goodAnswerNatural=goodAnswerNatural+1;
            else
                goodAnswerArtificial=goodAnswerArtificial+1;
            end
            
            if stimulusConditionCongruence=="C"
                goodAnswerCongruent=goodAnswerCongruent+1;
            else
                goodAnswerIncongruent=goodAnswerIncongruent+1;
            end
            
        end
    
    end

end


recognitionRates(1)=goodAnswer/144;
recognitionRates(2)=goodAnswerAnimal/72;
recognitionRates(3)=goodAnswerObject/72;
recognitionRates(4)=goodAnswerNatural/72;
recognitionRates(5)=goodAnswerArtificial/72;
recognitionRates(6)=goodAnswerCongruent/72;
recognitionRates(7)=goodAnswerIncongruent/72;


end

