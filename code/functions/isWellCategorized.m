function result = isWellCategorized(categorization,foreground)
% isWellCategorized determines if a given stimulus is well categorized
%
% Y = isRecognized(X1, X2) reads X1 and see if it corresponds to the
% foreground X2 of the stimulus 
% 
% output:
% Y : categorization status ('G' = well categorized, 'B' = bad categorized, 'N/A' = not available)
% 
% inputs:
% Y1 : categorization reponse
% Y2 : stimulus foreground object type
%
% Code created on September 4, 2023 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University

    if categorization==foreground
        result="G";
    else
        result="B";
    end

    if categorization=="N/A" 
        result="N/A";
    end

end

