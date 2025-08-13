function categorization = howIsItCategorized(txtPath)
% howIsItCategorized extracts participant categorization responses from
% behavioral data files, determining whether each stimulus was classified
% as "Animal" or "Object" during the experimental task.
%
% Y = howIsItCategorized(X) reads and processes behavioral response data
% to create a vector of participant categorization decisions. 
% 
% outputs:
% Y : response vector containing 'ani' (animal) or 'obj' (object, i.e. furniture) classifications
% 
% inputs:
% X : file path to the behavioral data file containing participant responses
%
% Code created on September 4, 2023 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University


categorization=repmat("N/A", 1, 72);

if exist(txtPath)
    txt=importdata(txtPath);
    cpt=0;
    for line=5:length(txt)

        buffer=txt{line};
        buffer2=split(string(buffer), ';');
        
        if length(buffer2) >= 3
            cpt=cpt+1;
            categorization_=buffer2(2);

            if categorization_=="Mobilier"
                categorization(cpt)='obj';
            elseif categorization_=="Animal"
                categorization(cpt)='ani';
            end
            
        end
    end
end

categorization=categorization';

end

