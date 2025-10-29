function [pValueBuffer] = uncorrectedSignificanceBetweenTwoMatrices(matrix1, matrix2, dependentOrIndependent, tail)
% uncorrectedSignificanceBetweenTwoMatrices computes uncorrected p-values for pointwise comparisons between two datasets.
%
%   Y = uncorrectedSignificanceBetweenTwoMatrices(X1, X2, X3, X4) performs 
%   either independent or paired statistical tests at each 
%   time point between two input matrices and returns the resulting 
%   uncorrected p-values.
%
%   output:
%     Y : numeric vector (nTimepoints × 1) containing the uncorrected p-values
%         for each pointwise comparison between X1 and X2.
%
%   inputs:
%     X1 : numeric matrix (nTimepoints × nSubjects1) for condition/group 1
%     X2 : numeric matrix (nTimepoints × nSubjects2) for condition/group 2
%     X3 : string, either "independent" (unpaired t-test) or "dependent" (paired t-test)
%     X4 : string specifying the test tail:
%             'both'   - two-tailed test
%             'right'  - right-tailed test
%             'left'   - left-tailed test
%
%  PS: After reviewing all my functions, I realized that the ttest and
%  ttest2 functions basically work on matrices, so this function may be
%  unnecessary, or even slower... My apologies...  
%
%   Code created on August 24, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University


[numberOfPoints numberOfSubject1]=size(matrix1);
[numberOfPoints numberOfSubject2]=size(matrix2);


if dependentOrIndependent=="independent" %unpaired

    for cptPoint=1:numberOfPoints

        for cptSubject=1:numberOfSubject1
            distributionAtChosenPoint1(cptSubject)=matrix1(cptPoint , cptSubject);
        end
                        
        for cptSubject=1:numberOfSubject2
            distributionAtChosenPoint2(cptSubject)=matrix2(cptPoint , cptSubject);
        end

        %We remove nan columns
        distributionAtChosenPoint1=distributionAtChosenPoint1(~isnan(distributionAtChosenPoint1));
        distributionAtChosenPoint2=distributionAtChosenPoint2(~isnan(distributionAtChosenPoint2));

        significanceBuffer(cptPoint)=0;
        pValueBuffer(cptPoint)=NaN;
                        
        if ~ (isempty(distributionAtChosenPoint1) || isempty(distributionAtChosenPoint2))
%             [p h]=ranksum(distributionAtChosenPoint1,distributionAtChosenPoint2, 'tail', tail);
            [h p]=ttest2(distributionAtChosenPoint1,distributionAtChosenPoint2, 'tail', tail);
            significanceBuffer(cptPoint)=h;
            pValueBuffer(cptPoint)=p;
        end
                        
    end
    

elseif dependentOrIndependent=="dependent" %paired
 
    numberOfSubject=numberOfSubject1;
    
    distributionOfCondition1AtChosenPoint=[];
    distributionOfCondition2AtChosenPoint=[]; 
    
      for cptPoint=1:numberOfPoints
       


        for cptSubject=1:numberOfSubject1
            distributionOfCondition1AtChosenPoint(cptSubject)=matrix1(cptPoint , cptSubject);
            distributionOfCondition2AtChosenPoint(cptSubject)=matrix2(cptPoint , cptSubject);
        end
                
        %we compute the significance
        significanceBuffer(cptPoint)=0;
        pValueBuffer(cptPoint)=NaN;
                          
        
        if ~ (isempty(distributionOfCondition1AtChosenPoint) || isempty(distributionOfCondition2AtChosenPoint))

            %We remove nan columns
            whoIsNanCondition1=isnan(distributionOfCondition1AtChosenPoint);
            whoIsNanCondition2=isnan(distributionOfCondition2AtChosenPoint);

            whoIsNan=or(whoIsNanCondition1, whoIsNanCondition2);
            distributionOfCondition1AtChosenPoint=distributionOfCondition1AtChosenPoint(~whoIsNan);
            distributionOfCondition2AtChosenPoint=distributionOfCondition2AtChosenPoint(~whoIsNan);
                                       
            if isequal( size(distributionOfCondition1AtChosenPoint) , size(distributionOfCondition2AtChosenPoint) )
%                 [p h]=signrank(distributionOfCondition1AtChosenPoint,distributionOfCondition2AtChosenPoint, 'tail', tail );
                [h p]=ttest(distributionOfCondition1AtChosenPoint,distributionOfCondition2AtChosenPoint, 'tail', tail );
            else
                error(sprintf("Sizes not matching"))
            end

            significanceBuffer(cptPoint)=h;
            pValueBuffer(cptPoint)=p;

        end
    end
    
else
    error("Only valid arguments : ""dependent"" or ""indepedent"" ")
end
    
    
end

