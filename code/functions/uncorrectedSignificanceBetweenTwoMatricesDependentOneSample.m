function [pValueBuffer] = uncorrectedSignificanceBetweenTwoMatricesDependentOneSample(matrix1)
% uncorrectedSignificanceBetweenTwoMatricesDependentOneSample computes uncorrected p-values for one-sample pointwise comparisons against zero.
%
%   Y = uncorrectedSignificanceBetweenTwoMatricesDependentOneSample(X) performs
%   a one-sample statistical test at each time point of the input matrix
%   and returns the resulting uncorrected p-values.
%
%   output:
%     Y : numeric vector (nTimepoints × 1) containing the uncorrected p-values
%         for each pointwise comparison of X against zero.
%
%   input:
%     X : numeric matrix (nTimepoints × nSubjects) representing measurements
%          for a single condition or group.
%
%  PS: After reviewing all my functions, I realized that the ttest function
%  basically work on matrices, so this function may be 
%  unnecessary, or even slower... My apologies...  
%
%   Code created on August 8, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University


[numberOfPoints numberOfSubject]=size(matrix1);

distributionOfCondition1AtChosenPoint=[];

    for cptPoint=1:numberOfPoints

        for cptSubject=1:numberOfSubject
            distributionOfCondition1AtChosenPoint(cptSubject)=matrix1(cptPoint , cptSubject);
        end

        %we compute the significance
        significanceBuffer(cptPoint)=0;
        pValueBuffer(cptPoint)=NaN;

        %We remove nan columns
            whoIsNan=isnan(distributionOfCondition1AtChosenPoint);

            distributionOfCondition1AtChosenPoint=distributionOfCondition1AtChosenPoint(~whoIsNan);
            
        if ~ isempty(distributionOfCondition1AtChosenPoint) 
           

%             [p h]=signrank(distributionOfCondition1AtChosenPoint);
            [h p]=ttest(distributionOfCondition1AtChosenPoint);


            significanceBuffer(cptPoint)=h;
            pValueBuffer(cptPoint)=p;
        end
    

    
    end

end

