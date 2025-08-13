function [permutedMatrix1 permutedMatrix2]= subjectMatrixPermutation(matrix1, matrix2, dependentOrIndependent)
% subjectMatrixPermutation permutes data for subsequent cluster based permutation testing.
%
%   [Y1, Y2] = subjectMatrixPermutation(X1, X2, X3) performs a permutation
%   of subjects in two matrices according to the specified
%   permutation type, for use in permutation tests.
%
%   outputs:
%     Y1 : permuted version of X1 after column reallocation
%     Y2 : permuted version of X2 after column reallocation
%
%   inputs:
%     X1 : numeric matrix (nTimepoints × nSubjects1) containing data from group 1
%     X2 : numeric matrix (nTimepoints × nSubjects2) containing data from group 2
%     X3 : permutation type, specified as:
%           - "independent" : unpaired permutation. Columns from X1 and X2 
%                             are combined, shuffled, and split back into Y1 and Y2.
%           - "dependent"   : paired permutation. For each subject (same subject count
%                             in X1 and X2 required), columns are swapped between X1
%                             and X2 with 50% probability.
%
%   Code created on July 24, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University

    [nTimepoints nSubjects1]=size(matrix1);
    [nTimepoints nSubjects2]=size(matrix2);

    
    if dependentOrIndependent=="independent"
        %permutation for unpaired comparison 
        combinedMatrix = [matrix1 matrix2];
        [nTimepoints nSubjects]=size(combinedMatrix);

        permutedCombinedMatrix = combinedMatrix(:, randperm( nSubjects ));

        permutedMatrix1 = permutedCombinedMatrix(:, 1:nSubjects1);
        permutedMatrix2 = permutedCombinedMatrix(:, nSubjects1+1:end);
        
    elseif dependentOrIndependent=="dependent"        
        %permutation for paired comparison 
        if nSubjects1 ~= nSubjects2
            error("Number of subjects are different")
        end
        
        indexSwapping = rand(1,nSubjects1)>0.5;
        matrix1buffer=matrix1;
        matrix2buffer=matrix2;
        
        matrix1buffer(:,indexSwapping)=matrix2(:,indexSwapping);
        matrix2(:,indexSwapping)= matrix1(:,indexSwapping);
        matrix1=matrix1buffer;
        

        permutedMatrix1=matrix1;
        permutedMatrix2=matrix2;    
    
    else
        error("Only valid arguments : ""dependent"" or ""indepedent"" ")
    end

end

