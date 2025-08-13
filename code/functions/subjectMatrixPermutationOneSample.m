function [permutedMatrix1]=subjectMatrixPermutationOneSample(matrix1)
% subjectMatrixPermutationOneSample applies sign-flipping permutation for one-sample tests.
%
%   Y = subjectMatrixPermutationOneSample(X) randomly flips the sign of each
%   subject's data column with 50% probability. 
%
%   output:
%     Y : numeric matrix (nTimepoints × nSubjects) where each subject's column
%         has been multiplied by either +1 or -1, chosen at random.
%
%   input:
%     X : numeric matrix (nTimepoints × nSubjects) containing the original data.
%
%   See also subjectMatrixPermutation
%
%   Code created on August 8, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University

    [nTimepoints nSubjects]=size(matrix1);
    %permutation for paired comparison
    permutationCoefficients = -1+2*(rand(1,nSubjects)>0.5);
    permutedMatrix1=permutationCoefficients.*matrix1;

end

