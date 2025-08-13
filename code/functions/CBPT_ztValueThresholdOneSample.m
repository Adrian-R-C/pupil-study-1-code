function [uncorrectedSignificance correctedSignificance] = CBPT_ztValueThresholdOneSample(matrix1, pThreshold, nPermutation, tail, parametricOrNonParametric)
% CBPT_ztValueThresholdOneSample computes original clusters in non-permuted data, then 
% finds clusters at each permutation in newly permuted data. Unlike
% CBPT_ztValueThreshold, here the data are first subtracted before performing one-sample tests, but the results are the same, for a sufficiently high number of permutations.  
%
%   [Y1, Y2] = CBPT_ztValueThresholdOneSample(X1, X2, X3, X4, X5) conducts
%   statistical comparison between datasets X1 (subtracted data between two conditions) and a null vector 
%   (one sample testing) to identify significant clusters based on the given
%   p-value threshold (X2), test configuration (only dependent), 
%   directionality (X4), and statistical method (X5).
%   This analysis is then repeated across X3 iterations to establish
%   the permutation distribution.
%   
%   outputs:
%   Y1 : significance vector from unpermuted data
%   Y2 : significance vector after correction
%
%   inputs:
%   X1 : data matrix for condition 1 minus condition 2 (within-subject analysis) 
%   X2 : significance threshold for data comparison
%   X3 : number of permutations
%   X4 : two-tailed or one-tailed test
%   X5 : parametric or non-parametric test
%
%   See also CBPT_ztValueThreshold.
%
%   Code created on September 12, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

%Matrix has to be of size [ nPoints, nSub ]
[numberOfPoints1 numberOfSubject1]=size(matrix1);
numberOfPoints=numberOfPoints1;

significance = zeros(numberOfPoints,1);    

numberOfTimepoints=numberOfPoints1;

% Degrees of freedom (N1 - 1)
degreeOfFreedom = size(matrix1, 2) - 1; 

if string(tail) == "both" 
    tThreshold = tinv(1 - pThreshold, degreeOfFreedom);
    zThreshold = norminv(1 - pThreshold);
elseif (string(tail) == "right")
     tThreshold = tinv(1 - pThreshold, degreeOfFreedom);
     zThreshold = norminv(1 - pThreshold);   
elseif (string(tail) == "left")
    tThreshold = -tinv(1 - pThreshold, degreeOfFreedom);
    zThreshold = -norminv(1 - pThreshold);
else
    error("Only both, right, or left")
end


    %%%%%%%%%%%%%%%%%%%%%
    %%Clusters initiaux%%
    %%%%%%%%%%%%%%%%%%%%%

    [tValueInitial pValueInitialPARAMETRIC zValueInitial pValueInitialNONPARAMETRIC...
        positiveClusterInitialPARAMETRIC negativeClusterInitialPARAMETRIC...
        positiveClusterInitialNONPARAMETRIC negativeClusterInitialNONPARAMETRIC] = extractClusterFromTwoMatricesOneSample(matrix1, tail, abs(tThreshold), abs(zThreshold) );
    
    uncorrectedSignificancePARAMETRIC=zeros(1,numberOfTimepoints);    
    uncorrectedSignificanceNONPARAMETRIC=zeros(1,numberOfTimepoints);
    
    if string(tail) == "both" 
        uncorrectedSignificancePARAMETRIC=abs(tValueInitial)>=tThreshold;
        uncorrectedSignificanceNONPARAMETRIC=abs(zValueInitial)>=zThreshold;
    elseif (string(tail) == "right")
        uncorrectedSignificancePARAMETRIC=tValueInitial>=tThreshold;
        uncorrectedSignificanceNONPARAMETRIC=zValueInitial>=zThreshold;  
    elseif (string(tail) == "left")
        uncorrectedSignificancePARAMETRIC=tValueInitial<=tThreshold;
        uncorrectedSignificanceNONPARAMETRIC=zValueInitial<=zThreshold;
    else
        error("Only both, right, or left")
    end
    
        
    if isempty(positiveClusterInitialPARAMETRIC.length) & isempty(negativeClusterInitialPARAMETRIC.length) & isempty(positiveClusterInitialNONPARAMETRIC.length) & isempty(negativeClusterInitialNONPARAMETRIC.length)
        return
    end
    
    %Preallocation to faster computing
    maxSizeClustersPositivePARAMETRIC=zeros(1,nPermutation);
    maxSizeClustersNegativePARAMETRIC=maxSizeClustersPositivePARAMETRIC;

    numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC=zeros(positiveClusterInitialPARAMETRIC.number,nPermutation);
    numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC=zeros(negativeClusterInitialPARAMETRIC.number,nPermutation);

    significantClusterLengthsPARAMETRIC = [];
    significantPositiveClusterLengthsPARAMETRIC = [];
    significantNegativeClusterLengthsPARAMETRIC = [];

    significantClusterMassesPARAMETRIC = [];
    significantPositiveClusterMassesPARAMETRIC = [];
    significantNegativeClusterMassesPARAMETRIC = [];
    
    maxSizeClustersPositiveNONPARAMETRIC=zeros(1,nPermutation);
    maxSizeClustersNegativeNONPARAMETRIC=maxSizeClustersPositiveNONPARAMETRIC;

    numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC=zeros(positiveClusterInitialNONPARAMETRIC.number,nPermutation);
    numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC=zeros(negativeClusterInitialNONPARAMETRIC.number,nPermutation);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Clusters de permutations%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for cptPermutation=1:nPermutation
        fprintf('Permutation number %d out of %d\n', cptPermutation,nPermutation);
        [permutedMatrix1]=subjectMatrixPermutationOneSample(matrix1);
        
        [tValuePermutation pValuePermutationPARAMETRIC zValuePermutation pValuePermutationNONPARAMETRIC...
         positiveClusterPermutationPARAMETRIC negativeClusterPermutationPARAMETRIC positiveClusterPermutationNONPARAMETRIC negativeClusterPermutationNONPARAMETRIC]...
         = extractClusterFromTwoMatricesOneSample(permutedMatrix1, tail, abs(tThreshold), abs(zThreshold));

     
         maxSizeClustersPositivePARAMETRIC(cptPermutation) = max(positiveClusterPermutationPARAMETRIC.length);
         maxSizeClustersNegativePARAMETRIC(cptPermutation) = max(negativeClusterPermutationPARAMETRIC.length);
         maxSizeClustersPARAMETRIC(cptPermutation) = max(maxSizeClustersPositivePARAMETRIC(cptPermutation), maxSizeClustersNegativePARAMETRIC(cptPermutation));

         maxMassClustersPositivePARAMETRIC(cptPermutation) = max(positiveClusterPermutationPARAMETRIC.tValueSum);
         maxMassClustersNegativePARAMETRIC(cptPermutation) = max(negativeClusterPermutationPARAMETRIC.tValueSum);
         maxMassClustersPARAMETRIC(cptPermutation) = max( abs(maxMassClustersPositivePARAMETRIC(cptPermutation)), abs(maxMassClustersNegativePARAMETRIC(cptPermutation)) );

         maxSizeClustersPositiveNONPARAMETRIC(cptPermutation) = max(positiveClusterPermutationNONPARAMETRIC.length);
         maxSizeClustersNegativeNONPARAMETRIC(cptPermutation) = max(negativeClusterPermutationNONPARAMETRIC.length);
         maxSizeClustersNONPARAMETRIC(cptPermutation) = max(maxSizeClustersPositiveNONPARAMETRIC(cptPermutation), maxSizeClustersNegativeNONPARAMETRIC(cptPermutation));

         maxMassClustersPositiveNONPARAMETRIC(cptPermutation) = max(positiveClusterPermutationNONPARAMETRIC.zValueSum);
         maxMassClustersNegativeNONPARAMETRIC(cptPermutation) = max(negativeClusterPermutationNONPARAMETRIC.zValueSum);
         maxMassClustersNONPARAMETRIC(cptPermutation) = max( abs(maxMassClustersPositiveNONPARAMETRIC(cptPermutation)), abs(maxMassClustersNegativeNONPARAMETRIC(cptPermutation)) );

    end

 
    significancePositivePARAMETRIC=zeros(1,numberOfTimepoints);
    significanceNegativePARAMETRIC=zeros(1,numberOfTimepoints);
 
    significancePositiveNONPARAMETRIC=zeros(1,numberOfTimepoints);
    significanceNegativeNONPARAMETRIC=zeros(1,numberOfTimepoints);
    
    %We define the significance of initial clusters from unpermuted data
    threshold=0.05;
    

if string(tail) == "both"

        %PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialPARAMETRIC.number
            numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersPARAMETRIC > positiveClusterInitialPARAMETRIC.tValueSum(cptInitialClusterPositive) );
    %         numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)=sum( maxSizeClustersPARAMETRIC > positiveClusterInitialPARAMETRIC.length(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)/nPermutation;
            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
%                 pValueFDRPositive(cptInitialClusterPositive)
                significancePositivePARAMETRIC(positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
            end
        end 

        %NON PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialNONPARAMETRIC.number
            numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersNONPARAMETRIC > positiveClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)/nPermutation;
            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
%                 pValueFDRPositive(cptInitialClusterPositive)
                significancePositiveNONPARAMETRIC(positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
            end
        end 

        %PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersPARAMETRIC > abs(negativeClusterInitialPARAMETRIC.tValueSum(cptInitialClusterNegative)) );
    %         numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxSizeClustersPARAMETRIC > negativeClusterInitialPARAMETRIC.length(cptInitialClusterNegative) );
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)/nPermutation;
            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
%                 pValueFDRNegative(cptInitialClusterNegative)
                significanceNegativePARAMETRIC(negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
            end
        end 

        %NON PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialNONPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersNONPARAMETRIC > abs(negativeClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterNegative) ));
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)/nPermutation;
            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
%                 pValueFDRNegative(cptInitialClusterNegative)
                significanceNegativeNONPARAMETRIC(negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
            end
        end            
        
elseif (string(tail) == "right")   

        %PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialPARAMETRIC.number
            %mass
            numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersPositivePARAMETRIC > positiveClusterInitialPARAMETRIC.tValueSum(cptInitialClusterPositive) );
            %length
    %         numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)=sum( maxSizeClustersPositivePARAMETRIC > positiveClusterInitialPARAMETRIC.length(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)/nPermutation;
            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
%                 pValueFDRPositive(cptInitialClusterPositive)
                significancePositivePARAMETRIC(positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
            end 
        end

        %NON PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialNONPARAMETRIC.number
            numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersPositiveNONPARAMETRIC > positiveClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)/nPermutation;
            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
%                 pValueFDRPositive(cptInitialClusterPositive)
                significancePositiveNONPARAMETRIC(positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
            end
        end 

        

elseif (string(tail) == "left")   

        %PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersNegativePARAMETRIC < negativeClusterInitialPARAMETRIC.tValueSum(cptInitialClusterNegative) );
    %         numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxSizeClustersNegativePARAMETRIC > negativeClusterInitialPARAMETRIC.length(cptInitialClusterNegative) );
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)/nPermutation;
            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
%                 pValueFDRNegative(cptInitialClusterNegative)
                significanceNegativePARAMETRIC(negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
            end
        end 


        %NON PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialNONPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersNegativeNONPARAMETRIC < negativeClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterNegative) );
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)/nPermutation;
            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
%                 pValueFDRNegative(cptInitialClusterNegative)
                significanceNegativeNONPARAMETRIC(negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
            end
        end 
    
end
    

    if string(tail) == "both" 
        correctedSignificancePARAMETRIC = or(significancePositivePARAMETRIC,significanceNegativePARAMETRIC);
        correctedSignificanceNONPARAMETRIC = or(significancePositiveNONPARAMETRIC,significanceNegativeNONPARAMETRIC);
    elseif (string(tail) == "right")
        correctedSignificancePARAMETRIC = significancePositivePARAMETRIC;
        correctedSignificanceNONPARAMETRIC = significancePositiveNONPARAMETRIC;
    elseif (string(tail) == "left")
        correctedSignificancePARAMETRIC = significanceNegativePARAMETRIC;
        correctedSignificanceNONPARAMETRIC = significanceNegativeNONPARAMETRIC;
    else
        error("Only both, right, or left")
    end
    
    
    if parametricOrNonParametric=="parametric"
        uncorrectedSignificance = uncorrectedSignificancePARAMETRIC;
        correctedSignificance   = correctedSignificancePARAMETRIC;
    elseif parametricOrNonParametric=="nonparametric"
        uncorrectedSignificance = uncorrectedSignificanceNONPARAMETRIC;
        correctedSignificance   = correctedSignificanceNONPARAMETRIC;
    else %parametric by default
        uncorrectedSignificance = uncorrectedSignificancePARAMETRIC;
        correctedSignificance   = correctedSignificancePARAMETRIC;
    end    
end

