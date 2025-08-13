function [uncorrectedSignificance correctedSignificance clusterCharacteristics] = CBPT_ztValueThreshold(matrix1, matrix2, pThreshold, nPermutation, dependentOrIndependent, tail, parametricOrNonParametric)
% CBPT_ztValueThreshold computes original clusters in non-permuted data, then 
% finds clusters at each permutation in newly permuted data.
%
%   [Y1, Y2, Y3] = CBPT_ztValueThreshold(X1, X2, X3, X4, X5, X6, X7) conducts
%   statistical comparison between datasets X1 (condition/group 1) and X2 
%   (condition/group 2) to identify significant clusters based on the given
%   p-value threshold (X3), test configuration (X5), 
%   directionality (X6), and statistical method (X7).
%   This analysis is then repeated across X4 iterations to establish
%   the permutation distribution.
%   
%   outputs:
%   Y1 : significance vector from unpermuted data
%   Y2 : significance vector after correction
%   Y3 : characteristics of each cluster
%
%   inputs:
%   X1 : data matrix for condition 1 (within-subject analysis) or group 1
%   (between-subject analysis) 
%   X2 : data matrix for condition 2 (within-subject analysis) or group 2
%   (between-subject analysis) 
%   X3 : significance threshold for data comparison
%   X4 : number of permutations
%   X5 : string that determines whether to perform a dependent
%   (within-subject) or independent (between-subject) test 
%   X6 : two-tailed or one-tailed test
%   X7 : parametric or non-parametric test
%
%   See also CBPT_ztValueThresholdOneSample.
%
%   Code created on September 10, 2024 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

%Matrix has to be of size [ nPoints, nSub ]
[numberOfPoints1 numberOfSubject1]=size(matrix1);
[numberOfPoints2 numberOfSubject2]=size(matrix2);

numberOfPoints=numberOfPoints1;
significance = zeros(numberOfPoints,1);    

if numberOfPoints1 ~= numberOfPoints2
    error("Number of timepoints are different")
end

if dependentOrIndependent=="dependent"
    if numberOfSubject1 ~= numberOfSubject2
        error("Number of subjects are different, not possible for a paired sign rank")
    end
end   

if dependentOrIndependent~="dependent" && dependentOrIndependent~="independent"
    error("Only valid arguments : ""dependent"" or ""indepedent"" ")
end

 
numberOfTimepoints=numberOfPoints1;

if dependentOrIndependent=="dependent"
    % Degrees of freedom (N1 - 1)
    degreeOfFreedom = size(matrix1, 2) - 1; 
    
else
    % Degrees of freedom (N1 + N2 - 2)
    degreeOfFreedom = size(matrix1, 2) + size(matrix2, 2) - 2; 
    
end   

%pThreshold fixed at 0.025
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

%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    %%Initial clusters %%
    %%%%%%%%%%%%%%%%%%%%%

    [tValueInitial pValueInitialPARAMETRIC sizeEffectPARAMETRIC zValueInitial pValueInitialNONPARAMETRIC...
        positiveClusterInitialPARAMETRIC negativeClusterInitialPARAMETRIC...
        positiveClusterInitialNONPARAMETRIC negativeClusterInitialNONPARAMETRIC] = extractClusterFromTwoMatrices(matrix1, matrix2, dependentOrIndependent, tail, abs(tThreshold), abs(zThreshold));
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Clusters from permuted data%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for cptPermutation=1:nPermutation
        fprintf('Permutation number %d out of %d\n', cptPermutation,nPermutation);
        [permutedMatrix1 permutedMatrix2]=subjectMatrixPermutation(matrix1, matrix2, dependentOrIndependent);
        
        [tValuePermutation pValuePermutationPARAMETRIC sizeEffectPARAMETRIC zValuePermutation pValuePermutationNONPARAMETRIC...
         positiveClusterPermutationPARAMETRIC negativeClusterPermutationPARAMETRIC positiveClusterPermutationNONPARAMETRIC negativeClusterPermutationNONPARAMETRIC]...
         = extractClusterFromTwoMatrices(permutedMatrix1, permutedMatrix2, dependentOrIndependent, tail, abs(tThreshold), abs(zThreshold) );

         
         maxSizeClustersPositivePARAMETRIC(cptPermutation) = max(positiveClusterPermutationPARAMETRIC.length);
         maxSizeClustersNegativePARAMETRIC(cptPermutation) = max(negativeClusterPermutationPARAMETRIC.length);
         maxSizeClustersPARAMETRIC(cptPermutation) = max(maxSizeClustersPositivePARAMETRIC(cptPermutation), maxSizeClustersNegativePARAMETRIC(cptPermutation));
        
         maxMassClustersPositivePARAMETRIC(cptPermutation) = max(positiveClusterPermutationPARAMETRIC.tValueSum);
         maxMassClustersNegativePARAMETRIC(cptPermutation) = min(negativeClusterPermutationPARAMETRIC.tValueSum);
         maxMassClustersPARAMETRIC(cptPermutation) = max( abs(maxMassClustersPositivePARAMETRIC(cptPermutation)), abs(maxMassClustersNegativePARAMETRIC(cptPermutation)) );
         
%             significantPositiveClusterLengthsPARAMETRIC = [significantPositiveClusterLengthsPARAMETRIC, positiveClusterPermutationPARAMETRIC.length];
%             significantNegativeClusterLengthsPARAMETRIC = [significantNegativeClusterLengthsPARAMETRIC, negativeClusterPermutationPARAMETRIC.length];
% 
%             significantPositiveClusterMassesPARAMETRIC = [];
%             significantNegativeClusterMassesPARAMETRIC = [];


         maxSizeClustersPositiveNONPARAMETRIC(cptPermutation) = max(positiveClusterPermutationNONPARAMETRIC.length);
         maxSizeClustersNegativeNONPARAMETRIC(cptPermutation) = max(negativeClusterPermutationNONPARAMETRIC.length);
         maxSizeClustersNONPARAMETRIC(cptPermutation) = max(maxSizeClustersPositiveNONPARAMETRIC(cptPermutation), maxSizeClustersNegativeNONPARAMETRIC(cptPermutation));

         maxMassClustersPositiveNONPARAMETRIC(cptPermutation) = max(positiveClusterPermutationNONPARAMETRIC.zValueSum);
         maxMassClustersNegativeNONPARAMETRIC(cptPermutation) = min(negativeClusterPermutationNONPARAMETRIC.zValueSum);
         maxMassClustersNONPARAMETRIC(cptPermutation) = max( abs(maxMassClustersPositiveNONPARAMETRIC(cptPermutation)), abs(maxMassClustersNegativeNONPARAMETRIC(cptPermutation)) );

%           tValueMatrix(:,cptPermutation)=tValuePermutation;
          
    end
%             significantClusterLengthsPARAMETRIC = [significantPositiveClusterLengthsPARAMETRIC, significantNegativeClusterLengthsPARAMETRIC];
%             significantClusterMassesPARAMETRIC =  [significantPositiveClusterMassesPARAMETRIC, significantNegativeClusterMassesPARAMETRIC];

 
    significancePositivePARAMETRIC=zeros(1,numberOfTimepoints);
    significanceNegativePARAMETRIC=zeros(1,numberOfTimepoints);
 
    significancePositiveNONPARAMETRIC=zeros(1,numberOfTimepoints);
    significanceNegativeNONPARAMETRIC=zeros(1,numberOfTimepoints);
    
    iteratorSignificantClusterPARAMETRIC=0;
    clusterPARAMETRIC_characteristics.positionCluster    =[];
    clusterPARAMETRIC_characteristics.meanTValue         =[];
    clusterPARAMETRIC_characteristics.meanSizeEffect     =[];

    iteratorSignificantClusterNONPARAMETRIC=0;
    clusterNONPARAMETRIC_characteristics.positionCluster    =[];
    clusterNONPARAMETRIC_characteristics.meanZValue         =[];
    clusterNONPARAMETRIC_characteristics.meanSizeEffect     =[];    
    

    %We define the significance of initial clusters from unpermuted data
    threshold=0.05;
    
if string(tail) == "both"

            %PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialPARAMETRIC.number
            numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersPARAMETRIC > positiveClusterInitialPARAMETRIC.tValueSum(cptInitialClusterPositive) );
    %         numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)=sum( maxSizeClustersPARAMETRIC > positiveClusterInitialPARAMETRIC.length(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterPositive)/nPermutation;

            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
                significancePositivePARAMETRIC(positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
                iteratorSignificantClusterPARAMETRIC=iteratorSignificantClusterPARAMETRIC+1;
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).positionCluster = positiveClusterInitialPARAMETRIC.position ;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanTValue =      positiveClusterInitialPARAMETRIC.tValueMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanSizeEffect =  positiveClusterInitialPARAMETRIC.sizeEffectMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).pValue =          pValueFDRPositive(cptInitialClusterPositive);    

            end
        end 


        %NON PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialNONPARAMETRIC.number
            numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersNONPARAMETRIC > positiveClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)/nPermutation;

            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
                significancePositiveNONPARAMETRIC(positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
                iteratorSignificantClusterNONPARAMETRIC=iteratorSignificantClusterNONPARAMETRIC+1;
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).positionCluster = positiveClusterInitialNONPARAMETRIC.position ;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanZValue =      positiveClusterInitialNONPARAMETRIC.zValueMean;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanSizeEffect =  positiveClusterInitialNONPARAMETRIC.sizeEffectMean;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).pValue =          pValueFDRPositive(cptInitialClusterPositive);                   
            end
        end 

    
        %PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersPARAMETRIC > abs(negativeClusterInitialPARAMETRIC.tValueSum(cptInitialClusterNegative)) );
    %         numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxSizeClustersPARAMETRIC > negativeClusterInitialPARAMETRIC.length(cptInitialClusterNegative) );
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)/nPermutation;

            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
                significanceNegativePARAMETRIC(negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
                iteratorSignificantClusterPARAMETRIC=iteratorSignificantClusterPARAMETRIC+1;
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).positionCluster = negativeClusterInitialPARAMETRIC.position ;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanTValue =      negativeClusterInitialPARAMETRIC.tValueMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanSizeEffect =  negativeClusterInitialPARAMETRIC.sizeEffectMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).pValue =          pValueFDRNegative(cptInitialClusterNegative);   
            end
        end 


        %NON PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialNONPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersNONPARAMETRIC > abs(negativeClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterNegative) ));
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)/nPermutation;

            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
                significanceNegativeNONPARAMETRIC(negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
                iteratorSignificantClusterNONPARAMETRIC=iteratorSignificantClusterNONPARAMETRIC+1;
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).positionCluster = negativeClusterInitialNONPARAMETRIC.position ;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanZValue =      negativeClusterInitialNONPARAMETRIC.zValueMean;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanSizeEffect =  negativeClusterInitialNONPARAMETRIC.sizeEffectMean; 
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).pValue =          pValueFDRNegative(cptInitialClusterNegative); 
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
                significancePositivePARAMETRIC(positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
                iteratorSignificantClusterPARAMETRIC=iteratorSignificantClusterPARAMETRIC+1;
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).positionCluster = positiveClusterInitialPARAMETRIC.position ;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanTValue =      positiveClusterInitialPARAMETRIC.tValueMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanSizeEffect =  positiveClusterInitialPARAMETRIC.sizeEffectMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).pValue =          pValueFDRPositive(cptInitialClusterPositive);   
            end
        end 


        %NON PARAMETRIC
        for cptInitialClusterPositive=1:positiveClusterInitialNONPARAMETRIC.number
            numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)=sum( maxMassClustersPositiveNONPARAMETRIC > positiveClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterPositive) );
            pValueFDRPositive(cptInitialClusterPositive)=numberOfPositiveClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterPositive)/nPermutation;

            if pValueFDRPositive(cptInitialClusterPositive) <= threshold
                significancePositiveNONPARAMETRIC(positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(1):positiveClusterInitialNONPARAMETRIC.position{cptInitialClusterPositive}(2))=1;
                iteratorSignificantClusterNONPARAMETRIC=iteratorSignificantClusterNONPARAMETRIC+1;
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).positionCluster = positiveClusterInitialNONPARAMETRIC.position ;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanZValue =      positiveClusterInitialNONPARAMETRIC.zValueMean;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanSizeEffect =  positiveClusterInitialNONPARAMETRIC.sizeEffectMean;
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).pValue =          pValueFDRPositive(cptInitialClusterPositive);   
            end
        end 

        

elseif (string(tail) == "left")   

        %PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersNegativePARAMETRIC < negativeClusterInitialPARAMETRIC.tValueSum(cptInitialClusterNegative) );
    %         numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)=sum( maxSizeClustersNegativePARAMETRIC > negativeClusterInitialPARAMETRIC.length(cptInitialClusterNegative) );
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOnePARAMETRIC(cptInitialClusterNegative)/nPermutation;

            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
                significanceNegativePARAMETRIC(negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
                iteratorSignificantClusterPARAMETRIC=iteratorSignificantClusterPARAMETRIC+1;
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).positionCluster = negativeClusterInitialPARAMETRIC.position ;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanTValue =      negativeClusterInitialPARAMETRIC.tValueMean;   
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).meanSizeEffect =  negativeClusterInitialPARAMETRIC.sizeEffectMean;  
                clusterPARAMETRIC_characteristics(iteratorSignificantClusterPARAMETRIC).pValue =          pValueFDRNegative(cptInitialClusterNegative);   
            end
        end 


        %NON PARAMETRIC
        for cptInitialClusterNegative=1:negativeClusterInitialNONPARAMETRIC.number
            numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)=sum( maxMassClustersNegativeNONPARAMETRIC < negativeClusterInitialNONPARAMETRIC.zValueSum(cptInitialClusterNegative) );
            pValueFDRNegative(cptInitialClusterNegative)=numberOfNegativeClusterFoundBiggerThanInitialOneNONPARAMETRIC(cptInitialClusterNegative)/nPermutation;

            if pValueFDRNegative(cptInitialClusterNegative) <= threshold
                significanceNegativeNONPARAMETRIC(negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(1):negativeClusterInitialNONPARAMETRIC.position{cptInitialClusterNegative}(2))=1;
                iteratorSignificantClusterNONPARAMETRIC=iteratorSignificantClusterNONPARAMETRIC+1;
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).positionCluster = negativeClusterInitialNONPARAMETRIC.position ;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanZValue =      negativeClusterInitialNONPARAMETRIC.zValueMean;   
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).meanSizeEffect =  negativeClusterInitialNONPARAMETRIC.sizeEffectMean; 
                clusterNONPARAMETRIC_characteristics(iteratorSignificantClusterNONPARAMETRIC).pValue =          pValueFDRNegative(cptInitialClusterNegative);
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
        clusterCharacteristics  = clusterPARAMETRIC_characteristics;
    elseif parametricOrNonParametric=="nonparametric"
        uncorrectedSignificance = uncorrectedSignificanceNONPARAMETRIC;
        correctedSignificance   = correctedSignificanceNONPARAMETRIC;
        clusterCharacteristics  = clusterNONPARAMETRIC_characteristics;
    else %parametric by default
        uncorrectedSignificance = uncorrectedSignificancePARAMETRIC;
        correctedSignificance   = correctedSignificancePARAMETRIC;
        clusterCharacteristics  = clusterPARAMETRIC_characteristics;
    end    
    
end

