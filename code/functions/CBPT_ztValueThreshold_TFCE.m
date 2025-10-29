function [uncorrectedSignificance correctedSignificance clusterCharacteristics] = CBPT_ztValueThreshold_TFCE(matrix1, matrix2, H, E, dH, nPermutation, dependentOrIndependent, tail, parametricOrNonParametric, alphaThreshold)
% CBPT_ztValueThreshold_TFCE performs a cluster-based permutation test (CBPT) with Threshold-Free
%   Cluster Enhancement (TFCE) on 1D time series data (here pupillometry), using both
%   parametric (t-tests) and nonparametric (rank tests) approaches.
%
%   This function compares two conditions (within-subject or between-group)
%   across multiple timepoints, computes t/z values, applies TFCE, and
%   corrects significance through permutation testing.
%
% -------------------------------------------------------------------------
%   INPUTS:
%
%   matrix1 : [nTimepoints x nSubjects1]  
%       Data matrix for condition 1 (within-subject) or group 1 (between-subject).
%
%   matrix2 : [nTimepoints x nSubjects2]  
%       Data matrix for condition 2 (within-subject) or group 2 (between-subject).
%
%   H : scalar  
%       TFCE parameter controlling the weight of the cluster height.
%
%   E : scalar  
%       TFCE parameter controlling the weight of the cluster extent.
%
%   dH : scalar  
%       TFCE step size for the integration of cluster heights.
%
%   nPermutation : integer  
%       Number of random permutations to compute the null distribution.
%
%   dependentOrIndependent : string, {"dependent","independent"}  
%       Defines whether the comparison is paired (within-subject) or unpaired (between-groups).
%
%   tail : string, {"both","right","left"}  
%       Specifies the type of test: two-tailed, right-tailed, or left-tailed.
%
%   parametricOrNonParametric : string, {"parametric","nonparametric"}  
%       Defines which test to use for final significance evaluation:
%       - "parametric" : uses t-tests
%       - "nonparametric" : uses Wilcoxon/Mann-Whitney
%
%   alphaThreshold : scalar (0 < alpha < 1)  
%       Significance level threshold (e.g., 0.05).
%
% -------------------------------------------------------------------------
%   OUTPUTS:
%
%   uncorrectedSignificance : [nTimepoints x 1] logical  
%       Vector marking significant timepoints without permutation correction.
%
%   correctedSignificance : [nTimepoints x 1] logical  
%       Vector marking significant timepoints after cluster-based permutation correction.
%
%   clusterCharacteristics : struct  
%
%       Information about detected clusters, including:
%       .positionCluster : [startIdx endIdx] indices of cluster boundaries
%       .meanTValue      : average t-value within the cluster (parametric mode)
%       .meanZValue      : average z-value within the cluster (non-parametric mode)
%       .meanSizeEffect  : average Cohen’s d within the cluster
%       .pValue          : corrected p-value (0 by default, actually here
%       the p-value is not relevant, it's just more convenient for the rest
%       of the code)
%
%   Created on September 16, 2025
%   Author: Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition (CerCo)
%   CNRS / Toulouse University


%Matrix has to be of size [ nPoints, nSub ]
[numberOfPoints1 numberOfSubject1]=size(matrix1);
[numberOfPoints2 numberOfSubject2]=size(matrix2);

numberOfTimepoints=numberOfPoints1;

if numberOfPoints1 ~= numberOfPoints2
    error("Number of timepoints are different")
end

if dependentOrIndependent=="dependent"
    if numberOfSubject1 ~= numberOfSubject2
        error("Number of subjects are different, not possible for a paired comparison")
    end
end   

if dependentOrIndependent~="dependent" && dependentOrIndependent~="independent"
    error("Only valid arguments : ""dependent"" or ""independent"" ")
end


if dependentOrIndependent=="dependent"
    % Degrees of freedom (N1 - 1)
    degreeOfFreedom = size(matrix1, 2) - 1; 
    
else
    % Degrees of freedom (N1 + N2 - 2)
    degreeOfFreedom = size(matrix1, 2) + size(matrix2, 2) - 2; 
end   


% Degrees of freedom 
% Since we use TFCE, this will be used only for uncorrected t/p time
% courses
if string(tail) == "both" 
    tThreshold = tinv(1 - alphaThreshold/2, degreeOfFreedom);
    zThreshold = norminv(1 - alphaThreshold/2);
elseif (string(tail) == "right")
     tThreshold = tinv(1 - alphaThreshold, degreeOfFreedom);
     zThreshold = norminv(1 - alphaThreshold);   
elseif (string(tail) == "left")
    tThreshold = -tinv(1 - alphaThreshold, degreeOfFreedom);
    zThreshold = -norminv(1 - alphaThreshold);
else
    error("Only both, right, or left")
end

    %%%%%%%%%%%%%%%%%%%%%
    %%Initial clusters %%
    %%%%%%%%%%%%%%%%%%%%%

    tValueInitial=zeros(1,numberOfTimepoints);  
    zValueInitial=zeros(1,numberOfTimepoints);  
    sizeEffect=zeros(1,numberOfTimepoints);
    uncorrectedSignificancePARAMETRIC=zeros(1,numberOfTimepoints);    
    uncorrectedSignificanceNONPARAMETRIC=zeros(1,numberOfTimepoints);


    for cptTimepoint = 1:numberOfTimepoints

        distributionAtChosenPoint1 = matrix1(cptTimepoint, :);
        distributionAtChosenPoint2 = matrix2(cptTimepoint, :);

        %We remove nan columns
        distributionAtChosenPoint1=distributionAtChosenPoint1(~isnan(distributionAtChosenPoint1));
        distributionAtChosenPoint2=distributionAtChosenPoint2(~isnan(distributionAtChosenPoint2));

        if numel(distributionAtChosenPoint1) > 1 && numel(distributionAtChosenPoint2) > 1
                if dependentOrIndependent=="independent" %unpaired


                        [~,~,~,statsPARAMETRIC]     = ttest2(distributionAtChosenPoint1, distributionAtChosenPoint2);
                        tValueInitial(cptTimepoint)    = statsPARAMETRIC.tstat;
                        sizeEffect(cptTimepoint)=computeCohen_d(distributionAtChosenPoint1, distributionAtChosenPoint2);

                        [~,~,statsNONPARAMETRIC]  = ranksum(distributionAtChosenPoint1, distributionAtChosenPoint2);
                        if isfield(statsNONPARAMETRIC, "zval")
                            zValueInitial(cptTimepoint)=statsNONPARAMETRIC.zval;
                        else 
                            zValueInitial(cptTimepoint)=0;
                        end

                elseif dependentOrIndependent=="dependent" %paired

                        [~,~,~,statsPARAMETRIC]     = ttest(distributionAtChosenPoint1, distributionAtChosenPoint2);
                        tValueInitial(cptTimepoint)    = statsPARAMETRIC.tstat;
                        sizeEffect(cptTimepoint)=computeCohen_d(distributionAtChosenPoint1, distributionAtChosenPoint2);

                        [~,~,statsNONPARAMETRIC]  = signrank(distributionAtChosenPoint1, distributionAtChosenPoint2);
                        if isfield(statsNONPARAMETRIC, "zval")
                            zValueInitial(cptTimepoint)=statsNONPARAMETRIC.zval;
                        else 
                            zValueInitial(cptTimepoint)=0;
                        end

                else
                        error("Only valid arguments : ""dependent"" or ""independent"" ")
                end
        

        else 

            tValueInitial(cptTimepoint)=0;
            zValueInitial(cptTimepoint)=0;
            sizeEffect(cptTimepoint)=0;

    end

    end
    
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

        
    %We then convert the t-values/z-values to TFCE score. This score
    %replaces the mass or length in conventional CBPT
    TFCEscoreInitialPARAMETRIC      =  tfceScore_1D(tValueInitial, H, E, dH);
    TFCEscoreInitialNONPARAMETRIC   =  tfceScore_1D(zValueInitial, H, E, dH);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Clusters from permuted data%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for cptPermutation=1:nPermutation
        fprintf('Permutation number %d out of %d\n', cptPermutation,nPermutation);
        [permutedMatrix1 permutedMatrix2]=subjectMatrixPermutation(matrix1, matrix2, dependentOrIndependent);
        
        tValuePermutation = zeros(1,numberOfTimepoints);
        zValuePermutation = zeros(1,numberOfTimepoints);

        for cptTimepoint = 1:numberOfTimepoints

            distributionAtChosenPoint1 = permutedMatrix1(cptTimepoint, :);
            distributionAtChosenPoint2 = permutedMatrix2(cptTimepoint, :);

            %We remove nan columns
            distributionAtChosenPoint1=distributionAtChosenPoint1(~isnan(distributionAtChosenPoint1));
            distributionAtChosenPoint2=distributionAtChosenPoint2(~isnan(distributionAtChosenPoint2));

            if numel(distributionAtChosenPoint1) > 1 && numel(distributionAtChosenPoint2) > 1
                if dependentOrIndependent=="independent" %unpaired

                [~,~,~,statsPARAMETRIC]     = ttest2(distributionAtChosenPoint1, distributionAtChosenPoint2);
                tValuePermutation(cptTimepoint)    = statsPARAMETRIC.tstat;

                [~,~,statsNONPARAMETRIC]  = ranksum(distributionAtChosenPoint1, distributionAtChosenPoint2);
                    if isfield(statsNONPARAMETRIC, "zval")
                    zValuePermutation(cptTimepoint)=statsNONPARAMETRIC.zval;
                    else 
                    zValuePermutation(cptTimepoint)=0;
                    end

                elseif dependentOrIndependent=="dependent" %paired

                [~,~,~,statsPARAMETRIC]     = ttest(distributionAtChosenPoint1, distributionAtChosenPoint2);
                tValuePermutation(cptTimepoint)    = statsPARAMETRIC.tstat;

                [~,~,statsNONPARAMETRIC]  = signrank(distributionAtChosenPoint1, distributionAtChosenPoint2);
                    if isfield(statsNONPARAMETRIC, "zval")
                    zValuePermutation(cptTimepoint)=statsNONPARAMETRIC.zval;
                    else 
                    zValuePermutation(cptTimepoint)=0;
                    end

                else
                error("Only valid arguments : ""dependent"" or ""independent"" ")
                end


            else 

            tValuePermutation(cptTimepoint)=0;
            zValuePermutation(cptTimepoint)=0;

            end

        end



         TFCEscoreClustersPARAMETRIC      =  tfceScore_1D(tValuePermutation, H, E, dH);
         TFCEscoreClustersNONPARAMETRIC   =  tfceScore_1D(zValuePermutation, H, E, dH);
     
        switch tail
            case "both"
                thresholdTFCE_PARAMETRIC(cptPermutation) = max(abs(TFCEscoreClustersPARAMETRIC)); % two-tailed
                thresholdTFCE_NONPARAMETRIC(cptPermutation) = max(abs(TFCEscoreClustersNONPARAMETRIC)); % two-tailed
            case "right"
                thresholdTFCE_PARAMETRIC(cptPermutation) = max(TFCEscoreClustersPARAMETRIC); % one-tailed positive
                thresholdTFCE_NONPARAMETRIC(cptPermutation) = max(TFCEscoreClustersNONPARAMETRIC); % one-tailed positive
            case "left"
                thresholdTFCE_PARAMETRIC(cptPermutation) = min(TFCEscoreClustersPARAMETRIC); % one-tailed negative
                thresholdTFCE_NONPARAMETRIC(cptPermutation) = min(TFCEscoreClustersNONPARAMETRIC); % one-tailed negative
            otherwise
                error("Only both, right, or left");

        end
    end

% Computation of the new p-value vector (i.e. corrected significance)
    correctedPValuePARAMETRIC     = zeros(1,numberOfTimepoints); 
    correctedPValueNONPARAMETRIC  = zeros(1,numberOfTimepoints); 
    for cptTimepoint = 1:numberOfTimepoints
        switch tail
            case "both"
                correctedPValuePARAMETRIC(cptTimepoint)      = ( sum(thresholdTFCE_PARAMETRIC      >= abs(TFCEscoreInitialPARAMETRIC(cptTimepoint))))    / nPermutation;
                correctedPValueNONPARAMETRIC(cptTimepoint)   = ( sum(thresholdTFCE_NONPARAMETRIC   >= abs(TFCEscoreInitialNONPARAMETRIC(cptTimepoint)))) / nPermutation;
            case "right"
                correctedPValuePARAMETRIC(cptTimepoint)      = ( sum(thresholdTFCE_PARAMETRIC      >= TFCEscoreInitialPARAMETRIC(cptTimepoint)))         / nPermutation;
                correctedPValueNONPARAMETRIC(cptTimepoint)   = ( sum(thresholdTFCE_NONPARAMETRIC   >= TFCEscoreInitialNONPARAMETRIC(cptTimepoint)))      / nPermutation;
            case "left"
                correctedPValuePARAMETRIC(cptTimepoint)      = ( sum(thresholdTFCE_PARAMETRIC      <= TFCEscoreInitialPARAMETRIC(cptTimepoint)))         / nPermutation;
                correctedPValueNONPARAMETRIC(cptTimepoint)   = ( sum(thresholdTFCE_NONPARAMETRIC   <= TFCEscoreInitialNONPARAMETRIC(cptTimepoint)))      / nPermutation;
        end
    end


correctedSignificancePARAMETRIC    = correctedPValuePARAMETRIC < alphaThreshold;
correctedSignificanceNONPARAMETRIC = correctedPValueNONPARAMETRIC < alphaThreshold;

% Now we compute the cluster characteristics and return the corrected
% clusters

    if parametricOrNonParametric=="parametric"
        correctedPValue=correctedPValuePARAMETRIC;
    elseif parametricOrNonParametric=="nonparametric"
        correctedPValue=correctedPValueNONPARAMETRIC;
    else %parametric by default
        correctedPValue=correctedPValuePARAMETRIC;
    end   
    clusterCharacteristics = buildClustersCharacteristics_TFCE(correctedPValue, tValueInitial, zValueInitial, sizeEffect, alphaThreshold, parametricOrNonParametric);


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

    uncorrectedSignificance(1)=0;
    correctedSignificance(1)=0;
    
end

