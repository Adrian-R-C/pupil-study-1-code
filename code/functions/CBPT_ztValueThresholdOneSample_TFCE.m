function [uncorrectedSignificance correctedSignificance] = CBPT_ztValueThresholdOneSample_TFCE(matrix1, H, E, dH, nPermutation, tail, parametricOrNonParametric, alphaThreshold)
% CBPT_ztValueThresholdOneSample_TFCE performs a one-sample cluster-based 
%   permutation test (CBPT) using Threshold-Free Cluster Enhancement (TFCE) 
%   on 1D time series data, supporting both parametric (t-tests) and 
%   nonparametric (signed-rank tests) approaches.
%
%   This function tests whether a single condition differs from zero across 
%   multiple timepoints, computes t/z values at each timepoint, applies TFCE, 
%   and evaluates significance via permutation testing.
%
% -------------------------------------------------------------------------
%   INPUTS:
%
%   matrix1 : [nTimepoints x nSubjects]  
%       Data matrix for the single condition or contrast to test.
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
%   tail : string, {"both","right","left"}  
%       Specifies the type of test: two-tailed, right-tailed, or left-tailed.
%
%   parametricOrNonParametric : string, {"parametric","nonparametric"}  
%       Defines which test to use for final significance evaluation:
%       - "parametric" : uses t-tests
%       - "nonparametric" : uses Wilcoxon signed-rank tests
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
%       Vector marking significant timepoints after TFCE-based permutation 
%       correction.
%
%   Created on September 16, 2025
%   Author: Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition (CerCo)
%   CNRS / Toulouse University



%Matrix has to be of size [ nPoints, nSub ]
[numberOfTimepoints numberOfSubject]=size(matrix1);

% Degrees of freedom (N1 - 1)
degreeOfFreedom = size(matrix1, 2) - 1; 

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
    %%Clusters initiaux%%
    %%%%%%%%%%%%%%%%%%%%%

    tValueInitial=zeros(1,numberOfTimepoints);  
    zValueInitial=zeros(1,numberOfTimepoints);
    sizeEffect=zeros(1,numberOfTimepoints);
    uncorrectedSignificancePARAMETRIC=zeros(1,numberOfTimepoints);    
    uncorrectedSignificanceNONPARAMETRIC=zeros(1,numberOfTimepoints);
     
    
    distributionAtChosenPoint1=[];
    
    for cptTimepoint=1:numberOfTimepoints

        distributionAtChosenPoint1 = matrix1(cptTimepoint, :);

        %We remove nan columns
        whoIsNan=isnan(distributionAtChosenPoint1);

        distributionAtChosenPoint1=distributionAtChosenPoint1(~whoIsNan);

        if numel(distributionAtChosenPoint1)>1
                
                    %PARAMETRIC
                    [~, pValuePARAMETRIC(cptTimepoint), ~, statsPARAMETRIC] = ttest(distributionAtChosenPoint1, zeros(1, numel(distributionAtChosenPoint1)),'tail', tail);
                    tValueInitial(cptTimepoint) = statsPARAMETRIC.tstat;
                    sizeEffect(cptTimepoint)=computeCohen_d(distributionAtChosenPoint1, zeros(1, numel(distributionAtChosenPoint1)));

                    %NON PARAMETRIC
                    [pValueNONPARAMETRIC(cptTimepoint), ~, statsNONPARAMETRIC]=signrank(distributionAtChosenPoint1, zeros(1, numel(distributionAtChosenPoint1)), 'tail', tail);
                    if isfield(statsNONPARAMETRIC, "zval")
                        zValueInitial(cptTimepoint)=statsNONPARAMETRIC.zval;
                    else 
                        zValueInitial(cptTimepoint)=0;
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
        [permutedMatrix1]=subjectMatrixPermutationOneSample(matrix1);

        tValuePermutation = zeros(1,numberOfTimepoints);
        zValuePermutation = zeros(1,numberOfTimepoints);

        for cptTimepoint=1:numberOfTimepoints

            distributionAtChosenPoint1 = permutedMatrix1(cptTimepoint, :);

            %We remove nan columns
            whoIsNan=isnan(distributionAtChosenPoint1);

            distributionAtChosenPoint1=distributionAtChosenPoint1(~whoIsNan);

            if numel(distributionAtChosenPoint1)>1

                        %PARAMETRIC
                        [~, pValuePARAMETRIC(cptTimepoint), ~, statsPARAMETRIC] = ttest(distributionAtChosenPoint1, zeros(1, numel(distributionAtChosenPoint1)),'tail', tail);
                        tValuePermutation(cptTimepoint) = statsPARAMETRIC.tstat;

                        %NON PARAMETRIC
                        [pValueNONPARAMETRIC(cptTimepoint), ~, statsNONPARAMETRIC]=signrank(distributionAtChosenPoint1, zeros(1, numel(distributionAtChosenPoint1)), 'tail', tail);
                        if isfield(statsNONPARAMETRIC, "zval")
                            zValuePermutation(cptTimepoint)=statsNONPARAMETRIC.zval;
                        else 
                            zValuePermutation(cptTimepoint)=0;
                        end


            else
                tValuePermutation(cptTimepoint)=0;
                zValuePermutation(cptTimepoint)=0;
                pValuePARAMETRIC(cptTimepoint)=NaN;
                pValueNONPARAMETRIC(cptTimepoint)=NaN;
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
    for timepoint = 1:numberOfTimepoints
        switch tail
            case "both"
                correctedPValuePARAMETRIC(timepoint)      = ( sum(thresholdTFCE_PARAMETRIC      >= abs(TFCEscoreInitialPARAMETRIC(timepoint))))    / nPermutation;
                correctedPValueNONPARAMETRIC(timepoint)   = ( sum(thresholdTFCE_NONPARAMETRIC   >= abs(TFCEscoreInitialNONPARAMETRIC(timepoint)))) / nPermutation;
            case "right"
                correctedPValuePARAMETRIC(timepoint)      = ( sum(thresholdTFCE_PARAMETRIC      >= TFCEscoreInitialPARAMETRIC(timepoint)))         / nPermutation;
                correctedPValueNONPARAMETRIC(timepoint)   = ( sum(thresholdTFCE_NONPARAMETRIC   >= TFCEscoreInitialNONPARAMETRIC(timepoint)))      / nPermutation;
            case "left"
                correctedPValuePARAMETRIC(timepoint)      = ( sum(thresholdTFCE_PARAMETRIC      <= TFCEscoreInitialPARAMETRIC(timepoint)))         / nPermutation;
                correctedPValueNONPARAMETRIC(timepoint)   = ( sum(thresholdTFCE_NONPARAMETRIC   <= TFCEscoreInitialNONPARAMETRIC(timepoint)))      / nPermutation;
        end
    end


correctedSignificancePARAMETRIC    = correctedPValuePARAMETRIC < alphaThreshold;
correctedSignificanceNONPARAMETRIC = correctedPValueNONPARAMETRIC < alphaThreshold;

% % Now we compute the cluster characteristics
% clusterCharacteristics = buildClustersCharacteristics_TFCE(correctedPValuePARAMETRIC, tValueInitial, zValueInitial, sizeEffect, alphaThreshold, parametricOrNonParametric);

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

