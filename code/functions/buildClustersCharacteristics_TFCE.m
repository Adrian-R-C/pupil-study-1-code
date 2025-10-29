function clusters = buildClustersCharacteristics_TFCE(pValue, tValue, zValue, dValue, alpha, parametricOrNonParametric)
% buildClustersCharacteristics_TFCE identifies significant temporal clusters 
% using threshold free cluster enhancement and summarizes their statistical properties.
%
%   clusters = buildClustersCharacteristics_TFCE(pValue, tValue, zValue, ...
%               dValue, alpha, parametricOrNonParametric) finds contiguous
%   time regions (clusters) where the p-values are below the significance 
%   threshold alpha. For each cluster, the function computes the mean 
%   statistical value (t or z, depending on the method) and the mean effect 
%   size (Cohen’s d).
%
%   OUTPUT:
%   clusters : structure array with the following fields:
%       .positionCluster : [startIdx endIdx] indices of cluster boundaries
%       .meanTValue      : average t-value within the cluster (parametric mode)
%       .meanZValue      : average z-value within the cluster (non-parametric mode)
%       .meanSizeEffect  : average Cohen’s d within the cluster
%       .pValue          : corrected p-value (0 by default, actually here
%       the p-value is not relevant, it's just more convenient for the rest
%       of the code)
%
%   INPUTS:
%   pValue                 : vector of p-values 
%   tValue                 : vector of t-values [used if parametric]
%   zValue                 : vector of z-values [used if non-parametric]
%   dValue                 : vector of Cohen’s d effect sizes 
%   alpha                  : significance threshold (e.g., 0.05)
%   parametricOrNonParametric : string, either "parametric" or "nonparametric"
%                               to specify which statistic is summarized
%
%   Example:
%       clusters = buildClustersCharacteristics_TFCE(pvals, tvals, zvals, ...
%                   dvals, 0.05, "parametric");
%
%   See also: CBPT_ztValueThreshold_TFCE.
%
%   Code created on September 16, 2025 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.



    significanceMask = pValue < alpha;

    significanceDiff = diff([0 significanceMask 0]); % padding pour capter début et fin
    clusterOnset = find(significanceDiff == 1);
    clusterOffset   = find(significanceDiff == -1) - 1;

    if (~isempty(clusterOnset)) && (clusterOnset(1)==1) && (clusterOffset(1)==1)
        clusterOnset(1)=[];
        clusterOffset(1)=[];
    end

    numberOfClusters = numel(clusterOnset);

        if parametricOrNonParametric=="parametric"

            clusters.positionCluster    =[];
            clusters.meanTValue         =[];
            clusters.meanSizeEffect     =[];
        else
            clusters.positionCluster    =[];
            clusters.meanZValue         =[];
            clusters.meanSizeEffect     =[]; 
        end


    for c = 1:numberOfClusters
        idxRange = clusterOnset(c):clusterOffset(c);

        clusters(c).positionCluster     = [clusterOnset(c) clusterOffset(c)] ;
        if parametricOrNonParametric=="parametric"
                    clusters(c).meanTValue          = mean(tValue(idxRange));
        else
                    clusters(c).meanZValue          = mean(zValue(idxRange));
        end
        clusters(c).meanSizeEffect      = mean(dValue(idxRange));
        clusters(c).pValue              = 0;
    end
end
