function [tValue pValuePARAMETRIC sizeEffectPARAMETRIC zValue pValueNONPARAMETRIC positiveClusterPARAMETRIC negativeClusterPARAMETRIC positiveClusterNONPARAMETRIC negativeClusterNONPARAMETRIC] = extractClusterFromTwoMatrices(matrix1, matrix2, dependentOrIndependent, tail, tThreshold, zThreshold)
% extractClusterFromTwoMatrices performs point-by-point statistical comparison and extract clusters of significant time points.
%
% [Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9] = extractClusterFromTwoMatrices(X1, X2, X3, X4, X5, X6)
% compares two time-by-subject matrices across all time points using both parametric
% and non-parametric tests, and groups consecutive significant points into clusters.
% 
% outputs:
% Y1 : vector of t-statistics at each time point (parametric test)
% Y2 : vector of p-values from parametric test (paired or unpaired t-test)
% Y3 : vector of Cohen's d effect sizes (parametric)
% Y4 : vector of z-statistics from non-parametric test
% Y5 : vector of p-values from non-parametric test (signrank or ranksum test)
% Y6 : structure describing positive clusters (parametric)
% Y7 : structure describing negative clusters (parametric)
% Y8 : structure describing positive clusters (non-parametric)
% Y9 : structure describing negative clusters (non-parametric)
% 
% inputs:
% X1 : data matrix for condition 1 (within-subject analysis) or group 1
%       (between-subject analysis) 
% X2 : data matrix for condition 2 (within-subject analysis) or group 2
%       (between-subject analysis) 
% X3 : string that determines whether to perform a dependent
%       (within-subject) or independent (between-subject) test 
% X4 : two-tailed or one-tailed test
% X5 : threshold for t-values to define clusters (parametric)
% X6 : threshold for z-values to define clusters (non-parametric)
% 
% Code created on September 10, 2024 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University


[numberOfPoints numberOfSubject1]=size(matrix1);
[numberOfPoints numberOfSubject2]=size(matrix2);

sizeEffectPARAMETRIC=zeros(1,numberOfPoints);

if dependentOrIndependent=="independent" %unpaired

    for cptTimepoint=1:numberOfPoints

            distributionAtChosenPoint1 = matrix1(cptTimepoint, :);
            distributionAtChosenPoint2 = matrix2(cptTimepoint, :);

            %We remove nan columns
            distributionAtChosenPoint1=distributionAtChosenPoint1(~isnan(distributionAtChosenPoint1));
            distributionAtChosenPoint2=distributionAtChosenPoint2(~isnan(distributionAtChosenPoint2));

            if ~ (isempty(distributionAtChosenPoint1) || isempty(distributionAtChosenPoint2))

                %PARAMETRIC
                [~, pValuePARAMETRIC(cptTimepoint), ~, statsPARAMETRIC] = ttest2(distributionAtChosenPoint1, distributionAtChosenPoint2, 'tail', tail);
                tValue(cptTimepoint) = statsPARAMETRIC.tstat;
                sizeEffectPARAMETRIC(cptTimepoint)=computeCohen_d(distributionAtChosenPoint1, distributionAtChosenPoint2);
                
                %NON PARAMETRIC
                [pValueNONPARAMETRIC(cptTimepoint), ~, statsNONPARAMETRIC]=ranksum(distributionAtChosenPoint1, distributionAtChosenPoint2, 'tail', tail);
                if isfield(statsNONPARAMETRIC, "zval")
                    zValue(cptTimepoint)=statsNONPARAMETRIC.zval;
                else 
                    zValue(cptTimepoint)=0;
                end
                
            else
                tValue(cptTimepoint)=0;
                zValue(cptTimepoint)=0;
                pValueNONPARAMETRIC(cptTimepoint)=1;    
            end
    end
    
elseif dependentOrIndependent=="dependent" %paired
 
    numberOfSubject=numberOfSubject1;
    
    distributionAtChosenPoint1=[];
    distributionAtChosenPoint2=[]; 
    
    for cptTimepoint=1:numberOfPoints

        distributionAtChosenPoint1 = matrix1(cptTimepoint, :);
        distributionAtChosenPoint2 = matrix2(cptTimepoint, :);

        %if ~ (isempty(distributionOfCondition1AtChosenPoint) || isempty(distributionOfCondition2AtChosenPoint))
         %We remove nan columns
        whoIsNanCondition1=isnan(distributionAtChosenPoint1);
        whoIsNanCondition2=isnan(distributionAtChosenPoint2);

        whoIsNan=or(whoIsNanCondition1, whoIsNanCondition2);
        distributionAtChosenPoint1=distributionAtChosenPoint1(~whoIsNan);
        distributionAtChosenPoint2=distributionAtChosenPoint2(~whoIsNan);

        if ~ (isempty(distributionAtChosenPoint1) || isempty(distributionAtChosenPoint2))
                
                if isequal( size(distributionAtChosenPoint1) , size(distributionAtChosenPoint2) )

                    %PARAMETRIC
                    [~, pValuePARAMETRIC(cptTimepoint), ~, statsPARAMETRIC] = ttest(distributionAtChosenPoint1, distributionAtChosenPoint2, 'tail', tail);
                    tValue(cptTimepoint) = statsPARAMETRIC.tstat;
                    sizeEffectPARAMETRIC(cptTimepoint)=computeCohen_d(distributionAtChosenPoint1, distributionAtChosenPoint2);

                    %NON PARAMETRIC
                    [pValueNONPARAMETRIC(cptTimepoint), ~, statsNONPARAMETRIC]=signrank(distributionAtChosenPoint1, distributionAtChosenPoint2, 'tail', tail);
                    if isfield(statsNONPARAMETRIC, "zval")
                        zValue(cptTimepoint)=statsNONPARAMETRIC.zval;
                    else 
                        zValue(cptTimepoint)=0;
                    end
                else
                    error(sprintf("Sizes not matching"))
                    
                end
                
        else
            tValue(cptTimepoint)=0;
            zValue(cptTimepoint)=0;
            pValuePARAMETRIC(cptTimepoint)=NaN;
            pValueNONPARAMETRIC(cptTimepoint)=NaN;
        end
        %end

    end
    
else
    error("Only valid arguments : ""dependent"" or ""indepedent"" ")
end

    
    %PARAMETRIC
    tValuePositive=zeros(size(tValue));
    tValueNegative=zeros(size(tValue));

    
        %We look we timepoints are significant
        tValuePositive(tValue>0)=tValue(tValue>0);
        tValueNegative(tValue<0)=tValue(tValue<0);

        initiallyFoundSignificantTimepointPositivePARAMETRIC=tValuePositive>=tThreshold;
        initiallyFoundSignificantTimepointNegativePARAMETRIC=tValueNegative<=-tThreshold;
        
        %We make clusters out of connected significant timepoints
        [initiallyFoundClusterPositivePARAMETRIC, numberOfPositiveClustersPARAMETRIC] = bwlabel(initiallyFoundSignificantTimepointPositivePARAMETRIC);
        [initiallyFoundClusterNegativePARAMETRIC, numberOfNegativeClustersPARAMETRIC] = bwlabel(initiallyFoundSignificantTimepointNegativePARAMETRIC);

        %We store sum and mean of t value, and lengths of clusters
        clusterSumOfTValuePositive = zeros(1, numberOfPositiveClustersPARAMETRIC); 
        clusterSumOfTValueNegative = zeros(1, numberOfNegativeClustersPARAMETRIC);

        clusterMeanOfTValuePositive = zeros(1, numberOfPositiveClustersPARAMETRIC); 
        clusterMeanOfTValueNegative = zeros(1, numberOfNegativeClustersPARAMETRIC);

        clusterLengthOfTValuePositive = zeros(1, numberOfPositiveClustersPARAMETRIC); 
        clusterLengthOfTValueNegative = zeros(1, numberOfNegativeClustersPARAMETRIC);

        clusterPositionPositivePARAMETRIC = {}; 
        clusterPositionNegativePARAMETRIC = {}; 
    
        
        for cptCluster = 1:numberOfPositiveClustersPARAMETRIC        
            [clusterLengthOfTValuePositive(cptCluster) clusterPositionPositivePARAMETRIC{cptCluster} ~] = findSizeOfClustersOfOnes(initiallyFoundClusterPositivePARAMETRIC == cptCluster); % Length of cluster cptCluster
            clusterMeanOfTValuePositive(cptCluster)   = mean(tValuePositive(initiallyFoundClusterPositivePARAMETRIC == cptCluster)); % Mean t-value of cluster cptCluter
            clusterSumOfTValuePositive(cptCluster)    = sum(tValuePositive(initiallyFoundClusterPositivePARAMETRIC == cptCluster)); % Sum of t-value of cluster cptCluter
            clusterMeanOfSizeEffectPositive(cptCluster)   = mean(sizeEffectPARAMETRIC(initiallyFoundClusterPositivePARAMETRIC == cptCluster)); % Mean Cohen's d value of cluster cptCluter
        end

        for cptCluster = 1:numberOfNegativeClustersPARAMETRIC        
            [clusterLengthOfTValueNegative(cptCluster) clusterPositionNegativePARAMETRIC{cptCluster} ~] = findSizeOfClustersOfOnes(initiallyFoundClusterNegativePARAMETRIC == cptCluster); % Length of cluster cptCluster
            clusterMeanOfTValueNegative(cptCluster)   = mean(tValueNegative(initiallyFoundClusterNegativePARAMETRIC == cptCluster)); % Mean t-value of cluster cptCluter
            clusterSumOfTValueNegative(cptCluster)    = sum(tValueNegative(initiallyFoundClusterNegativePARAMETRIC == cptCluster)); % Sum of t-value of cluster cptCluter
            clusterMeanOfSizeEffectNegative(cptCluster)   = mean(sizeEffectPARAMETRIC(initiallyFoundClusterNegativePARAMETRIC == cptCluster)); % Mean Cohen's d value of cluster cptCluter
        end

    
    
    %NON PARAMETRIC
    zValuePositive=zeros(size(zValue));
    zValueNegative=zeros(size(zValue));
    
        %We look we timepoints are significant
        zValuePositive(zValue>0)=zValue(zValue>0);
        zValueNegative(zValue<0)=zValue(zValue<0);

        initiallyFoundSignificantTimepointPositiveNONPARAMETRIC=zValuePositive>=zThreshold;
        initiallyFoundSignificantTimepointNegativeNONPARAMETRIC=zValueNegative<=-zThreshold;

         %We make clusters out of connected significant timepoints
        [initiallyFoundClusterPositiveNONPARAMETRIC, numberOfPositiveClustersNONPARAMETRIC] = bwlabel(initiallyFoundSignificantTimepointPositiveNONPARAMETRIC);
        [initiallyFoundClusterNegativeNONPARAMETRIC, numberOfNegativeClustersNONPARAMETRIC] = bwlabel(initiallyFoundSignificantTimepointNegativeNONPARAMETRIC);

        %We store sum and mean of z value, and lengths of clusters
        clusterSumOfZValuePositive = zeros(1, numberOfPositiveClustersNONPARAMETRIC); 
        clusterSumOfZValueNegative = zeros(1, numberOfNegativeClustersNONPARAMETRIC);

        clusterMeanOfZValuePositive = zeros(1, numberOfPositiveClustersNONPARAMETRIC); 
        clusterMeanOfZValueNegative = zeros(1, numberOfNegativeClustersNONPARAMETRIC);

        clusterLengthOfZValuePositive = zeros(1, numberOfPositiveClustersNONPARAMETRIC); 
        clusterLengthOfZValueNegative = zeros(1, numberOfNegativeClustersNONPARAMETRIC);

        clusterPositionPositiveNONPARAMETRIC = {}; 
        clusterPositionNegativeNONPARAMETRIC = {}; 

        for cptCluster = 1:numberOfPositiveClustersNONPARAMETRIC        
            [clusterLengthOfZValuePositive(cptCluster) clusterPositionPositiveNONPARAMETRIC{cptCluster} ~] = findSizeOfClustersOfOnes(initiallyFoundClusterPositiveNONPARAMETRIC == cptCluster); % Length of cluster cptCluster
            clusterMeanOfZValuePositive(cptCluster)   = mean(zValuePositive(initiallyFoundClusterPositiveNONPARAMETRIC == cptCluster)); % Mean t-value of cluster cptCluter
            clusterSumOfZValuePositive(cptCluster)    = sum(zValuePositive(initiallyFoundClusterPositiveNONPARAMETRIC == cptCluster)); % Sum of t-value of cluster cptCluter
        end

        for cptCluster = 1:numberOfNegativeClustersNONPARAMETRIC        
            [clusterLengthOfZValueNegative(cptCluster) clusterPositionNegativeNONPARAMETRIC{cptCluster} ~] = findSizeOfClustersOfOnes(initiallyFoundClusterNegativeNONPARAMETRIC == cptCluster); % Length of cluster cptCluster
            clusterMeanOfZValueNegative(cptCluster)   = mean(zValueNegative(initiallyFoundClusterNegativeNONPARAMETRIC == cptCluster)); % Mean t-value of cluster cptCluter
            clusterSumOfZValueNegative(cptCluster)    = sum(zValueNegative(initiallyFoundClusterNegativeNONPARAMETRIC == cptCluster)); % Sum of t-value of cluster cptCluter
        end
    
        
        
        
        
    
    if ~isempty(clusterLengthOfTValuePositive)
        positiveClusterPARAMETRIC.length=clusterLengthOfTValuePositive;
        positiveClusterPARAMETRIC.tValueMean=clusterMeanOfTValuePositive;
        positiveClusterPARAMETRIC.tValueSum=clusterSumOfTValuePositive;
        positiveClusterPARAMETRIC.sizeEffectMean=clusterMeanOfSizeEffectPositive;
        positiveClusterPARAMETRIC.position=clusterPositionPositivePARAMETRIC;
        positiveClusterPARAMETRIC.number=numberOfPositiveClustersPARAMETRIC;
    else
        positiveClusterPARAMETRIC.length=0;
        positiveClusterPARAMETRIC.tValueMean=0;
        positiveClusterPARAMETRIC.tValueSum=0;
        positiveClusterPARAMETRIC.sizeEffectMean=0;
        positiveClusterPARAMETRIC.position=[0 0];
        positiveClusterPARAMETRIC.number=0;
    end
    
    if ~isempty(clusterLengthOfTValueNegative)
        negativeClusterPARAMETRIC.length=clusterLengthOfTValueNegative;
        negativeClusterPARAMETRIC.tValueMean=clusterMeanOfTValueNegative;
        negativeClusterPARAMETRIC.tValueSum=clusterSumOfTValueNegative;
        negativeClusterPARAMETRIC.sizeEffectMean=clusterMeanOfSizeEffectNegative;
        negativeClusterPARAMETRIC.position=clusterPositionNegativePARAMETRIC;
        negativeClusterPARAMETRIC.number=numberOfNegativeClustersPARAMETRIC;
    else
        negativeClusterPARAMETRIC.length=0;
        negativeClusterPARAMETRIC.tValueMean=0;
        negativeClusterPARAMETRIC.tValueSum=0;
        negativeClusterPARAMETRIC.sizeEffectMean=0;
        negativeClusterPARAMETRIC.position=[0 0];
        negativeClusterPARAMETRIC.number=0;
    end

        
    
    
    if ~isempty(clusterLengthOfZValuePositive)
        positiveClusterNONPARAMETRIC.length=clusterLengthOfZValuePositive;
        positiveClusterNONPARAMETRIC.zValueMean=clusterMeanOfZValuePositive;
        positiveClusterNONPARAMETRIC.zValueSum=clusterSumOfZValuePositive;
        positiveClusterNONPARAMETRIC.sizeEffectMean=[];
        positiveClusterNONPARAMETRIC.position=clusterPositionPositiveNONPARAMETRIC;
        positiveClusterNONPARAMETRIC.number=numberOfPositiveClustersNONPARAMETRIC;
    else
        positiveClusterNONPARAMETRIC.length=0;
        positiveClusterNONPARAMETRIC.zValueMean=0;
        positiveClusterNONPARAMETRIC.zValueSum=0;
        positiveClusterNONPARAMETRIC.sizeEffectMean=[];
        positiveClusterNONPARAMETRIC.position=[0 0];
        positiveClusterNONPARAMETRIC.number=0;
    end
    
    if ~isempty(clusterLengthOfZValueNegative)
        negativeClusterNONPARAMETRIC.length=clusterLengthOfZValueNegative;
        negativeClusterNONPARAMETRIC.zValueMean=clusterMeanOfZValueNegative;
        negativeClusterNONPARAMETRIC.zValueSum=clusterSumOfZValueNegative;
        negativeClusterNONPARAMETRIC.sizeEffectMean=[];
        negativeClusterNONPARAMETRIC.position=clusterPositionNegativeNONPARAMETRIC;
        negativeClusterNONPARAMETRIC.number=numberOfNegativeClustersNONPARAMETRIC;
    else
        negativeClusterNONPARAMETRIC.length=0;
        negativeClusterNONPARAMETRIC.zValueMean=0;
        negativeClusterNONPARAMETRIC.zValueSum=0;
        negativeClusterNONPARAMETRIC.sizeEffectMean=[];
        negativeClusterNONPARAMETRIC.position=[0 0];
        negativeClusterNONPARAMETRIC.number=0;
    end


end