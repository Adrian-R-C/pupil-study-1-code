function [clusterLengths clusterPositions numberOfClusterFound] = findSizeOfClustersOfOnes(vectorOfOnesAndZeros)
% findSizeOfClustersOfOnes finds clusters of 1 and returns their length, onset and offset, and number of clusters  
%
% [Y1, Y2, Y3] = findSizeOfClustersOfOnes(X) 
% 
% outputs:
% Y1 : clusters length
% Y2 : onset and offset of the clusters
% Y3 : number of clusters
% 
% input:
% X : binary vector
% 
% Code created on August 17, 2024 by
% Adrian RUIZ CHIAPELLO
% Centre de Recherche Cerveau et Cognition
% CNRS / Toulouse University



derivativeVector = diff([0 vectorOfOnesAndZeros 0]); % On calcule les différences
startIndices = find(derivativeVector == 1); % On regarde quand les clusters débutent
endIndices = find(derivativeVector == -1) - 1; % On regarde quand les clusters se terminent

clusterLengths = endIndices - startIndices + 1;


if ~isempty(clusterLengths)
    
    for cpt=1:length(clusterLengths)
        clusterPositions(cpt,1)=startIndices(cpt);
        clusterPositions(cpt,2)=endIndices(cpt);
    end
    
    numberOfClusterFound=length(clusterLengths);

else
    clusterPositions=[];
    clusterLengths=0;
    numberOfClusterFound=0;
end

end