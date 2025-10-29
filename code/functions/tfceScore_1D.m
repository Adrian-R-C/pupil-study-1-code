function tfceScore = tfceScore_1D(data, H, E, dH)

    tfceScore = zeros(size(data));
    h_values = 0:dH:max(abs(data));
    
    for h = h_values
        supra = abs(data) > h;
        clusters = bwconncomp(supra); 
        for c = 1:clusters.NumObjects
            idx = clusters.PixelIdxList{c};
            clusterSize = length(idx);
            contribution = (clusterSize^E) * (h^H) * dH;
            tfceScore(idx) = tfceScore(idx) + contribution;
        end
    end
    
    % We keep the sign
    tfceScore = tfceScore .* sign(data);
end
