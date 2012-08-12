function result = classifyNeuralEvents(timestamps, timelag, neural, cutoff)
    neural = zscore(neural');
    [~, score, ~, ~] = princomp(neural);
    neural = score';
    result = struct('Z', 0, 'C', 0);
    dist_matrix = zeros(1, length(timestamps)*(length(timestamps)-1)/2);
    for neuron = 1 : size(neural, 1)
        all_neural = zeros(length(timestamps), 2*timelag + 1);
        for idx = 1 : length(timestamps)
            all_neural(idx, :) = neural(neuron, timestamps(idx)-timelag: timestamps(idx)+timelag);
        end
        dist_matrix = dist_matrix + pdist(all_neural, 'correlation').^2;
    end
    dist_matrix = dist_matrix.^0.5;
    
    Z = linkage(dist_matrix, 'ward');
    C = cluster(Z, 'cutoff', cutoff, 'criterion', 'distance');
    result.Z = Z;
    result.C = C;
    clust_size = arrayfun(@(ind) length(find(C == ind)), 1:max(C));
    h = bar(1:max(C), clust_size);
    h_ch = get(h, 'Children');
    set(h_ch, 'CData', 1:max(C));
    title('Set of an event cluster');
end

