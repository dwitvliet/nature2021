function [pairs, nonpairs] = countPairs(res)
    pair_info = csvread("node_pairs.csv");
    res = table2array(res(:, [2 3]));
    pairs = 0;
    nonpairs = 0;
    for i = 1: size(pair_info, 1)
        neuron1 = pair_info(i, 1);
        neuron2 = pair_info(i, 2);

        if ismember(neuron1, res(:, 1)) && ismember(neuron2, res(:, 1))
            if res(res(:, 1) == neuron1, 2) == res(res(:, 1) == neuron2, 2)
                pairs = pairs + 1;
            else
                nonpairs = nonpairs + 1;
            end
        end
    end
end