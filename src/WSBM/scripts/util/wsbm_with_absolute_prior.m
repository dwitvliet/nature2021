% making a separate wrapper function using wsbm using the grouping
% specified by labels
function [label, model] = wsbm_with_absolute_prior(graphIdx, label, outputDir, varargin)
    disp(numel(graphIdx));
    if numel(graphIdx) == 1
        E = load(sprintf('all_edge_list_%d.mat', graphIdx), 'arr');
        E = E.arr;
    elseif numel(graphIdx) > 1
        E = graphIdx;
    end
    if isempty(varargin)
        muIter = 100;
        mainIter = 100;
    else
        for ii = 1:2:length(varargin)
            if ischar(varargin{ii})
                switch lower(varargin{ii})
                    case "muIter"
                        muIter = varargin{ii};
                    case "mainIter"
                        mainIter = varargin{ii};
                    otherwise
                        muIter = 100;
                        mainIter = 100;
                end
            end
        end
    end
    disp("here");
    mu_0 = assignPrior(label);
    k = numel(unique(label));
    tic
    [label, model] = wsbm(E, k, 'W_Distr', 'Exponential', 'muMaxIter', muIter, 'mainMaxIter', mainIter, 'mu_0', mu_0);
    s = struct('model', model, 'label', label);
    save(outputDir,'s');
    toc
end

function mu_0 = assignPrior(label)
    n = size(label, 1);
    disp(n);
    mu_0 = zeros(max(unique(label)), n);
    disp(mu_0)
    mu_0(sub2ind(size(mu_0), label', linspace(1, n, n))) = 1;
end