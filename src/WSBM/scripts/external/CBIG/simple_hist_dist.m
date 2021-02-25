function [emd_dist , ks_stat] = simple_hist_dist(x1,x2)
% calculate the earth mover's distance as the area between the csfs
% 1-dim, with all obs having equal weight

% from the 'evaluate_generative_model code...
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

emd_dist  =  sum(abs(sampleCDF1 - sampleCDF2));
ks_stat = max(abs(sampleCDF1 - sampleCDF2));



