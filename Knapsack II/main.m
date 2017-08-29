%% NSGA-II
% clear variables
close all

params.items = readtable('items.csv');
params.nIndividuals = 30;
params.nGenes = size(params.items, 1);
params.nMaxGenerations = 200;
params.mutationRate = 0.1;
params.visualize = 1;
params.writeFrames = 0;

[results] = nsga_ii(params);

exampleInd = results.elite(1:10);
fitness = reshape([exampleInd.fitness], 2, 10)';
value = fitness(:,1);
[~, ind] = sort(value);
exampleInd = exampleInd(ind);
fid = fopen('output.txt', 'w');
for i = 1:size(exampleInd, 2)
    items = params.items.item(exampleInd(i).genes)';
    fprintf(fid, '%s ', items{:});
    fprintf(fid, '\n');
end
fclose(fid);