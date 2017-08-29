%% NSGA-II
% clear variables
close all

params.nIndividuals = 20;
params.nGenes = 16;
params.nMaxGenerations = 60;
params.mutationRate = 0.1;
params.visualize = 1;
params.writeFrames = 0;

[results, stats] = nsga_ii(params);

plot([stats.bestFitness]);
hold on;
grid on;
plot([stats.medianFitness]);
title('NSGA-II');
xlabel('Generation');
ylabel('Fitness');
legend('Best', 'Median');

