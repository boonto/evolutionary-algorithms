%% Knapsack problem
clear variables
close all

geneticParameters.items = readtable('items.csv');
%geneticParameters.items = geneticParameters.items([1:5], :);
geneticParameters.nGenes = size(geneticParameters.items, 1);
geneticParameters.nIndividuals = 10;
geneticParameters.nMaxGenerations = 400;
geneticParameters.mutationRate = 1/geneticParameters.nGenes;
geneticParameters.crossoverRate = 1.0;
geneticParameters.maximumWeight = 400;

% Solve the problem
[best, nGenerations, maxFitnessOverGenerations, medianFitnessOverGenerations] = geneticSolveKnapsack(geneticParameters);

maxFitnessOverGenerations(nGenerations)
% Plot the results
figure(1);
ax1 = axes('Position', [0 0 1 1], 'Visible', 'off');
ax2 = axes('Position', [0.35 0.1 0.6 0.85]);
plot(0:nGenerations, maxFitnessOverGenerations, 0:nGenerations, medianFitnessOverGenerations);

% Label everything
title('Knapsack Problem');
xlabel('Generations');
ylabel('Fitness');
legend('max', 'median', 'Location', 'southeast');

% Additional information
axes(ax1);
description = {'Parameters:';
               ['nGenes=' num2str(geneticParameters.nGenes)];
               ['nIndividuals=' num2str(geneticParameters.nIndividuals)];
               ['crossoverRate=' num2str(geneticParameters.crossoverRate)];
               ['mutationRate=' num2str(geneticParameters.mutationRate)];
               ['nMaxGenerations=' num2str(geneticParameters.nMaxGenerations)];
               ['maximumWeight=' num2str(geneticParameters.maximumWeight)];};
text(0.02, 0.5, description);

output = cat(2, geneticParameters.items.item(best), num2cell(geneticParameters.items.weight(best)), num2cell(geneticParameters.items.value(best)));
output = cat(1, output, ['sum', num2cell(sum(geneticParameters.items.weight(best))), num2cell(sum(geneticParameters.items.value(best)))])
display(output);

