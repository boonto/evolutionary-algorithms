%% Traveling Salesman problem
clear variables;
close all;

% Load dataset 
cities = table2cell(readtable('cities.csv'));
% Precalculate distances
distances = squareform(pdist([cell2mat(cities(:,3)) cell2mat(cities(:,4))]));

% Standard Usage
if 0
    params.table = distances;
    params.nGenes = size(params.table, 1);
    params.nIndividuals = 20;
    params.nMaxGenerations = 5000;
    params.crossoverMethod = 'OrderOne';
%     params.crossoverMethod = 'AlternatingPosition';
    params.crossoverRate = 0.8;
%     params.mutationMethod = 'Swap';
    params.mutationMethod = 'Inversion';
%     params.mutationRate = 1/params.nGenes;
    params.mutationRate = 0.1;

    [results, stats] = geneticSolveTravelingSalesman(params);

    % Plot the results
    figure(1);
    ax1 = axes('Position', [0 0 1 1], 'Visible', 'off');
    ax2 = axes('Position', [0.35 0.1 0.6 0.85]);
    lines = plot(0:params.nMaxGenerations-1, stats.bestFitness, ...
         0:params.nMaxGenerations-1, stats.medianFitness);
%     ax3 = gca;
    lines(1).LineWidth = 2.0;
    lines(2).LineWidth = 1.0;
    axis([0 5000, 0 200]);
    grid on;
    
    % Label everything
    title('Traveling Salesman Problem');
    xlabel('Generation');
    ylabel('Distance');
    legend('min', 'median', 'Location', 'northeast');

    % Additional information
    axes(ax1);
    description = {'Parameters:';
                   ['nGenes=' num2str(params.nGenes)];
                   ['nIndividuals=' num2str(params.nIndividuals)];
                   ['crossoverMethod=' params.crossoverMethod];
                   ['crossoverRate=' num2str(params.crossoverRate)];
                   ['mutationMethod=' params.mutationMethod];
                   ['mutationRate=' num2str(params.mutationRate)];};
    text(0.02, 0.5, description);

    figure(2);
    geoshow([cell2mat(cities(results.bestIndividual(1,:), 3)); ...
             cell2mat(cities(results.bestIndividual(1,1), 3))], ...
            [cell2mat(cities(results.bestIndividual(1,:), 4)); ...
             cell2mat(cities(results.bestIndividual(1,1), 4))]);
    
    min(stats.bestFitness)
end

% Examine effects of different crossover and mutation methods
if 0
    params.table = distances;
    params.nGenes = size(params.table, 1);
    params.nIndividuals = 20;
    params.nMaxGenerations = 5000;
    params.crossoverRate = 0.8;
    params.mutationRate = 0.1;
    
    nSamples = 100;
    
    params.crossoverMethod = 'OrderOne';
    params.mutationMethod = 'Swap';
    parfor iSample = 1:nSamples
        [results{1, iSample}, stats{1, iSample}] = geneticSolveTravelingSalesman(params);
    end
    
    params.crossoverMethod = 'OrderOne';
    params.mutationMethod = 'Inversion';
    parfor iSample = 1:nSamples
        [results{2, iSample}, stats{2, iSample}] = geneticSolveTravelingSalesman(params);
    end
    
    params.crossoverMethod = 'AlternatingPosition';
    params.mutationMethod = 'Swap';
    parfor iSample = 1:nSamples
        [results{3, iSample}, stats{3, iSample}] = geneticSolveTravelingSalesman(params);
    end
    
    params.crossoverMethod = 'AlternatingPosition';
    params.mutationMethod = 'Inversion';
    parfor iSample = 1:nSamples
        [results{4, iSample}, stats{4, iSample}] = geneticSolveTravelingSalesman(params);
    end
    
    save('TSP_data', 'results', 'stats');
end

if 1
    % post process data
    tspData = load('TSP_data.mat');
    a = cell2mat(tspData.stats);
    
    bestFitness = zeros(size([a.bestFitness], 1) * size(a, 2), size(a, 1));
    bestFitness = [];
    for iConfiguration = 1:size(a, 1)
        b = min([a(iConfiguration, :).bestFitness]);
        bestFitness = [bestFitness, b(:)];
    end
    
    % visualize data
    figure(3);
    ax1 = axes('Position', [0 0 1 1], 'Visible', 'off');
    ax2 = axes('Position', [0.35 0.1 0.6 0.85]);
    
    yEntries = {'OO & S', 'OO & I', 'AP & S', 'AP & I'};    
   
    boxplot(bestFitness);
    set(gca, 'XTick', 1:4, 'XTickLabel', yEntries);
    title('Examination of different crossover & mutation methods');
    xlabel('Configuration');
    ylabel('Distance');
    grid on;
    
    axes(ax1);
    description = {'Crossover methods:';
                   ' Order One(OO)';
                   ' Alternating Position(AP)';
                   '';
                   'Mutation methods:';
                   ' Swap(S)';
                   ' Inversion(I)';
                   '';
                   'Samples = 100'};
    text(0.02, 0.5, description);
end