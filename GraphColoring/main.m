%% Graph coloring problem
clear variables;
close all;

% adjacency = [   0,1,0,1,0;
%                 0,0,1,0,0;
%                 0,0,0,1,1;
%                 0,0,0,0,1;
%                 0,0,0,0,0];
adjacency = bucky;

if 1
    params.table = adjacency;
    params.nColors = 4;
    params.nGenes = size(params.table, 1);
    params.nIndividuals = 10;
    params.nMaxGenerations = 200;
    params.crossoverRate = 1.0;
    params.mutationRate = 1/params.nGenes;

    [results, stats] = geneticSolveGraphColoring(params);

    % Plot the results
    figure(1);
    plot(0:params.nMaxGenerations-1, stats.bestFitness, ...
    0:params.nMaxGenerations-1, stats.medianFitness);
    axis([0, params.nMaxGenerations, -5, params.nGenes]);
    grid on;

    colors = jet(params.nColors);
    g = graph(adjacency,'upper');

    fig = figure(2);
    h = plot(g,'MarkerSize',12,'LineWidth',3);
    for i=1:params.nColors;highlight(h,find(results.bestIndividual==i),'NodeColor',colors(i,:));end
    title('Colored bucky graph');
    ax = gca;
    ax.XTick = []; ax.YTick = [];

    min(stats.bestFitness)
end

if 0
    params.table = adjacency;
    params.nColors = 4;
    params.nGenes = size(params.table, 1);
    params.nIndividuals = 10;
    params.nMaxGenerations = 200;
    params.crossoverRate = 1;
    params.mutationRate = 1/params.nGenes;
    
    nSamples = 100;
    
    parfor iSample = 1:nSamples
        [results{1, iSample}, stats{1, iSample}] = geneticSolveGraphColoring(params);
    end
    
    save('GC_data', 'results', 'stats');
end

if 0
    % post process data
    tspData = load('GC_data.mat');
    a = cell2mat(tspData.stats);
    
    bestFitness = [];
    for iConfiguration = 1:size(a, 1)
        b = [a(iConfiguration, :).bestFitness];
        bestFitness = [bestFitness, b(:)];
    end
    
    % visualize data
    fig = figure(3);
    boxplot(bestFitness);
    title('Performance over 100 samples');
    ylabel('Penalty');
    grid on;
    ax = gca;
    ax.XTick =[];
    save_figure(fig, ['boxplot'], 18); 
end
