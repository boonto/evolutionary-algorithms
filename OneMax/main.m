%% One Max problem
clear variables
close all

% Standard usage
if 1
    geneticParameters.nGenes = 16;
    geneticParameters.nIndividuals = 10;
    geneticParameters.nMaxGenerations = 1000;
    geneticParameters.mutationRate = 1/geneticParameters.nGenes;
    geneticParameters.crossoverRate = 1.0;
    geneticParameters.targetGene = true(1, geneticParameters.nGenes);
    
    % Solve the problem
    [best, nGenerations, maxFitnessOverGenerations, meanFitnessOverGenerations] = geneticSolveOneMax(geneticParameters);

    % Plot the results
    figure(1);
    ax1 = axes('Position', [0 0 1 1], 'Visible', 'off');
    ax2 = axes('Position', [0.35 0.1 0.6 0.85]);
    plot(0:nGenerations, maxFitnessOverGenerations, 0:nGenerations, meanFitnessOverGenerations);

    % Label everything
    title('One Max Problem');
    xlabel('Generation');
    ylabel('Fitness');
    legend('max', 'mean', 'Location', 'southeast');

    % Additional information
    axes(ax1);
    description = {'Parameters:';
                   ['nGenes=' num2str(geneticParameters.nGenes)];
                   ['nIndividuals=' num2str(geneticParameters.nIndividuals)];
                   ['crossoverRate=' num2str(geneticParameters.crossoverRate)];
                   ['mutationRate=' num2str(geneticParameters.mutationRate)];
                   '';
                   'Reached in';
                   ['Generation ' num2str(nGenerations)]};
    text(0.02, 0.5, description);
end


% Examine effect of mutation rate
if 0
    geneticParameters.nGenes = 16;
    geneticParameters.nIndividuals = 10;
    geneticParameters.nMaxGenerations = 1000;
    geneticParameters.crossoverRate = 1.0;
    geneticParameters.targetGene = true(1, geneticParameters.nGenes);
    
    meanGenerations = [];
    mutationRate = 0:0.05:1;
    for iRate = 1:length(mutationRate)
        geneticParameters.mutationRate = mutationRate(iRate);
        nGensPerSample = [];
        for iSample = 1:100    
            [~, nGenerations, ~, ~] = geneticSolveOneMax(geneticParameters);
            nGensPerSample = cat(1, nGensPerSample, nGenerations);
        end
        meanGenerations = cat(1, meanGenerations, mean(nGensPerSample));
    end
    figure(2);
    plot(mutationRate, meanGenerations);
    title('Effect of mutation rate');
    xlabel('mutationRate');
    ylabel('mean nGenerations');
end

% Examine effect of crossover rate
if 0
    geneticParameters.nGenes = 16;
    geneticParameters.nIndividuals = 10;
    geneticParameters.nMaxGenerations = 1000;
    geneticParameters.mutationRate = 1/geneticParameters.nGenes;
    geneticParameters.targetGene = true(1, geneticParameters.nGenes);
    
    meanGenerations = [];
    crossoverRate = 0:0.05:1;
    for iRate = 1:length(crossoverRate)
        geneticParameters.crossoverRate = crossoverRate(iRate);
        nGensPerSample = [];
        for iSample = 1:10000    
            [~, nGenerations, ~, ~] = geneticSolveOneMax(geneticParameters);
            nGensPerSample = cat(1, nGensPerSample, nGenerations);
        end
        meanGenerations = cat(1, meanGenerations, mean(nGensPerSample));
    end
    figure(3);
    plot(crossoverRate, meanGenerations);
    title('Effect of crossover rate');
    xlabel('crossoverRate');
    ylabel('mean nGenerations');
end

% Examine different target genes
if 0
    geneticParameters.nGenes = 8;
    geneticParameters.nIndividuals = 10;
    geneticParameters.nMaxGenerations = 1000;
    geneticParameters.mutationRate = 1/geneticParameters.nGenes;
    geneticParameters.crossoverRate = 1.0;
    
    nSamples = 1000;
    meanGenerations = [];
    
    geneticParameters.targetGene = [1 1 1 1 1 1 1 1];
    nGensPerSample = [];
    for iSample = 1:nSamples    
        [~, nGenerations, ~, ~] = geneticSolveOneMax(geneticParameters);
        nGensPerSample = cat(1, nGensPerSample, nGenerations);
    end
    meanGenerations = cat(1, meanGenerations, mean(nGensPerSample));
    
    geneticParameters.targetGene = [1 0 0 0 0 0 0 0];
    nGensPerSample = [];
    for iSample = 1:nSamples    
        [~, nGenerations, ~, ~] = geneticSolveOneMax(geneticParameters);
        nGensPerSample = cat(1, nGensPerSample, nGenerations);
    end
    meanGenerations = cat(1, meanGenerations, mean(nGensPerSample));
    
    geneticParameters.targetGene = [1 1 1 1 0 0 0 0];
    nGensPerSample = [];
    for iSample = 1:nSamples    
        [~, nGenerations, ~, ~] = geneticSolveOneMax(geneticParameters);
        nGensPerSample = cat(1, nGensPerSample, nGenerations);
    end
    meanGenerations = cat(1, meanGenerations, mean(nGensPerSample));
    
    geneticParameters.targetGene = [1 0 1 0 1 0 1 0];
    nGensPerSample = [];
    for iSample = 1:nSamples    
        [~, nGenerations, ~, ~] = geneticSolveOneMax(geneticParameters);
        nGensPerSample = cat(1, nGensPerSample, nGenerations);
    end
    meanGenerations = cat(1, meanGenerations, mean(nGensPerSample));
    
    figure(4);
    yEntries = {'[11111111]', '[10000000]', '[11110000]', '[10101010]'};    
   
    bar(meanGenerations);
    set(gca, 'XTick', 1:4, 'XTickLabel', yEntries);
    title('Different target genes');
    xlabel('targetGene');
    ylabel('mean nGenerations');
end