function [results, stats] = geneticSolveGraphColoring(params)
    % GENETICSOLVEGRAPHCOLORING genetically solve the graph coloring problem.
    % Requires:
    % params.table Graph that should be colored
    % params.nColors Number of colors used
    % params.nGenes Number of genes
    % params.nIndividuals Number of initial individuals
    % params.nMaxGenerations Maximum number of possible generations
    % params.mutationRate Rate of mutation
    % params.crossoverRate Rate of crossover
    
    % Remap parameters
    table = params.table;
    nColors = params.nColors;
    nGenes = params.nGenes;
    nIndividuals = params.nIndividuals;
    nMaxGenerations = params.nMaxGenerations;
    mutationRate = params.mutationRate;
    crossoverRate = params.crossoverRate;
    
    % Initialize output Parameters
    results.bestIndividual = [];
    stats.bestFitness = zeros(nMaxGenerations, 1);
    stats.medianFitness = zeros(nMaxGenerations, 1);
    
    % Initialize population randomly between 1 and nGenes
    population = randi(nColors, nIndividuals, nGenes);
    
    % Evolution loop
    for iGeneration = 1:nMaxGenerations
        % Evaluate population
        [bestIndividualId, ...
            stats.bestFitness(iGeneration), ...
            stats.medianFitness(iGeneration)] = evaluatePopulation(population, table);
        
        % Elite
        results.bestIndividual = population(bestIndividualId, :);
        
        % Selection
        parents = selection(population, table);        

        % Crossover
        children = crossover(parents, crossoverRate);

        % Mutation
        mutants = mutation(children, mutationRate, nColors);

        % Replace worst mutant with elite
        [~, worstIndividual] = max(calcFitness(mutants, table));
        mutants(worstIndividual, :) = results.bestIndividual;

        % New generation for next evolution cycle
        population = mutants;
    end
end

%% Evaluate the population for fitness
function [bestIndividual, bestFitness, medianFitness] = evaluatePopulation(population, table)
    fitness = calcFitness(population, table);
    [bestFitness, bestIndividual] = min(fitness); % minimization problem
    medianFitness = median(fitness);
end

%% Fitness function
function fitness = calcFitness(population, table)
    % Shift population by 1
    nIndividuals = size(population, 1);
    nGenes = size(population, 2);
    
    % Calculate fitness for each individual 
    penalty = [];
    for iIndividual = 1:nIndividuals
        indPenalty = 0;
        for node=1:nGenes
            for neigh=find(table(node,:)==1)
                if population(iIndividual, node) == population(iIndividual, neigh)
                    indPenalty = indPenalty + 1;
                end
            end
        end
        penalty(iIndividual) = indPenalty;
    end

    fitness = penalty';
end

%% Selection process
% Tournament style selection
function selected = selection(population, table)
    nIndividuals = size(population, 1);

    % Generate matchup pairs randomly
    matchupsLeft = population(randperm(nIndividuals), :);
    matchupsRight = population(randperm(nIndividuals), :);
    
    % Calculate respective fitnesses
    fitnessLeft = calcFitness(matchupsLeft, table);
    fitnessRight = calcFitness(matchupsRight, table);
    
    % Designate the winners
    winnersLeft = matchupsLeft(fitnessLeft <= fitnessRight, :);
    winnersRight = matchupsRight(fitnessLeft > fitnessRight, :);
    
    selected = cat(1, winnersLeft, winnersRight);
end

%% Crossover process
% Crossover parents by simply combining them in the middle
function children = crossover(parents, crossoverRate)
    nIndividuals = size(parents, 1);
    nGenes = size(parents, 2);
    
    % Generate parent pairs randomly
    parentsLeft = parents(randperm(nIndividuals), :);
    parentsRight = parents(randperm(nIndividuals), :);
    
    % Apply crossoverRate
    nParents = round(crossoverRate * nIndividuals);
    uncrossed = parentsLeft(nParents+1:nIndividuals, :);
    parentsLeft = parentsLeft(1:nParents, :);
    parentsRight = parentsRight(1:nParents, :);
       
    % Extract genes from parents
    extractedGenesLeft = parentsLeft(:, 1:round(0.5 * nGenes));
    extractedGenesRight = parentsRight(:, (round(0.5 * nGenes)+1):nGenes);
    
    % Combine genes
    children = cat(2, extractedGenesLeft, extractedGenesRight);
    
    % Combine with uncrossed
    children = [children; uncrossed];
end

%% Mutation process
% Mutate population by incrementing genes
function mutants = mutation(population, mutationRate, nColors)
    nIndividuals = size(population, 1);
    nGenes = size(population, 2);
    
    % Select genes which should be swapped
    genesToBeMutated = rand(nIndividuals, nGenes) < mutationRate;
    
    mutants = [];
    for iIndividual = 1:nIndividuals
        mutant = population(iIndividual, :);
        for iGene = 1:nGenes
            if genesToBeMutated(iIndividual, iGene)
                % Increment gene in the specified range
                mutant(iGene) = mutant(iGene)+1;
                if mutant(iGene) > nColors
                    mutant(iGene) = 1;
                end
            end
        end
        mutants = [mutants; mutant];
    end
end
