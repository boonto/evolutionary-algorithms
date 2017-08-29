function [results, stats] = geneticSolveTravelingSalesman(params)
    % GENETICSOLVETRAVELINGSALESMAN genetically solve the traveling salesman problem.
    % Requires:
    % params.nGenes Number of genes
    % params.nIndividuals Number of initial individuals
    % params.nMaxGenerations Maximum number of possible generations
    % params.mutationMethod
    % params.mutationRate Rate of mutation
    % params.crossoverMethod
    % params.crossoverRate Rate of crossover
    
    % Remap parameters
    table = params.table;
    nGenes = params.nGenes;
    nIndividuals = params.nIndividuals;
    nMaxGenerations = params.nMaxGenerations;
    mutationMethod = params.mutationMethod;
    mutationRate = params.mutationRate;
    crossoverMethod = params.crossoverMethod;
    crossoverRate = params.crossoverRate;
    
    % Initialize output Parameters
    results.bestIndividual = [];
    stats.bestFitness = zeros(nMaxGenerations, 1);
    stats.medianFitness = zeros(nMaxGenerations, 1);
    
    % Initialize population randomly between 1 and nGenes
    [~, population] = sort(rand(nIndividuals, nGenes), 2);
    
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
        children = [];
        if (strcmp(crossoverMethod, 'OrderOne'))
            children = orderOneCrossover(parents, crossoverRate);
        elseif (strcmp(crossoverMethod, 'AlternatingPosition'))
            children = alternatingPositionCrossover(parents, crossoverRate);
        end

        % Mutation
        mutants = [];
        if (strcmp(mutationMethod, 'Swap'))
            mutants = swapMutation(children, mutationRate);
        elseif (strcmp(mutationMethod, 'Inversion'))
            mutants = inversionMutation(children, mutationRate);
        end

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
    [bestFitness, bestIndividual] = min(fitness);
    medianFitness = median(fitness);
end

%% Fitness function
function fitness = calcFitness(population, table)
    % Shift population by 1
    shiftedPopulation = circshift(population, -1, 2);
    distanceIndices = sub2ind(size(table), population, shiftedPopulation);
    fitness = sum(table(distanceIndices), 2);
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
% Combines parents by dropping a random gene sequence from parentA to
% parentB
% crossoverRate defines the amount of crossover happening
function children = orderOneCrossover(parents, crossoverRate)
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
       
    children = zeros(nParents, nGenes);
    for iParent = 1:nParents
        % Extract gene from parentA
        parentLeft = parentsLeft(iParent, :);
        parentRight = parentsRight(iParent, :);
        
        rangeBounds = sort(randi(nGenes, 1, 2));
        rangeLeft = 1:rangeBounds(1)-1;
        rangeMiddle = rangeBounds(1):rangeBounds(2);
        
        extractedGenesLeft = parentLeft(rangeMiddle);
        
        % Remove extracted genes from parentB
        child = parentRight(~ismember(parentRight,extractedGenesLeft));
        
        % Combine extracted gene with leftover genes from parentB
        child = [child(rangeLeft) extractedGenesLeft child(rangeBounds(1):size(child,2))];
        children(iParent, :) = child;
    end
    
    % Combine with uncrossed
    children = [children; uncrossed];
end

% Combines parents by alternatingly picking genes from them
% crossoverRate defines the amount of crossover happening
function children = alternatingPositionCrossover(parents, crossoverRate)
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
    
    % Pick alternating genes from parents
    combinedParents = zeros(nParents, nGenes*2);
    combinedParents(:, 1:2:end-1) = parentsLeft;
    combinedParents(:, 2:2:end) = parentsRight;
    
    children = zeros(size(parentsLeft));
    
    % Remove duplicates while preserving order
    for iChild = 1:nParents
        [~, ind] = unique(combinedParents(iChild,:));
        children(iChild, :) = combinedParents(iChild, sort(ind));
    end
    
    % Combine with uncrossed
    children = cat(1, children, uncrossed);
end

%% Mutation process
% Mutate population by swapping single genes
% mutationRate defines the amount of swapped genes
function mutants = swapMutation(population, mutationRate)
    nIndividuals = size(population, 1);
    nGenes = size(population, 2);
    
    % Correct mutation rate to achieve average mutationRate of genes
    % mutated; correct?
    mutationRate = mutationRate * 0.545;
    % Select genes which should be swapped
    genesToBeMutated = rand(nIndividuals, nGenes) < mutationRate;
    
    mutants = population;
    for iIndividual = 1:nIndividuals
        % Get indices of to be mutated genes
        ix1 = find(genesToBeMutated(iIndividual, :) == 1);
        % Get swapping partner
        ix2 = ix1 + 1;
        for iIndex = 1:size(ix1, 2)
            if ix1(iIndex) == nGenes
                ix2(iIndex) = 1;
            end
            % Swap chosen gene
            mutants(iIndividual, [ix1(iIndex) ix2(iIndex)]) = mutants(iIndividual, [ix2(iIndex) ix1(iIndex)]);
        end
    end
    
    % Used to test average mutation rate
%     mean(sum(mutants ~= population, 2) / nGenes) 
end

% Mutate population by inversing gene sequences
% mutationRate defines the length of the inversed sequence
function mutants = inversionMutation(population, mutationRate)
    nIndividuals = size(population, 1);
    nGenes = size(population, 2);

    % This mutation method mutates a gene by 1/3 on average, amplify
    % mutationRate; correct?
    mutationRate = 3 * mutationRate;
    
    % Shuffle and apply crossoverRate
    population = population(randperm(nIndividuals), :);
    nMutants = round(mutationRate * nIndividuals);
    unmutated = population(nMutants+1:nIndividuals, :);
    mutants = population(1:nMutants, :);
    
    for iMutant = 1:nMutants
        % Select two random inversion points
        inversionPoints = sort(randi(nGenes,1,2));
        range = inversionPoints(1):inversionPoints(2);
        % Inverse range
        individual = population(iMutant, :);
        mutants(iMutant, range) = fliplr(individual(range));
    end
    
    % Combine with unmutated
    mutants = cat(1, mutants, unmutated);
    
    % Used to test average mutation rate
%     mean(sum(mutants ~= population, 2) / nGenes) 
end
