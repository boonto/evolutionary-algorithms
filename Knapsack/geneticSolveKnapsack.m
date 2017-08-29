function [best, nGenerations, maxFitnessOverGenerations, medianFitnessOverGenerations] = ...
geneticSolveKnapsack(geneticParameters)
    % GENETICSOLVEONEMAX genetically solve the one max problem.
    % Requires:
    % geneticParameters.nGenes Number of genes
    % geneticParameters.nIndividuals Number of initial individuals
    % geneticParameters.nMaxGenerations Maximum number of possible generations
    % geneticParameters.mutationRate Rate of mutation
    % geneticParameters.crossoverRate Rate of crossover
    % geneticParameters.targetGene The target gene which should be reached
    
    % Remap parameters
    items = geneticParameters.items;
    nGenes = geneticParameters.nGenes;
    nIndividuals = geneticParameters.nIndividuals;
    nMaxGenerations = geneticParameters.nMaxGenerations;
    mutationRate = geneticParameters.mutationRate;
    crossoverRate = geneticParameters.crossoverRate;
    maximumWeight = geneticParameters.maximumWeight;

    % Initialize output Parameters
    best = [];
    nGenerations = 0;
    maxFitnessOverGenerations = [];
    medianFitnessOverGenerations = [];
    
    % Initialize population randomly
    population = logical(round(rand(nIndividuals, nGenes)));
    
    % Counter for evolutions
    iGeneration = 0;
    
    % Evolution loop
    while 1
        % Evaluate population
        [maxFitness, idBest, medianFitness] = evaluatePopulation(population, items, maximumWeight);
        maxFitnessOverGenerations = cat(1, maxFitnessOverGenerations, maxFitness);
        medianFitnessOverGenerations = cat(1, medianFitnessOverGenerations, medianFitness);
        best = population(idBest, :);
        
        % Break loop when target is reached
        if (iGeneration == nMaxGenerations)
            break;
        end

        % Get best solution so far
        elite = population(idBest, :);
        
        % Selection
        parents = selection(population, items, maximumWeight);        

        % Crossover
        children = crossover(parents, crossoverRate);

        % Mutation
        mutants = mutation(children, mutationRate);

        % Replace worst mutant with elite
        [~, idWorst] = min(calcFitness(mutants, items, maximumWeight));
        mutants(idWorst, :) = elite;

        % New generation for next evolution cycle
        population = mutants;

        iGeneration = iGeneration + 1;
    end
    
    nGenerations = iGeneration;
end

%% Evaluate the population for fitness
function [maxFitness, id, medianFitness] = evaluatePopulation(population, items, maximumWeight)
    fitness = calcFitness(population, items, maximumWeight);
    [maxFitness, id] = max(fitness);
    medianFitness = median(fitness);
end

%% Fitness function
function fitness = calcFitness(population, items, maximumWeight)
    fitness = [];
    for iIndividual = 1:size(population, 1);
        weight = sum(items.weight(population(iIndividual, :)));
        if weight <= maximumWeight
            fitness = [fitness; sum(items.value(population(iIndividual, :)))];
        %elseif weight > (maximumWeight * 1.05)
        %    fitness = [fitness; sum(items.value(population(iIndividual, :))) * 0.75];
        else
            fitness = [fitness; 0];
        end
    end
end

%% Selection process
% Tournament style selection
function selected = selection(population, items, maximumWeight)
    nIndividuals = size(population, 1);

    % Generate matchup pairs randomly
    matchupsLeft = population(randperm(nIndividuals), :);
    matchupsRight = population(randperm(nIndividuals), :);
    
    % Calculate respective fitnesses
    fitnessLeft = calcFitness(matchupsLeft, items, maximumWeight);
    fitnessRight = calcFitness(matchupsRight, items, maximumWeight);
    
    % Designate the winners
    winnersLeft = matchupsLeft(fitnessLeft >= fitnessRight, :);
    winnersRight = matchupsRight(fitnessLeft < fitnessRight, :);
    
    selected = cat(1, winnersLeft, winnersRight);
end

%% Crossover process
% Combines parents based on crossoverRate
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
    children = cat(1, children, uncrossed);
end

%% Mutation process
% Mutate population based on the mutationRate
function mutants = mutation(population, mutationRate)
    nIndividuals = size(population, 1);
    nGenes = size(population, 2);
    
    % Generate multiple permutations 1:nGenes
    [~, out] = sort(rand(nIndividuals, nGenes), 2);
    % Pick genes to be mutated from permutations
    mutations = out <= nGenes*mutationRate;

    % Apply mutations
    mutants = xor(population, mutations);
end
