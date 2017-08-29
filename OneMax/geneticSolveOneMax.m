function [best, nGenerations, maxFitnessOverGenerations, meanFitnessOverGenerations] = ...
geneticSolveOneMax(geneticParameters)
    % GENETICSOLVEONEMAX genetically solve the one max problem.
    % Requires:
    % geneticParameters.nGenes Number of genes
    % geneticParameters.nIndividuals Number of initial individuals
    % geneticParameters.nMaxGenerations Maximum number of possible generations
    % geneticParameters.mutationRate Rate of mutation
    % geneticParameters.crossoverRate Rate of crossover
    % geneticParameters.targetGene The target gene which should be reached
    
    % Remap parameters
    nGenes = geneticParameters.nGenes;
    nIndividuals = geneticParameters.nIndividuals;
    nMaxGenerations = geneticParameters.nMaxGenerations;
    mutationRate = geneticParameters.mutationRate;
    crossoverRate = geneticParameters.crossoverRate;
    targetGene = geneticParameters.targetGene;
    targetFitness = calcFitness(targetGene, targetGene);

    % Initialize output Parameters
    best = [];
    nGenerations = 0;
    maxFitnessOverGenerations = [];
    meanFitnessOverGenerations = [];
    
    % Initialize population randomly
    population = logical(round(rand(nIndividuals, nGenes)));
    
    % Counter for evolutions
    iGeneration = 0;
    
    % Evolution loop
    while 1
        % Evaluate population
        [maxFitness, idBest, meanFitness] = evaluatePopulation(population, targetGene);
        maxFitnessOverGenerations = cat(1, maxFitnessOverGenerations, maxFitness);
        meanFitnessOverGenerations = cat(1, meanFitnessOverGenerations, meanFitness);
        best = population(idBest, :);
        
        % Break loop when target is reached
        if (maxFitness == targetFitness || iGeneration == nMaxGenerations)
            break;
        end

        % Get best solution so far
        elite = population(idBest, :);
        
        % Selection
        parents = selection(population, targetGene);        

        % Crossover
        children = crossover(parents, crossoverRate);

        % Mutation
        mutants = mutation(children, mutationRate);

        % Replace worst mutant with elite
        [~, idWorst] = min(calcFitness(mutants, targetGene));
        mutants(idWorst, :) = elite;

        % New generation for next evolution cycle
        population = mutants;

        iGeneration = iGeneration + 1;
    end
    
    nGenerations = iGeneration;
end

%% Evaluate the population for fitness
function [maxFitness, id, meanFitness] = evaluatePopulation(population, targetGene)
    fitness = calcFitness(population, targetGene);
    [maxFitness, id] = max(fitness);
    meanFitness = mean(fitness);
end

%% Fitness function
function fitness = calcFitness(population, targetGene)
    %fitness = sum(population, 2); %wrong fitness function for targets != all ones
    
    % modified fitness function to support different targetGenes
    mask = bsxfun(@eq, population, targetGene);
    fitness = sum(mask, 2);
end

%% Selection process
% Tournament style selection
function selected = selection(population, targetGene)
    nIndividuals = size(population, 1);

    % Generate matchup pairs randomly
    matchupsLeft = population(randperm(nIndividuals), :);
    matchupsRight = population(randperm(nIndividuals), :);
    
    % Calculate respective fitnesses
    fitnessLeft = calcFitness(matchupsLeft, targetGene);
    fitnessRight = calcFitness(matchupsRight, targetGene);
    
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
