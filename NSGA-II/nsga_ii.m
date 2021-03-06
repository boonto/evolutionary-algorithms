%% NSGA-II
function [ results, stats ] = nsga_ii( params )

% Remap parameters
nIndividuals = params.nIndividuals;
nGenes = params.nGenes;
nMaxGenerations = params.nMaxGenerations;
mutationRate = params.mutationRate;
visualize = params.visualize;
writeFrames = params.writeFrames;

% Define objective functions
objectives = {@leadingZeros, @trailingOnes};

%% Initialization
% Create unique individuals
population = [];
% TODO ugly hack
while size(population, 2) ~= nIndividuals
    for iFront = 1:nIndividuals
        population(iFront).genes = logical(round(rand(1, nGenes)));
        population(iFront).S = [];
        population(iFront).n = 0;
        population(iFront).rank = 0;
        population(iFront).dist = 0;
        for m = 1:size(objectives, 2)
            population(iFront).fitness(m) = 0;
        end
    end
    population = structUnion(population, population);
end

%% Generation loop
iGeneration = 0;
while iGeneration < nMaxGenerations
    iGeneration = iGeneration + 1;
    
    %% Offspring
    offspring = population;
    
    % Mutation
    for iMutant = 1:size(population, 2)
        offspring(iMutant).genes = xor(offspring(iMutant).genes, (rand(1, nGenes) < mutationRate)); 
    end
    
    %% Choose new generation
    combinedPop = structUnion(population, offspring);
    
    % Calculate Fitness
    combinedPop = fitnessAssignment(combinedPop, objectives);
    
    % Record stats
    fitness = [];
    for iInd = 1:size(combinedPop, 2)
        fitness = [fitness, sum(combinedPop(iInd).fitness)];
    end
    stats(iGeneration).medianFitness = median(fitness);
    stats(iGeneration).bestFitness = max(fitness);
    
    % Generate fronts
    fronts = dominationSort(combinedPop);
    
    % Start filling next population
    nextPopulation = [];
    iFront = 1;
    while size(nextPopulation, 2) + size(fronts{iFront}, 2) <= nIndividuals
        nextPopulation = structUnion(nextPopulation, fronts{iFront});
        iFront = iFront + 1;
    end
    
    % Get diverse individuals from leftover front
    fronts{iFront} = crowdingDistanceAssignment(fronts{iFront}, objectives);
    % Sort by crowding
    for n = size(fronts{iFront}, 2):-1:1
        for j = 1:n-1
            if crowdingComparison(fronts{iFront}(j+1), fronts{iFront}(j))
                fronts{iFront}([j j+1]) = fronts{iFront}([j+1 j]);
            end
        end
    end
    nextPopulation = structUnion(nextPopulation, fronts{iFront}(1:nIndividuals - size(nextPopulation, 2)));
    
    population = nextPopulation;
    results.elite = fronts{1};
    
    %% Plotting
    if visualize || writeFrames
        fig = figure(1);
        
        clf;
        hold on;
        grid on;
        axis([0 nGenes 0 nGenes]);
        cmap = hsv(size(fronts, 2));
        for f = 1:size(fronts, 2)
            points = [];
            for j = 1:size(fronts{f}, 2)
                points = [points; fronts{f}(j).fitness(1) fronts{f}(j).fitness(2)];
            end
            if ~isempty(points)
                plot(points(:,1), points(:,2), '.-', 'Color', cmap(f, :), 'MarkerSize', 40, 'LineWidth', 2);
                title('NSGA-II');
                xlabel('Leading zeros');
                ylabel('Trailing ones');
            end
        end
        if writeFrames
            saveas(fig, ['images/',num2str(iGeneration, '%02d'),'.png']);
        end
        if visualize
            pause(0.05);
        end
        hold off;
    end
end
end

%% Assign fitness to struct individuals
function population = fitnessAssignment(population, objectives)
for i = 1:size(population, 2)
    for m = 1:size(objectives, 2)
        population(i).fitness(m) = feval(objectives{m}, population(i));
    end
end
end

%% Create fronts from domination
function fronts = dominationSort(population)
fronts = {[]};
for p = 1:size(population, 2)
    population(p).S = [];
    population(p).n = 0;
    for q = 1:size(population, 2)
        if sum(population(p).fitness) > sum(population(q).fitness)
            population(p).S = union(population(p).S, q);
        elseif sum(population(q).fitness) > sum(population(p).fitness)
            population(p).n = population(p).n + 1;
        end
    end
    if population(p).n == 0
        population(p).rank = 1;
        fronts{1} = structUnion(fronts{1}, population(p));
    end
end
i = 1;
while ~isempty(fronts{i})
    Q = [];
    for p = 1:size(fronts{i}, 2)
        for q = 1:size(fronts{i}(p).S, 2)
            population(fronts{i}(p).S(q)).n = population(fronts{i}(p).S(q)).n - 1;
            if population(fronts{i}(p).S(q)).n == 0
                population(fronts{i}(p).S(q)).rank = i + 1;
                Q = structUnion(Q, population(fronts{i}(p).S(q)));
            end
        end
    end
    i = i + 1;
    fronts{i} = Q;
end
end

%% Crowding comparison operator
function ret = crowdingComparison(p, q)
ret = p.rank < q.rank || ((p.rank == q.rank) && (p.dist >= q.dist));
end

%% Assign crowding distance to struct individuals
function population = crowdingDistanceAssignment(population, k)
l = size(population, 2);
for i = 1:l
    population(i).dist = 0;
end

% TODO
f_max = size(population(i).genes,2);
f_min = 0;

for m = 1:size(k,1)
    % Sort by objective
    fitness = zeros(size(population));
    for i = 1:size(population, 2)
        fitness(i) = population(i).fitness(m);
    end
    [~, ind] = sort(fitness);
    population = population(ind);
    
    population(1).dist = inf;
    population(l).dist = inf;
    for i = 2:l-1
        population(i).dist = population(i).dist ...
            + (feval(k{m}, population(i+1)) ...
             - feval(k{m}, population(i-1)))/(f_max-f_min);
    end
end
end

%% Union of two sets of structs
function unified = structUnion(A,B)
unified = [A, B];
[~, ind] = unique(reshape([unified.genes], size(unified(1).genes, 2), size(unified,2))', 'rows');
unified = unified(ind);
end

%% Fitness functions

% Leading zeros
function fitness = leadingZeros(a)
fitness = find(a.genes, 1, 'first') - 1;
if isempty(fitness)
    fitness = size(a.genes, 2);
end
end

% Trailing ones
function fitness = trailingOnes(a)
fitness = find(fliplr(~a.genes), 1, 'first') - 1;
if isempty(fitness)
    fitness = size(a.genes, 2);
end
end