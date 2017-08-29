%% Assessment 3: Function Minimization with CMA-ES

% reset the workplace
clear variables;
close all;

if 0
    % Parametrization
    params.nMaxGenerations = 100;
    params.lambda = 10;
    params.mu = params.lambda/2;

    nSamples = 100;

    params.n = 2;
    params.function = @frosen;

    parfor iSample = 1:nSamples
        [results{1, iSample}, stats{1, iSample}] = cma_es(params);
    end

    params.n = 12;
    params.function = @frosen;

    parfor iSample = 1:nSamples
        [results{2, iSample}, stats{2, iSample}] = cma_es(params);
    end

    params.n = 2;
    params.function = @frastrigin;

    parfor iSample = 1:nSamples
        [results{3, iSample}, stats{3, iSample}] = cma_es(params);
    end

    params.n = 12;
    params.function = @frastrigin;

    parfor iSample = 1:nSamples
        [results{4, iSample}, stats{4, iSample}] = cma_es(params);
    end

    % Save data
    params.nSamples = nSamples;
    save('FM_data', 'results', 'stats', 'params');
end

%% Visualize data
if 1
    smData = load('FM_data.mat');
    a = cell2mat(smData.stats);
    
    % Extract relevant data
    medianFuncVal = [];
    for iConfiguration = 1:size(a, 1)
        b = median(reshape([a(iConfiguration,:).medianFuncVal], ...
            smData.params.nMaxGenerations, smData.params.nSamples), 2);
        medianFuncVal = [medianFuncVal, b(:)];
    end
     
    % Plot data
    figure(3);
    plot(medianFuncVal(:,1));
    grid on;
    hold on;
    plot(medianFuncVal(:,2));
    plot(medianFuncVal(:,3));
    plot(medianFuncVal(:,4));
    
    title('Function Minimization with CMA-ES');
    xlabel('Generation');
    ylabel('Function Value');
    legend('Rosenbrock 2D', 'Rosenbrock 12D', 'Rastrigin 2D', 'Rastrigin 12D');
    yticks(sort([-600:200:1200 -480 -80]));
    
    c = cell2mat(smData.results);
    bestRo2D = c(1,smData.params.nMaxGenerations).bestIndividual
    bestRo12D = c(2,smData.params.nMaxGenerations).bestIndividual
    bestRa2D = c(3,smData.params.nMaxGenerations).bestIndividual
    bestRa12D = c(4,smData.params.nMaxGenerations).bestIndividual
end

%% Rosenbrock function
function f=frosen(x)
  if size(x,1) < 2; error('dimension must be greater one'); end
  f = 1e2*sum((x(1:end-1,:).^2 - x(2:end,:)).^2,1) + sum((x(1:end-1,:)-1).^2,1);
end
  
%% Rastrigin function
function f=frastrigin(x)
  N = size(x,1);
  f = (N/20)*807.06580387678 - (10 * (N-sum(cos(2*pi*x),1)) + sum(x.^2,1));
  f(any(abs(x) > 5.12)) = 1e2*N;
  f = -f;
end
