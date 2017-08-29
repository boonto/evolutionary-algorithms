%% CMA-ES algorithm implementation
%  results  Struct with the best solution for every generation
%  stats    Struct with various stats over all generations
% 
%  params.nMaxGenerations  Termination criteria for generation loop
%  params.mu               Number of individuals 
%  params.lambda           Number of offspring
%  params.n                Number of dimensions
%  params.function         Objective function to use
function [results, stats] = cma_es(params)

% Remap parameters
nMaxGenerations = params.nMaxGenerations;
mu = params.mu;
lambda = params.lambda;
n = params.n;
func = params.function;

%% Initialize static variables
weights = linspace(mu, 1, mu);
weights = weights/sum(weights);
mu_eff = 1/(sum(weights.^2)); % Variance effective selection mass
c_c = 1/n; % Cumulation evolution path scale
c_1 = 2/(n^2); % Rank-1-update scale
c_mu = min(mu_eff/(n^2), 1-c_1); % Rank-mu-update scale
c_sigma = 1/n; % Sigma update scale
d_sigma = 1; % Damping factor for sigma update

%% Initialize dynamic variables
% Create mu individuals (consisting of n y values)
x = (rand(n, mu)-0.5) * 4; % range [-2 2]
% Sort individuals by function value
funcVal = func(x);
[~, sorted] = sort(funcVal);
x = x(:, sorted); % Minimization problem

m = sum(x .* weights, 2); % Weight mean of individuals
p_c = zeros(n, 1); % Cumulation evolution path
C = eye(n); % Covariance matrix
p_mu = zeros(n, 1); % Conjugate evolution path
sigma = 0.5; % step-size

%% Generation Loop
for t = 1:nMaxGenerations
    % Generate new generation
    new_x = m + sigma * mvnrnd(zeros(n, 1), C, lambda)';
    new_x(new_x > 2) = 2;
    new_x(new_x < -2) = -2;
    % Sort individuals by function value
    funcVal = func(new_x);
    [funcVal, sorted] = sort(funcVal);
    new_x = new_x(:, sorted(1:mu));
    
    % Update mean
    new_m = sum(new_x .* weights, 2);
        
    % Rank-1-update
    % Update evolution path
    new_p_c = (1 - c_c) * p_c + sqrt(c_c*(2-c_c)*mu_eff)*(new_m-m)/sigma;
    C_1 = new_p_c * new_p_c';
    % Rank-mu-update
    y = (new_x - m)/sigma;
    C_mu = zeros(n, n, mu);
    for i = 1:mu
        C_mu(:, :, i) = weights(i) .* (y(:,i) * y(:,i)');
    end
    C_mu = sum(C_mu, 3);
    % Update covariance matrix
    new_C = (1 - c_1 - c_mu) * C + c_1 * C_1 + c_mu * C_mu;
    
    % Update sigma
    [Q, D] = eig(C); % Eigendecomposition
    new_p_mu = (1 - c_sigma) * p_mu + sqrt(c_sigma * (2 - c_sigma) * mu_eff)...
        * (Q*sqrt(D)*Q') *... % Recombine matrix
        (new_m - m)/sigma;
    new_sigma = sigma * exp((c_sigma/d_sigma)*(norm(new_p_mu)/...
        mean(sqrt(sum(mvnrnd(zeros(n), eye(n)).^2, 1)))...
        -1));
    
    % Reassign dynamic variables
    x = new_x;
    m = new_m;
    p_c = new_p_c;
    C = new_C;
    p_mu = new_p_mu;
    sigma = new_sigma;
    
    % Record stats
    results.bestIndividual = x(:,1);
    results.bestIndividual;
    % minimization problem
    stats.bestFuncVal(t) = min(funcVal);
    stats.medianFuncVal(t) = median(funcVal);
end
end