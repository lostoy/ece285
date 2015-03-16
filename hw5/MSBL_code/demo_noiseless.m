%    Goal: Show how to use MSBL in noiseless case
%  Author: Zhilin Zhang (z4zhang@ucsd.edu)
%    Date: March 05, 2011
% Version: 1.0


clear all;  

% Experiment Assignment

iterNum   = 50;          % Trial number (i.e. number of repeating the experiment)

% Problem dimension
N = 25;                  % Row number of the dictionary matrix 
M = 100;                 % Column number of the dictionary matrix
L = 4;                   % Number of measurement vectors
K = 12;                  % Number of nonzero rows (i.e. source number) in the solution matrix
                    

for it = 1 : iterNum
    fprintf('\n\nTrial #%d:\n',it);

    % Generate dictionary matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(N,M);
    Phi = Phi./(ones(N,1)*sqrt(sum(Phi.^2)));
   
    % Generate the K nonzero rows     
    nonzeroW = randn(K,L);

    % Locations of nonzero rows are randomly chosen
    ind = randperm(M);
    indice = ind(1:K);
    Wgen = zeros(M,L);
    Wgen(indice,:) = nonzeroW;

    % Noiseless signal
    Y = Phi * Wgen;


    %============================ Run MSBL ==========================
    lambda = 1e-10;          % Initial value for the regularization parameter. 
                             %  In noiseless cases, you can set lambda =
                             %  1e-10 (or any other very small values) and 
                             %  set learn_Lambda = 0, which can lead to excellent performance.
    
    learn_Lambda = 0;        % Using its lambda learning rule to learn an (sub-)optimal lambda. 
                             %  When SNR < 20 dB, the lambda learning rule may not be
                             %  robust. In this case, you probabaly need to use
                             %  other methods to find a lambda value and set learn_Lambda = 0.
                             
    tic;
    [Weight3,gamma_est3,gamma_used3,count3] = MSBL(Phi,Y, lambda, learn_Lambda);
    time3 = toc;
    TIME3(it) = time3;

    % Failure rate
    F3 = perfSupp(Weight3,indice,'firstlargest', K);      
    fail_MSBL(it) = (F3~=1);      
    
    % MSE
    perf_MSBL(it) = (norm(Wgen - Weight3,'fro')/norm(Wgen,'fro'))^2;  
    
    fprintf(' MSBL: time = %5.2f; Findex = %3.2f, Ave-MSE = %3.2f%%; Ave-Fail_Rate = %4.3f%%; Ave-Time = %4.3f\n',...
        time3,F3,sum(perf_MSBL)/it*100,sum(fail_MSBL)/it*100,sum(TIME3)/it);

end




