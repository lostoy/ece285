%    Goal: show how to use MSBL in noisy case and show that its lambda
%          learning rule is not robust in such low SNR cases.
%
%  Author: Zhilin Zhang (z4zhang@ucsd.edu)
%    Date: March 05, 2011
% Version: 1.0


clear all;  

% Experiment Assignment

iterNum   = 1 ;          % Trial number (i.e. number of repeating the experiment)

% Problem dimension
N = 40;                  % Row number of the dictionary matrix 
M = 150;                 % Column number of the dictionary matrix
L = 4;                   % Number of measurement vectors
K = 10;                  % Number of nonzero rows (i.e. source number) in the solution matrix
SNR = 10;                % SNR                        


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
    signal = Phi * Wgen;

    % Observation noise   
    stdnoise = std(reshape(signal,N*L,1))*10^(-SNR/20);
    noise = randn(N,L) * stdnoise;

    % Noisy signal
    Y = signal + noise;


    %============= Run MSBL using the lambda learning rule ======================
    lambda0 = 1e-2;          % Initial value for the regularization parameter. 
                             %  In noiseless cases, you can set lambda =
                             %  1e-10 (or any other very small values) and 
                             %  set learn_Lambda = 0, which can lead to excellent performance.
    
    learn_Lambda0 = 1;        % Using its lambda learning rule to learn an (sub-)optimal lambda. 
                             %  When SNR < 20 dB, the lambda learning rule may not be
                             %  robust. In this case, you probabaly need to use
                             %  other methods to find a lambda value and set learn_Lambda = 0.
                             
    tic;
    [Weight0,gamma_est0,gamma_used0,count0] = MSBL(Phi,Y, lambda0, learn_Lambda0);
    time0 = toc;
    TIME0(it) = time0;

    % Failure rate
    F0 = perfSupp(Weight0,indice,'firstlargest', K);      
    fail_MSBL0(it) = (F0~=1);      
    
    % MSE
    perf_MSBL0(it) = (norm(Wgen - Weight0,'fro')/norm(Wgen,'fro'))^2;  
    
    fprintf('   MSBL(using lambda rule): time = %5.2f; Findex = %3.2f, Ave-MSE = %3.2f%%; Ave-Fail_Rate = %4.3f%%; Ave-Time = %4.3f\n',...
        time0,F0,sum(perf_MSBL0)/it*100,sum(fail_MSBL0)/it*100,sum(TIME0)/it);


    %============= Run MSBL but fix the lambda to the true noise variance ==============
    lambda1 = stdnoise^2;    % fix lambda to the true noise variance
    learn_Lambda1 = 0;       % do not use its lamda learning rule
                             
    tic;
    [Weight1,gamma_est1,gamma_used1,count1] = MSBL(Phi,Y, lambda1, learn_Lambda1);
    time1 = toc;
    TIME1(it) = time1;

    % Failure rate
    F1 = perfSupp(Weight1,indice,'firstlargest', K);      
    fail_MSBL1(it) = (F1~=1);      
    
    % MSE
    perf_MSBL1(it) = (norm(Wgen - Weight1,'fro')/norm(Wgen,'fro'))^2;  
    
    fprintf('MSBL(given noise variance): time = %5.2f; Findex = %3.2f, Ave-MSE = %3.2f%%; Ave-Fail_Rate = %4.3f%%; Ave-Time = %4.3f\n',...
        time1,F1,sum(perf_MSBL1)/it*100,sum(fail_MSBL1)/it*100,sum(TIME1)/it);
end

fprintf('\nNow you can see, in low SNR cases the lambda learning rule is not robust\n');
fprintf('In such cases you''d better fix lambda to a value, which is around the true noise variance\n');


