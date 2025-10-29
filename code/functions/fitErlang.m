function [fittedModel,modelParameters] = fitErlang(dataY, dataX)
%fitErlang fits an erlang model on pupil median response, but to be honest
%we stopped caring, so good luck if you wanna continue the job 

%The model
% custom_function = @(params, x) ( (params(1)^floor(params(2))) * exp(-params(1)*(x-nthroot(params(4)/params(3), floor(params(2))-1))) / factorial(floor(params(2))-1))...
%                                    .* ...
%                                ( (params(3)*(x-nthroot(params(4)/params(3), floor(params(2))-1)).^(floor(params(2))-1))-params(4) ) ;
% % Bornes inférieures et supérieures pour les paramètres (L, k, A, B)
% lb = [0.01, 2, 1, 0];  % bornes inférieures 
% ub = [Inf, 10, Inf, Inf];     % bornes supérieures 
%                            
% % Estimation initiale des paramètres (L, k, A, B)
% initial_guess = [0.01, 3, 1, 0];

custom_function = @(params, x) ( ( (params(1)* (x-params(2)).^params(3))-params(4) ).* exp(-params(5) *(x-params(2) ) ) )  ;
% Bornes inférieures et supérieures pour les paramètres
lb = [-10, 0, 2, 20, 0];  % bornes inférieures 
ub = [Inf, 1, 5, 30, Inf];     % bornes supérieures 
                           
% Estimation initiale des paramètres (A, Phi, n, B, lambda)
initial_guess = [1, 0.4, 2, 10, 1];




% Options pour le fit par lsq
opts = optimoptions('lsqcurvefit',...
                    'Algorithm','levenberg-marquardt',...
                    'OptimalityTolerance', 1E-8,...
                    'MaxIterations', 1E6);

% Ajustement de la courbe
modelParameters = lsqcurvefit(custom_function, initial_guess, dataX, dataY, lb, ub, opts);

%Le fit
fittedModel = custom_function(modelParameters, dataX);                           
end

