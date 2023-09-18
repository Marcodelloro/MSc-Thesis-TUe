% Trial to fit a the data from RIVM to the actual SIDARTHE model from the
% G.Giordano paper.
% 8 ODEs: \dot{S} = S*(\alfa*I + \beta*D + \gamma*A + \beta*R)
%         \dot{I} = S*(\alfa*I + \beta*D + \gamma*A + \beta*R) - (\epsilon + \zeta + \lambda)*I
%         \dot{D} = \epsilon*I - (\zeta + \lambda)*D
%         \dot{A} = \zeta*I - (\theta + \mu + \kappa)*A
%         \dot{R} = \zeta*D + \tehta*A -(\mu + \kappa)*R
%         \dot{T} = \mu*A + \mu*R - (\sigma + \tau)*T
%         \dot{H} = \lambda*I + \lambda*D + \kappa*A + \kappa*R + \sigma*T
%         \dot{E} = \tau*T
% 11 Parameters to be estimated
% NOTE THAT sigma and tau are depending on T 
% NOTE THAT all contagion parameters \alfa, \beta, \gamma will be then TIME DEPENDENT based on the VARIANTS

