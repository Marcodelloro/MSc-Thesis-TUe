% Trial to fit a the data from RIVM to the actual SIDARTHE model from the
% G.Giordano paper.
% 8 ODEs: \dot{S} = S*(\alfa*I + \beta*D + \gamma*A + \beta*R)
%         \dot{I} = S*(\alfa*I + \beta*D + \gamma*A + \beta*R) - (\epsilon + \zeta + \lambda)*I
%         \dot{D} = \epsilon*I - (\zeta + \lambda)*D
%         \dot{A} = \zeta*I - (\theta + \mu + \kappa)*A
%         \dot{R} = \zeta*D + \theta*A -(\mu + \kappa)*R
%         \dot{T} = \mu*A + \mu*R - (\sigma + \tau)*T
%         \dot{H} = \lambda*I + \lambda*D + \kappa*A + \kappa*R + \sigma*T
%         \dot{E} = \tau*T
% 11 Parameters to be estimated
% NOTE THAT sigma and tau are depending on T 
% NOTE THAT all contagion parameters \alfa, \beta, \gamma will be then TIME DEPENDENT based on the VARIANTS

clear all
close all
clc
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')

load("cleandata.mat")
import casadi.*

%% Model creation and variables

% The model has eight different states + 16 different parameters to be estimated

S = SX.sym('S');
I = SX.sym('I');
D = SX.sym('D');
A = SX.sym('A');
R = SX.sym('R');
T = SX.sym('T');
H = SX.sym('H');
E = SX.sym('E');
x =[S; I; D; A; R; T; H; E]; % vectors with all the states

alfa = SX.sym('alfa');
beta = SX.sym('beta');
gamma = SX.sym('gamma');
delta = SX.sym('delta');
epsi = SX.sym('epsilon');
zeta = SX.sym('zeta');
eta = SX.sym('eta');
theta = SX.sym('theta');
mu = SX.sym('mu');
nu = SX.sym('nu');
rho = SX.sym('rho');
sigma = SX.sym('sigma');
xi = SX.sym('xi');
kappa = SX.sym('kappa');
lambda = SX.sym('lambda');
tau = SX.sym('tau');
parms = [alfa; beta; gamma; delta; epsi; zeta; eta; theta; ...
         mu; nu; rho; sigma; xi; kappa; lambda; tau]; % vectors with all the parameters

T = 57; % Time horizon - 57 weeks

% Model equations
xdot = [    ( S*( alfa*I + beta*D + gamma*A + beta*R ) ); ...                               % \dot{S}
            ( S*(alfa*I + beta*D + gamma*A + beta*R) - (epsi + zeta + lambda)*I ); ...      % \dot{I}
            ( epsi*I - (zeta + lambda)*D ); ...                                             % \dot{D}
            ( zeta*I - (theta + mu + kappa)*A ); ...                                        % \dot{A}
            ( zeta*D + theta*A - (mu + kappa)*R ); ...                                      % \dot{R}
            ( mu*A+mu*R - (sigma + tau)*T ); ...                                            % \dot{T}
            ( lambda*I + lambda*D + kappa*A + kappa*R + sigma*T ); ...                      % \dot{H}
            ( tau*T) ];                                                                     % \dot{E}

% Implementation of the integrator
ode = struct('x',x,'p',parms,'ode',x_dot);     
F = integrator('F', 'idas', ode);



