function [model] = switching_pandemic(ODEs_0,coefs_0,ICUs)
    
    import casadi.*
    model = NosnocModel();
   
    %% Numer of ODE layers
    
    % Differential states
    S = MX.sym('S') ;
    I = MX.sym('I') ;
    D = MX.sym('D') ;
    T1 = MX.sym('T1') ;
    T2 = MX.sym('T2') ;
    H = MX.sym('H') ;
    E = MX.sym('E') ;

    % Coefficients
    alpha = MX.sym('alpha') ;
    beta = MX.sym('beta') ;
    gamma = MX.sym('gamma') ;
    delta1 = MX.sym('delta1') ;
    delta2 = MX.sym('delta2') ;
    epsi = MX.sym('epsi') ;
    tau1 = MX.sym('tau1') ;
    tau2 = MX.sym('tau2') ;
    sigma1 = MX.sym('sigma1') ;
    sigma2 = MX.sym('sigma2') ;
    lambda = MX.sym('lambda') ;

    ODEs = [ S; I; D; T1; T2; H; E ] ;
    coefs = [ alpha; beta; gamma; delta1; delta2; epsi; tau1; tau2; sigma1; sigma2; lambda ] ;          
    
    model.x0 = [ODEs_0];
    model.x = [ODEs];
    
    %% Switching Functions

    model.c =  delta2*D - (sigma2 + tau2)*T2 + epsi*T1 - ICUs ;
  
    % sign matrix for the modes
    model.S = [1;-1];

    %% Modes of the ODEs layers 
  
    f_1 = [ -S * (alpha*I + beta*D) ;...
            S * (alpha*I + beta*D) - (gamma+lambda)*I ;...
            I*gamma - D*(lambda + delta1 + delta2) ;...
            delta1*D - (sigma1 + tau1 - epsi)*T1 ;...
            delta2*D - (sigma2 + tau2)*T2 + epsi*T1 ;...
            (I+D)*lambda + T1*sigma1 + T2*sigma2 ;...
            T1*tau1 + T2*tau2
           ];
    
  
    f_2 = [ -S * (alpha*I + beta*D) ;...
            S * (alpha*I + beta*D) - (gamma+lambda)*I ;...
            I*gamma - D*(lambda + delta1 + delta2) ;...
            delta1*D - (sigma1 + tau1 - epsi)*T1 ;...
            delta2*D - (sigma2 + tau2)*T2 ;...
            (I+D)*lambda + T1*sigma1 + T2*sigma2 ;...
            T1*tau1 + T2*tau2 + epsi*T1
           ];


    F = [f_1 f_2];
    
    
    model.F = F ;
end