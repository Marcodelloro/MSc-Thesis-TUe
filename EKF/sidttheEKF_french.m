function [Xk_plus,Pk_plus,Ob,Abar,Bbar,Cbar,T,k,A,C] = sidttheEKF_2(Xaug,P,R,w_dist,Zmeas)
 
lambda = 0.05; % chosen value of lambda, kept constant
beta = 0.02;
gamma = 0.3;

% Linearizd matrices (augmented state)
% A = linearization of the dynamics function with augmented state

A = [ -(Xaug(8)*Xaug(2) + beta*Xaug(3))               -(Xaug(1)*Xaug(8))                   -(Xaug(1)*beta)                      0                       0                     0                       0              -(Xaug(1)*Xaug(2))          0                   0                   0                    0                   0                     0                     0                   ;
       (Xaug(8)*Xaug(2) + beta*Xaug(3))     (Xaug(1)*Xaug(8) - gamma - lambda)              (Xaug(1)*beta)                      0                       0                     0                       0               (Xaug(1)*Xaug(2))          0                   0                   0                    0                   0                     0                     0                   ;
                    0                                     gamma                    -(lambda + Xaug(9) + Xaug(10))               0                       0                     0                       0                      0                -Xaug(3)            -Xaug(3)               0                    0                   0                     0                     0                   ;  
                    0                                       0                                     Xaug(9)       -(Xaug(12) + Xaug(14) + Xaug(11))       0                     0                       0                      0                 Xaug(3)               0                -Xaug(4)              -Xaug(4)              0                 -Xaug(4)                  0                   ;
                    0                                       0                                     Xaug(10)                  Xaug(11)          -(Xaug(13)+Xaug(15))            0                       0                      0                   0                 Xaug(3)             Xaug(4)                0                -Xaug(5)                 0                  -Xaug(5)               ;     
                    0                                    lambda                                   lambda                    Xaug(12)                 Xaug(13)                 0                       0                      0                   0                   0                    0                  Xaug(4)            Xaug(5)                 0                     0                   ;
                    0                                       0                                        0                      Xaug(14)                 Xaug(15)                 0                       0                      0                   0                   0                    0                   0                   0                  Xaug(4)                Xaug(5)               ;
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ;
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ;
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ;
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ;
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ;                  
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ; 
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ;
                    0                                       0                                        0                          0                       0                     0                       0                      0                   0                   0                    0                   0                   0                     0                     0                   ];

                    
% L = linearization of the dynamics function wrt disturbances 
% 10 coefficients with disturbances - 17 total states

L = [ zeros(7,8); eye(8) ];

% L = linearization of the measurement function
% States measured are only D-T1-T2-H-E

C = [zeros(5,2) eye(5) zeros(5,8)];

Qp = diag([w_dist]);   % Process noise covariance

% Unmeasured estimate / State and Covariance propagation   
Xk_min = (eye(15)+A)*Xaug;
Pk_min = (eye(15)+A)*P*(eye(15)+A)' + L*Qp*L';

% Computation of covariance and Kalman gain
K = Pk_min*C'/(C*Pk_min*C' + R) ;
Pk_plus = (eye(15)- K*C)*Pk_min;

% Computation of state estimate update
Xk_plus = Xk_min + K*(Zmeas - C*Xk_min);

Ob = obsv(A,C);

B = zeros(15,1);
[Abar,Bbar,Cbar,T,k] = obsvf(A,B,C);


end