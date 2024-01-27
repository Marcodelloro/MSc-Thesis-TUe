function [Xk_plus,Pk_plus,Ob] = sidttheEKFobsv(Xaug,P,R,J,L,Qp,C,Zmeas,vars,static_coefs,guessed_vals,nstates)

J_substituted = subs( J, [vars,static_coefs], [Xaug; guessed_vals']' );
A = eval(J_substituted);

% Unmeasured estimate / State and Covariance propagation   
Xk_min = (eye(nstates)+A)*Xaug;
Pk_min = (eye(nstates)+A)*P*(eye(nstates)+A)' + L*Qp*L';

% Computation of covariance and Kalman gain
K = Pk_min*C'/(C*Pk_min*C' + R) ;
Pk_plus = (eye(nstates)- K*C)*Pk_min;

% Computation of state estimate update
Xk_plus = Xk_min + K*(Zmeas - C*Xk_min);

Ob = obsv(A,C);

end