function deriv=deriv(X,params)
% dS/dt = - beta * S * I
% dI/dt = + beta * S * I - gamma * I
% dR/dt =                + gamma * I
%S=X(1); I=X(2); R=X(3);
% params=[beta gamma];
beta=params(1); gamma=params(2); N=params(3); 
S=X(1); I=X(2); R=X(3);

dSdt = - beta/N * S * I; 
dIdt = + beta/N * S * I - gamma * I;
dRdt =                  + gamma * I;
% Para mejor performance, se podria reemplazar la 
% ecuacion para R por R = N - S - I

deriv=[dSdt dIdt dRdt];
       
