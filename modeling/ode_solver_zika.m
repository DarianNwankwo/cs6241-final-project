%% ODE Solver for Test Systems - Zika Problem

%Solves system of ODES in matlab using ODE 45

%% Model 1 - Populations of FR, MR, FW, MW

%Params 
alpha = 0.75;
r = 0.5; 
K = 500; 
delta = 0.05; 
b = 1;
AWF = 0.5; 
AWM = 0.5; 

times = [0, 50];
ics = [10, 100, 10, 10];

[t, mos] = ode45(@(t, X) odesolveModel1(t, X, r, alpha, K, delta, b, AWF, AWM), times, ics);

plot(t, mos(:,1), t, mos(:,2), t, mos(:,3), t, mos(:,4))
legend('FR','MR', 'FW', 'MW')

%% Function Definitions  
function dXdt = odesolveModel1(t, X, r, alpha, K, delta, b, AWF, AWM)
    FR = X(1); 
    MR = X(2); 
    FW = X(3);
    MW = X(4); 
    F = FR + FW;
    M = MR + MW; 
     
    frEqn = r * alpha * (1 - ((F + M) / K)) * X(1) * (MR / (b + M)) - delta * FR;
    mrEqn = r * (1 - alpha) * (1 - ((F + M) / K)) * FR * (MR / (b + M)) - delta * MR;
    fwEqn = r * alpha * (1 - ((F + M) / K)) * FW * (M / (b + M)) - delta * FW + AWF;
    mwEqn = r * (1 - alpha) * (1 - ((F + M) / K)) * FW * (M / (b + M)) - delta * MW + AWM; 
    dXdt = [ frEqn; mrEqn; fwEqn; mwEqn ];
end 
