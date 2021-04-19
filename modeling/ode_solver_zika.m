%% ODE Solver for Test Systems - Zika Problem

%Solves system of ODES in matlab using Forward Euler and ODE 45

%% Model 1 - Populations of FR, MR, FW, MW

%Params
alpha = 0.75;
r = 0.5;
K = 500;
delta = 0.05;
b = 1;
AWF = 5;
AWM = 5;

%ics
fr0 = 10;
mr0 = 100;
fw0 = 10;
mw0 = 10;
rpop0 = 110;

times = [0, 50];
%ics = [10, 100, 10, 10]; %ICs for 4d case
ics = [110, 10, 10];
%% Forward Euler Manual Solve
t0 = 0;
tf = 50;
dt = 1; %check DT
n = 50; %number of forward euler steps to take

tvals = linspace(t0,tf,n);

fr = zeros(n,1);
mr = zeros(n,1);
fw = zeros(n,1);
mw = zeros(n,1);
rpop = zeros(n,1);

%set initial values
fr(1) = fr0;
mr(1) = mr0;
fw(1) = fw0;
mw(1) = mw0;
rpop(1) = rpop0;

%Loop to compute via forward euler method
for i = 2:n
    %set up M and F
    
    %for 4d
    %M = mr(i-1) + mw(i-1);
    %F = fr(i-1) + fw(i-1);
    
    M = (1-alpha)* rpop(i-1) + mw(i-1);
    F = alpha * rpop(i-1) + fw(i-1);
    
    rpop(i) = rpop(i-1) + dt * (r * (1 - ((F + M) / K)) * ((alpha * rpop(i-1) * (1 - alpha) * rpop(i-1) / (b + M))) - delta * rpop(i-1));
    fw(i) = fw(i-1) + dt * (r * alpha * (1 - ((F + M) / K)) * fw(i-1) * (M / (b + M)) - delta * fw(i-1) + AWF);
    mw(i) = mw(i-1) + dt * (r * (1 - alpha) * (1 - ((F + M) / K)) * fw(i-1) * (M / (b + M)) - delta * mw(i-1) + AWM);
    
    %For 4d
    %fr(i) = fr(i-1) + dt * (r * alpha * (1 - ((F + M) / K)) * fr(i-1) * (mr(i-1) / (b + M)) - delta * fr(i-1));
    %mr(i) = mr(i-1) + dt * (r * (1 - alpha) * (1 - ((F + M) / K)) * fr(i-1) * (mr(i-1) / (b + M)) - delta * mr(i-1));
end

f1 = figure();
%plot(tvals, fr, tvals, mr, tvals, fw, tvals, mw) %for 4d
plot(tvals, rpop, tvals, fw, tvals, mw)
title('Mosquito Populations (Forward Euler Method)')
%legend('FR','MR', 'FW', 'MW') %for 4d
legend('R', 'FW', 'MW')

%% ODE 45 Check
[t, mos] = ode45(@(t, X) odesolveModel1(t, X, r, alpha, K, delta, b, AWF, AWM), times, ics);

f2 = figure();
%plot(t, mos(:,1), t, mos(:,2), t, mos(:,3), t, mos(:,4)) %for 4d
plot(t, mos(:,1), t, mos(:,2), t, mos(:,3))
title('Mosquito Populations (ODE45 Built-In)')
legend('R', 'FW', 'MW')

%% Function Definitions
function dXdt = odesolveModel1(t, X, r, alpha, K, delta, b, AWF, AWM)
    R = X(1);
    FW = X(2);
    MW = X(3);
    
    %For 4D
    %FR = X(1);
    %MR = X(2);
    %FW = X(3);
    %MW = X(4);
    %F = FR + FW;
    %M = MR + MW;
    
    F = alpha * R + FW;
    M = (1-alpha) * R + MW;
    
    rEqn = r * (1 - ((F + M) / K)) * ((alpha * R * (1 - alpha) * R / (b + M))) - delta * R;
    fwEqn = r * alpha * (1 - ((F + M) / K)) * FW * (M / (b + M)) - delta * FW + AWF;
    mwEqn = r * (1 - alpha) * (1 - ((F + M) / K)) * FW * (M / (b + M)) - delta * MW + AWM;
    
    %For 4D
    %frEqn = r * alpha * (1 - ((F + M) / K)) * X(1) * (MR / (b + M)) - delta * FR;
    %mrEqn = r * (1 - alpha) * (1 - ((F + M) / K)) * FR * (MR / (b + M)) - delta * MR;
    %dXdt = [ frEqn; mrEqn; fwEqn; mwEqn ];
    
    %Return soln vector
    dXdt = [ rEqn; fwEqn; mwEqn ];
end
