function R = keyrate(statesA,statesB,ta,tb,pdc,delta)
% calculate the asymptotic key rate of MDI-QKD protocol

% inputs:
    % statesA: states that Alice prepares
    % statesB: states that Bob prepares
    % ta: the transmittivity of AC channel
    % tb: the transmittivity of BC channel
    % pdc: dark count rate of the detectors
    % delta: the global phase misalignment
% output:
    % R: asymptotic key rate
    
    
%% number of states
na = size(statesA,1); % number of Alice's states
nb = size(statesB,1); % number of Bob's states

%% Gram matrix of the input states
lambda = GramInputStates(statesA,statesB);

%% simulate observed statistics
[statL, statR, statFail] = GenerateStatistics(statesA,statesB,ta,tb,pdc,delta);

% detection probability
pdet = 1/4 * (statL(1,1) + statL(1,2) + statL(2,1) + statL(2,2) + ...
              statR(1,1) + statR(1,2) + statR(2,1) + statR(2,2) ); 

% probability of wrong detection
err = 1/4 * (statL(1,2) + statL(2,1) + statR(1,1) + statR(2,2) );

% QBER
Q = err/pdet;

% reshape the arrays into vectors
statL = reshape(statL,[na*nb,1]);
statR = reshape(statR,[na*nb,1]);
statFail = reshape(statFail,[na*nb,1]);

%% phase-error rate estimation
eph = PhaseError(lambda,pdet,statL,statR,statFail,na,nb);

%% key rate calculation
R = max(0, pdet * ( 1 - h2(eph) - h2(Q) ) );