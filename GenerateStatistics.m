function [statL, statR, statFail] = GenerateStatistics(statesA,statesB,ta,tb,pdc,delta)
% simulates the observed statistics in the experiment
% computes the conditional probability P(z|a,b) that Charlie announces 'z'
% given that Alice chooses state 'a' and Bob chooses state 'b'

% inputs:
    % statesA: the coherent states' amplitudes prepared by Alice
    % statesB: the coherent states' amplitudes prepared by Bob
    % ta: the channel transmittivity between Alice and Charlie
    % tb: the channel transmittivity between Bob and Charlie
    % pdc: dark count rate of each detector
    % delta: mismatch in the global phase of Alice and Bob
% outputs:
    % statL: conditional probability of Charlie announcing 'L'
    % statR: conditional probability of Charlie announcing 'R'
    % statFail: conditional probability of Charlie announcing 'Fail'
    
    
%% apply delta phase-shift the states of Bob
statesB = exp(1i*delta) * statesB;

%% the states at Charlie's lab
% the amplitudes BEFORE the 50/50 BS
alpha = statesA * sqrt(ta);
beta = statesB * sqrt(tb);

% create a meshgrid for different pair of states
[a,b] = meshgrid(alpha,beta); 

% the amplitudes AFTER the 50/50 BS
aL = (a + b)/sqrt(2);
aR = (a - b)/sqrt(2);

% the intensities AFTER the 50/50 BS
muL = abs(aL).^2;
muR = abs(aR).^2;

%% observed statistics
% probability of click for each detector
pL = 1 - (1-pdc) .* exp(-muL);
pR = 1 - (1-pdc) .* exp(-muR);

% probability of Charlie announcing 'L','R', or 'Fail'
statL = pL .* (1-pR);
statR = (1-pL) .* pR;
statFail = pL .* pR + (1-pL) .* (1-pR);
