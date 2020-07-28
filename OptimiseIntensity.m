function [R,mua,mub] = OptimiseIntensity(ta,tb,pdc,delta,mua_range,mub_range,res)
% optimise the intensity of Alice and Bob
% we create a grid corresponding to different intensity configuration and
% find the maximum key rate in that grid

% inputs:
    % ta: transmittivity of AC channel
    % tb: transmittivity of BC channel
    % pdc: dark count rate of the detectors
    % delta: global phase misalignment
    % mua_range: the range for Alice's intensities
    % mub_range: the range for Bob's intensities
    % res: the resolution of the optimisation grid
% outputs:
    % R: optimal key rate
    % mua: optimal intensity for Alice
    % mub: optimal intensity for Bob

    
%% create possible combinations of intensities
mua_vec = linspace(mua_range(1),mua_range(2),res);
mub_vec = linspace(mub_range(1),mub_range(2),res);

k = CombVec(mua_vec,mub_vec);
ka = k(1,:);
kb = k(2,:);

nk = size(k,2);

%% input states
alpha = sqrt(ka);
beta = sqrt(kb);
statesA = [alpha; -alpha; 1i*alpha; -1i*alpha];
statesB = [beta; -beta; 1i*beta; -1i*beta];

%% calculate key rate
K = zeros(1,nk);
parfor i = 1:nk
    ai = statesA(:,i);
    bi = statesB(:,i);
    
    K(i) = keyrate(ai,bi,ta,tb,pdc,delta); % key rate calculation for each intensity config.
end

% finding the optimal key rate
[R, I] = max(K);

% optimal intensities
mua = ka(I);
mub = kb(I);