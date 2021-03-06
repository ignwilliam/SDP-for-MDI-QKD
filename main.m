% find optimal (asymptotic) key rate for phase-encoding MDI-QKD for fixed
% protocol parameters (we do not assume Alice and Bob are equidistant from
% the untrusted party, Charlie)

%% protocol parameters
dist = 20; % distance between Alice and Bob (in km)
tau = 1/2; % the ratio AC/BC distance
distAC = dist * tau; % distance between Alice and Charlie (in km)
distBC = dist * (1-tau); % distance between Bob and Charlie (in km)
det = 0.85;  % detector efficiency
ta = det * 10^( -0.2 * distAC / 10); % transmittivity of AC channel
tb = det * 10^( -0.2 * distBC / 10); % transmittivity of BC channel
pdc = 5E-8; % dark count rate
delta = pi/12; % global phase misalignment

%% optimise key rate
mu_range = [0.01, 0.1]; % range of intensity
res = 10; % resolution of the optimisation grids
[R, mua, mub] = OptimiseIntensity(ta,tb,pdc,delta,mu_range,mu_range,res);