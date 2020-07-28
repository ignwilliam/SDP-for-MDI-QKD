function lambda = GramInputStates(statesA,statesB)
% computes the Gram matrix of the input states (assuming coherent states)

% inputs:
    % statesA: the input states' amplitudes for Alice
    % statesB: the input states' amplitudes for Bob
% output: 
    % lambda: Gram matrix of the input states
    
    
%% dimensions
na = size(statesA,1); % number of states of Alice
nb = size(statesB,1); % number of states of Bob
D = na*nb; % total number of possible pairs

%% initialise the Gram matrix
lambda = zeros(D,D);

%% compute inner products
for i = 0:na-1
    for j = 0:nb-1
        for k = 0:na-1
            for l = 0:nb-1
                x = 1 + i*nb + j;
                y = 1 + k*nb + l;
                lambda(x,y) = IPcoherent(statesA(i+1,:),statesA(k+1,:)) * ...
                              IPcoherent(statesB(j+1,:), statesB(l+1,:));
            end
        end
    end
end