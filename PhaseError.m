function eph = PhaseError(lambda,pdet,statL,statR,statFail,na,nb)
% computes the phase-error rate
% inputs:
    % lambda: Gram matrix of the input states
    % pdet: detection probability
    % statL: conditional probability P(L|a,b)
    % statR: conditional probability P(R|a,b)
    % statFail: conditional probability P(fail|a,b)
    % na: number of Alice's states
    % nb: number of Bob's states
% output:
    % eph: phase-error rate

%% compute phase-error rate via SDP
N = na*nb; % number of possible pairs
D = 3*N; % dimension of the Gram matrix

cvx_begin sdp
    variable G(D,D) semidefinite hermitian
    expression overlap(N,N);
    expression pL(N);
    expression pR(N);
    expression pFail(N);
    expression ephase;
    
    %% define expressions for constraints
    % overlap
    for i = 1:N
        for j = 1:N
            for k = 1:3
                overlap(i,j) = overlap(i,j) + G(3*(i-1) + k, 3*(j-1) + k);
            end
        end
    end
    
    % observed statistics
    for i = 1:N
        pL(i) = G(3*(i-1) + 1, 3*(i-1) + 1);
        pR(i) = G(3*(i-1) + 2, 3*(i-1) + 2);
        pFail(i) = G(3*(i-1) + 3, 3*(i-1) + 3);
    end
    
    %% define expression for objective function
    % compute index for relevant elements
    % formula: index = 3*(a*nb + b) + z
    L00 = 1; % for a = 0, b = 0, z = 1
    R00 = 2; % for a = 0, b = 0, z = 2
    L01 = 4; % for a = 0, b = 1, z = 1
    R01 = 5; % for a = 0, b = 1, z = 2
    L10 = 3*nb + 1; % for a = 1, b = 0, z = 1
    R10 = 3*nb + 2; % for a = 1, b = 0, z = 2
    L11 = 3*(nb + 1) + 1; % for a = 1, b = 1, z = 1
    R11 = 3*(nb + 1) + 2; % for a = 1, b = 1, z = 2
    
    % define phase-error (c.f. Annex A of PRA 99, 062332)
    ephase = 1/2 + 1/(4*pdet) *  real( G(L00,L11) - G(R00,R11) - G(L01,L10) + G(R01,R10) );
    
    %% optimisation proper
    maximize ephase
    subject to
        ephase <= 1/2;
        
        % overlap constraints
        for i = 1:N
            for j = 1:N
                overlap(i,j) == lambda(i,j);
            end
        end
        
        % statistics constraints
        for i = 1:N
            pL(i) == statL(i);
            pR(i) == statR(i);
            pFail(i) == statFail(i);
        end
cvx_end

eph = ephase;