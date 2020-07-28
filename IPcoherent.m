function braket = IPcoherent(bra,ket)
% computes the inner product of (possibly multimode) coherent states

M = size(bra,2); % the number of modes

%% compute the inner product
braket = 1; % initialise

% take the inner product for each mode
for j = 1:M
    braket = braket * exp(-abs(bra(j)-ket(j))^2/2) * exp(1i*imag(conj(bra(j))*ket(j)));
end

