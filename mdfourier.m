% Fourier basis on [0,1]
% m:  # design points
% K:  # fPCs
% A:  orthonormal matrix, d-by-d
% b:  d-m-by-K matrix
% t:  design points
function [b,tt] = mdfourier(m,K,A,t)
    if nargin == 4
        tt = t;
    else
        tt = linspace(0,1,m);
    end
    d = size(A,1);
    b = zeros(d,m,K);
    L = ceil(K/d);
    bb = fourier(m,L);
    
    for j = 1:L
        for i = 1:d
            k = (j-1)*d + i;
            if k > K; break; end
            w = A(:,i) * bb(:,j)'; % d by m
            b(:,:,k) = w;
        end
    end
end
