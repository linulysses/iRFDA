% Fourier basis on [0,1]
% m:  # design points
% K:  # fPCs
% b:  m-by-K matrix
% t:  design points
function [b,t] = fourier(m,K)
    t = linspace(0,1,m);
    b = zeros(m,K);
    for jj = 1:K
        if mod(jj,2) == 1 % sin[(jj+1)/2*pit]
            b(:,jj) = sqrt(2)*sin((jj+1)*pi*t);
        else
            b(:,jj) = sqrt(2)*cos(jj*pi*t);
        end
    end
end
