function [R] = XrMSE_SPD(obj,n,X,XX,Ks)

    mfd = obj.mfd;
    matD = sqrt(mfd.D);
    m = size(X,2);
    K = length(Ks);
    
    
    RR = zeros(m,n,K);
    for i = 1:n
        for j = 1:m
            x = X(:,j,i);
            x0 = zeros(matD,matD);
            for u = 1:matD;
                for v = 1:matD;
                    x0(u,v) = x((u-1)*matD+v);
                end
            end

            for k = 1:K
                y = XX{Ks(k)}(:,j,i);
                y0 = zeros(matD,matD);
                for u = 1:matD;
                    for v = 1:matD;
                        y0(u,v) = y((u-1)*matD+v);
                    end
                end

                [u,s] = svd(x0);
                [v,w] = svd(y0);
                logx = u * diag(log(diag(s))) * u';
                logy = v * diag(log(diag(w))) * u';
                w = (logx - logy).^2;
                RR(j,i,k) = sum(w(:)); % square of frobenius norm
            end
        end
    end
    R = sqrt(squeeze(mean(mean(RR,1),2)))';
end