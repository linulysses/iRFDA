function [cp,sp] = co_convert_S2(pp)

    [~,n] = size(pp);
    
    sp = zeros(2,n);
    cp = zeros(3,n);
    for j = 1:n
        p = pp(:,j);
        if length(p) == 3  % the p is given in Cartesian coordiante
            cp(:,j) = p;
            sp(2,j) = acos(cp(3,j));
            if cp(3,j) < 1 && cp(3,j) > -1  % not north/south polar, cp(3)=cos(varphi)
                sin_varphi = sqrt(1-cp(3,j)^2); % sin(varphi), must be positive, as 0<varphi<pi
                if cp(2,j) >= 0
                    sp(1,j) = acos(cp(1,j)/sin_varphi);
                else 
                    sp(1,j) = acos(-cp(1,j)/sin_varphi) + pi;
                end
            else
                sp(1,j) = pi/2;
            end
        else  % p is given the math spherical coordinate
            sp(:,j) = p;
            cp(:,j) = [cos(sp(1,j))*sin(sp(2,j));  ...
                sin(sp(1,j))*sin(sp(2,j)); cos(sp(2,j))];
        end
    end
end