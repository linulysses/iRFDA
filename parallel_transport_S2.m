function [v] = parallel_transport_S2(pp,qq,u)
% pp: 3-by-n
% qq: 3-by-n
% v : 3-by-n
    [~,n] = size(pp);
    v = zeros(3,n);
    for i = 1:n
        v(:,i) = parallel_transport_S2_geo(pp(:,i),qq(:,i),u(:,i));
    end
end