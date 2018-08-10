function [v] = parallel_transport_S2_geo(p,q,u)
% paralell transport of u from p to q along minimizing geodesic on the unit
% sphere S2.

% p:    the coordinate of p, could be Cartesian or math spherical
% q:    the coordinate of q, could be Cartesian or math spherical
% u:    3-by-1 vector, the tangent vector at p to be transported

%% prepare coordinates
[cp,~] = co_convert_S2(p);
[cq,~] = co_convert_S2(q);

%% check if p=q or p=-q
if all(cp==cq) || all(cp==-cq)
    warning('p==q or p==-q');
    v = zeros(3,1);
    return;
end

l2norm = @(u) sqrt(sum(u.^2));


k = cross(cp,cq);
k = k./l2norm(k);
costheta = dot(cp,cq); % 0<theta<pi
sintheta = sqrt(1-costheta^2);
v = costheta*u + sintheta*cross(k,u) + (1-costheta)*dot(k,u)*k;

end
