function v = lorentz(x,y)

v = x(1,:).*y(1,:)+x(2,:).*y(2,:)-x(3,:).*y(3,:);
    