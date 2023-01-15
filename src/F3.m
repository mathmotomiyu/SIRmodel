%随伴方程式

function yout = F3(params,t,xi,u,pi)

y1 = (params.m  + params.c * xi(2) + u) * pi(1) - params.c * xi(2) * pi(2) - u * pi(3);
y2 = params.c * xi(1) * pi(1) + (params.m - params.c * xi(1) + params.d) * pi(2) - params.d * pi(3) - 1;
y3 = params.m * pi(3) + t * 0;

yout = [y1;y2;y3];