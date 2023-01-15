%制御系

function yout = F1(params,t,xi,u)

y1 = params.m - params.m * xi(1) - params.c * xi(1) * xi(2) - u * xi(1);
y2 = - params.m * xi(2) + params.c * xi(1) * xi(2) - params.d * xi(2);
y3 = - params.m * xi(3) + u * xi(1) + params.d * xi(2) + t * 0;

yout = [y1;y2;y3];