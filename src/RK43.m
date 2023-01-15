%随伴方程式を解くためのルンゲクッタ法

function yout = RK43(params,t,xi,uo,pi)
h1 = -params.h;
tm = t + h1/2;
k1 = h1 * F3(params,t,xi,uo,pi);
k2 = h1 * F3(params,tm,xi,uo,pi + k1/2);
k3 = h1 * F3(params,tm,xi,uo,pi + k2/2);
k4 = h1 * F3(params,t + h1,xi,uo,pi + k3);
yout = pi + (k1 + k4 + 2.0*(k2 + k3))/6.0;