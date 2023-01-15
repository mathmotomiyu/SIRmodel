%制御系を解くルンゲクッタ法

function yout = RK41(params,t,xi,uo)
h = params.h;
tm = t + h/2;
k1 = h * F1(params,t,xi,uo);
k2 = h * F1(params,tm,xi + k1/2,uo);
k3 = h * F1(params,tm,xi + k2/2,uo);
k4 = h * F1(params,t + h,xi + k3,uo);
yout = xi + (k1 + k4 + 2.0*(k2 + k3))/6.0;