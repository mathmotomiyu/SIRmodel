%アルミホ法の中で制御系を解くためのルンゲクッタ法、ここではRK41と同じ

function yout = RK42(params,t,x,uo)
h = params.h;
tm = t + h/2;
k1 = h * F1(params,t,x,uo);
k2 = h * F1(params,tm,x + k1/2,uo);
k3 = h * F1(params,tm,x + k2/2,uo);
k4 = h * F1(params,t + h,x + k3,uo);
yout = x + (k1 + k4 + 2.0*(k2 + k3))/6.0;