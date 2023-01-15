%射影

function y = Proj(params,u)
y = u ;
if u < params.u1
    y = params.u1;
end
if u > params.u2
    y = params.u2;
end