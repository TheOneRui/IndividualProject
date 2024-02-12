syms f(t)
f(t) = 50 - (-0.015)*[1+(6.011*exp(-0.438*t)*sin((0.348*t)-9.60))];


syms y(t)
y(t) = piecewise(t < 0,0,t >= 0,-0.3);

t = (-2:0.01:20)';

plot(t,f(t))
hold on
plot(t,af(t))
hold off

syms df(t)