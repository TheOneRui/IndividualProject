%Part 1: calculate GSFR(s)
model1 = "function1";
open_system(model1);

%assign the block (blk) to the one changing
blk1 = model1 + "/Pd";
%create initial Pd
Pd = -0.3;
%sets the value of constant block to the parameter Pd
set_param(blk1,"Value","Pd");

out = sim(model1);

time = out.tout;
data = out.simout;

hold on
for k = 0 : 1 : 6
    for j = length(time):-1:1
        if(time(j) <= k)
            data(j) = 0;
        end
    end
    plot(time,data);
end

hold off



%Part 2: defining f(t) and F(s)
%set functions symbolically
%syms f(t)
%syms df(t)
%define f(t) where f0 is 50 hz
%f(t) = 50 - ((RP)/(DR+K))(1 + a lot of stuff (referenced the older paper, copy values) );
%run f(t) to acquire a set of data in respect to time

%define df(t) as the function of f(t) and unit step function 
%(a substitutes for tau and will be assigned values)
%df(t) = f(t)*U(t-a);
%Might be able to replace f(t) with the data f(t)... maybe. doesn't work
%symbolically though.