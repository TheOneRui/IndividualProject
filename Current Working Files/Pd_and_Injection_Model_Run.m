%Part 1: calculate GSFR(s)
model1 = "GSFR_with_P_sigIn";
open_system(model1);

%create initial Pd
Pd = -0.3;
Injection_File = load("correction_point_8");
SigIn = Injection_File.y11;

out = sim(model1);



plot(out.simout);


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