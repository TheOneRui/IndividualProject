mdl = "function1";
open_system(mdl);


Pds = [-0.05 -0.10 -0.15 -0.20 -0.25 -0.3];

%assign the block (blk) to the one changing
blk = mdl + "/Pd";

%create initial Pd
Pd = 0;

%sets the value of constant block to the parameter Pd
set_param(blk,"Value","Pd");

for k = length(Pds):-1:1
    simIn(k) = Simulink.SimulationInput(mdl);
    simIn(k) = setVariable(simIn(k),"Pd",Pds(k));
end

out = sim(simIn,"UseFastRestart","on");
y1=out(1).yout;
y2=out(2).yout;
y3=out(3).yout;
y4=out(4).yout;
y5=out(5).yout;
y6=out(6).yout;

%this is ugly, and I would prefer to have preprepared the data for plotting
hold on
plot(y2{1}.Values)
plot(y3{1}.Values)
plot(y4{1}.Values)
plot(y5{1}.Values)
plot(y6{1}.Values)
hold off