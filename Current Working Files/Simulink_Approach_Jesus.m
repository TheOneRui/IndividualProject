%Part 1: calculate GSFR(s)
model1 = "GSFR_with_P_sigIn";
open_system(model1);

%create initial Pd, load an empty signal for Power injection (SigIn)
Pd = -0.3;
SigIn = [0 0];

%run the Simulink Model for GSFR, store output
GSFR_output = sim(model1);
%the Data is hidden in GSFR_output's sub section of "simout"
plot(GSFR_output.simout);



% Injection_File = load("Correction_half");
% SigIn = Injection_File.SigIn;