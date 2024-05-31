%Initialization of the experiment
%Correction Factor Value:
B = 1;
%Pd or P step loss in the case of an LoG event
Pd = -0.1;

%Part 1: calculate GSFR(s) and the BASE CASE
model1 = "GSFR_with_P_sigIn";
open_system(model1);
%load an empty signal for Power injection (SigIn)
SigIn = [0 0];

%run the Simulink Model for GSFR, store output
output = sim(model1);
GSFR_output = output.simout;
%the Data is hidden in GSFR_output's sub section of "simout" (Plotted for
%ease of visibility)
%plot(GSFR_output);



%Part 2: establishing the signals needed
%Need to establish Tau and Tau2, using the value of B (correction factor) 
%and the Dtr, Settling frequency, and frequency nadir.
%can then manually do the unit step function work with a for loop
Nadir = min(GSFR_output);
settling_Frequency = GSFR_output.Data(end, :);
Dtr = settling_Frequency - Nadir;
New_Nadir = B*Dtr + Nadir;

%Finding tau
tau = GSFR_output.Time(find(GSFR_output.Data < Dtr*B+Nadir,1,'first'));
tau_Index = find(GSFR_output.Data < Dtr*B+Nadir,1,'first');
%Finding tau2
tau2 = GSFR_output.Time(find(GSFR_output.Data < Dtr*B+Nadir,1,'last'));
tau2_Index = find(GSFR_output.Data < Dtr*B+Nadir,1,'last');


%creating the compensating injection frequency
dft = GSFR_output;
for i = 1 : length(dft.Data)
    if(i >= tau_Index && i <= tau2_Index)
        dft.Data(i) = New_Nadir - GSFR_output.Data(i);
    else
        dft.Data(i) = 0;
    end
end


%Part 3: running this through a Simulink Model which does the Laplace
%transform and Inversve Laplace transform using time series dft and GSFR
model2 = "inverse_GSFR";
open_system(model2);
pt_output = sim(model2);


hold on
plot(dft) %Just proving it works and looks right
xlim([0 20]);
ylim([0 0.2]);
plot(pt_output.simout);
hold off


% Injection_File = load("Correction_half");
% SigIn = Injection_File.SigIn;