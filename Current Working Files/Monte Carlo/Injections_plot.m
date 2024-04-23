%Part 0: Initialization of the experiment
%Initialize Global Variables
Percent_Correction = 80;
Pd = -0.3; %the LoG event
SigIn = [0 0]; %load an empty signal for Power injection (SigIn)

%Initialize Simulink Models
GSFR_model = "GSFR_Individual_Vars";
load_system(GSFR_model); %load the GSFR Model with Workspace Variables
inverse_GSFR_model = "inverse_GSFR_Individual_Vars";
load_system(inverse_GSFR_model); %load the inverse GSFR Model

%Progress Bar for the Longer Monte Carlo Sims
clf; %Clear Plot

%Part 1: Running The Assumed Values GSFR model
%Initialize variables of the ORIGINAL SFR model (assumed network values)
R = 0.075;
H = 5.5;
K = 0.95;
Fh = 0.25;
Tr = 9.0;
D = 1.0;
B = Percent_Correction/100; %Correction Factor Value:


%Part 1.1: calculate GSFR(s) and the BASE CASE
%calculate the extra variables
wn = sqrt((D*R+K)/(2*H*R*Tr));
c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;
%polynomial variables
a1 = (D*R+K)*2*c*wn;
b1 = wn^2;
c1 = (wn^2)*R*Tr;
d1 = (wn^2)*R;
x1 = D*R+K;

%run the Simulink Model for GSFR, store output
output = sim(GSFR_model);
GSFR_output = output.simout;


%Part 2: establishing the signals needed
%Need to establish Tau and Tau2, using the value of B (correction factor) 
%and the Dtr, Settling frequency, and frequency nadir.
%can then manually do the unit step function work with a for loop
Nadir = min(GSFR_output);
settling_Frequency = 1 + ((R*Pd)/(D*R+K));
Dtr = settling_Frequency - Nadir;
New_Nadir = B*Dtr + Nadir;

%Finding tau
tau_Index = find(GSFR_output.Data < Dtr*B+Nadir,1,'first');
%Finding tau2
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


%Part 2: running Delta f(t) through Inverse GSFR for injection p(t)
inverse_output = sim(inverse_GSFR_model);
Ideal_Injection = inverse_output.simout; %store the results for Monte Carlo


%Manipulate SigIn to use a "Simplified" Triangular Injection
injection_Max = max(Ideal_Injection);
injection_Index = find(Ideal_Injection.Data >= injection_Max,1,'first');
injection_Final = Ideal_Injection.Data(find(Ideal_Injection.Data >= 0,1,'last'));
final_Index = length(Ideal_Injection.Data);
%find the "Half" point, to not include second triangle
for i = injection_Index : length(Ideal_Injection.Data)-1
    if (Ideal_Injection.Data(i)-Ideal_Injection.Data(i+1)<0)
        final_Index = i;
        if (Ideal_Injection.Data(i)>0)
            injection_Final = Ideal_Injection.Data(i);
        end
        break
    end
end
injection_Steps = (injection_Max - injection_Final) /(final_Index - injection_Index);

%copy injection for overiding
Halved_injection = Ideal_Injection;

for i = injection_Index : length(Ideal_Injection.Data)
    value = injection_Max - ((i-injection_Index)*injection_Steps);
    if(value>=0)
        Halved_injection.Data(i) = value;
    else
        Halved_injection.Data(i) = 0;
    end
end

%Manipulate SigIn to use a "Simplified" Triangular Injection
injection_Max = max(Ideal_Injection);
injection_Index = find(Ideal_Injection.Data >= injection_Max,1,'first');
injection_Final = Ideal_Injection.Data(length(Ideal_Injection.Data));
injection_Steps = (injection_Max - injection_Final) /(length(Ideal_Injection.Data) - injection_Index);

simplified_injection = Ideal_Injection;

for i = injection_Index : length(Ideal_Injection.Data)
    simplified_injection.Data(i) = injection_Max - ((i-injection_Index)*injection_Steps);
end

%Plot and Save Ideal Injection for comparison
hold on
plot(simplified_injection,'DisplayName','Simplified Injection');
plot(Ideal_Injection,'DisplayName','Ideal Injection');
plot(Halved_injection,'DisplayName','Halved Injection');
legend
plot_title = sprintf('Power Injection for B = %f',B);
title(plot_title);
xlabel('Time (s) since LoG Event')
ylabel('Active Power (pu)')
hold off
file_title1 = sprintf('Power_Injection_B_%d.png',Percent_Correction);
saveas(gcf,file_title1);
