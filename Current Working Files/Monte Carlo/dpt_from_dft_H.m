%Part 0: Initialization of the experiment
%Initialize Global Variables
Percent_Correction = 40;
Pd = -0.3; %the LoG event
SigIn = [0 0]; %load an empty signal for Power injection (SigIn)

%Initialize Monte Carlo Conditions
resolution = 3;

%Initialize Storage Variables:
results = cell(1,6*resolution);
run_num = 0;

%Initialize variables of the SFR model (network values)
%can be modified in the future to run from input files
R = 0.05;
H = 4.0; %mean of 5.5 with dev 1.17
K = 0.95;
Fh = 0.3;
Tr = 8.0;
D = 1.0;
B = Percent_Correction/100; %Correction Factor Value:

%Initialize Simulink Models
GSFR_model = "GSFR_Individual_Vars";
load_system(GSFR_model); %load the GSFR Model with Workspace Variables

inverse_GSFR_model = "inverse_GSFR_Individual_Vars";
load_system(inverse_GSFR_model); %load the inverse GSFR Model


%Run the Actual Monte Carlo Simulation
%Uses a for loop determined by the global variables
%mean, standard deviation, and resolution, for Inertia (H)
for j = mean - 3*stand_dev : stand_dev/resolution : mean + 3*stand_dev
    %set the changing variable to J
    H = j;
    
    %calculate the extra variables
    wn = sqrt((D*R+K)/(2*H*R*Tr));
    c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;
    %polynomial variables
    a1 = (D*R+K)*2*c*wn;
    b1 = wn^2;
    c1 = (wn^2)*R*Tr;
    d1 = (wn^2)*R;
    x1 = D*R+K;


    %Part 1: calculate GSFR(s) and the BASE CASE
    %run the Simulink Model for GSFR, store output
    output = sim(GSFR_model);
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
    inverse_output = sim(inverse_GSFR_model);

    %Increment run number and add to the Results storage Cell Array
    run_num = run_num+1;
    results{1,run_num} = inverse_output.simout; %store the results for Plotting
end


%Plot Results for Saving
clf; %Clear Plot
hold on
for i = 1 : length(results) %Plot all the results in one Figure
    plot(results{1,i},':');
end
%Plot the Mean result last so it sits on top and is most obvious
plot(results{1,round(length(results)/2)},'k'); 
%Add appropriate title and axis labels
plot_title = sprintf('Envelope P(t) for B = %f, H = %d with standard deviation %d',B,mean,stand_dev);
title(plot_title);
xlabel('Time (s) since LoG Event')
ylabel('Active Power (pu)')
hold off

%Save results to figure file and PNG
file_title = sprintf('B_%d_H_%d.png',Percent_Correction,mean);
saveas(gcf,file_title);
file_title = sprintf('B_%d_H_%d',Percent_Correction,mean);
saveas(gcf,file_title);
