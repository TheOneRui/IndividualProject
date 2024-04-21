%---------------Initialize the Experiment and Models------------------
GSFR_model = "GSFR_Individual_Vars";
load_system(GSFR_model); %load the GSFR Model with Workspace Variables
inverse_GSFR_model = "inverse_GSFR_Individual_Vars";
load_system(inverse_GSFR_model); %load the inverse GSFR Model
%---------------Initialize the Experiment and Models------------------




%-----------------Initialize Experiment Variables---------------------
Total_Num_Runs = 10;
run_num = 2;
Percent_Correction = 100;
Pd = -0.1; %the LoG event
SigIn = [0 0]; %load an empty signal for Power injection (SigIn)

%Initialize Storage Variables:
results = cell(2,Total_Num_Runs+2);
run_num = 2;
%-----------------Initialize Experiment Variables---------------------





%-------------------Initialize Progress Tracker-----------------------
clf; %Clear Plot
progress = waitbar(0,'running simulation');
%-------------------Initialize Progress Tracker-----------------------




%----------------Part 1: Estimated GSFR Response----------------------
%Initialize Estimated Variables (See Anderson's Paper)
R = 0.05;
H = 4;
K = 0.95;
Fh = 0.3;
Tr = 8.0;
D = 1.0;
B = Percent_Correction/100; %Correction Factor Value:
%calculate the extra variables
wn = sqrt((D*R+K)/(2*H*R*Tr));
c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;
%polynomial variables
a1 = (D*R+K)*2*c*wn;
b1 = wn^2;
c1 = (wn^2)*R*Tr;
d1 = (wn^2)*R;
x1 = D*R+K;
%Run the Model
output = sim(GSFR_model);
GSFR_output = output.simout;
results{1,1} = GSFR_output; %store the results for Plotting Later
results{2,1} = [R H K Fh Tr D];
%----------------Part 1: Estimated GSFR Response----------------------




%---------Part 2: Finding Injection for Estimated Values--------------
%Find Nadir, Establish Tau and Tau2
Nadir = min(GSFR_output);
settling_Frequency = 1 + ((R*Pd)/(D*R+K));
Dtr = settling_Frequency - Nadir;
New_Nadir = B*Dtr + Nadir;

%Finding tau
tau = GSFR_output.Time(find(GSFR_output.Data < Dtr*B+Nadir,1,'first'));
tau_Index = find(GSFR_output.Data < Dtr*B+Nadir,1,'first');
%Finding tau2
tau2 = GSFR_output.Time(find(GSFR_output.Data < Dtr*B+Nadir,1,'last'));
tau2_Index = find(GSFR_output.Data < Dtr*B+Nadir,1,'last');

%Manually applying Unit Step Function
dft = GSFR_output;
for i = 1 : length(dft.Data)
    if(i >= tau_Index && i <= tau2_Index)
        dft.Data(i) = New_Nadir - GSFR_output.Data(i);
    else
        dft.Data(i) = 0;
    end
end

%Run the Delta f(t) through Inverse GSFR for injection p(t)
inverse_output = sim(inverse_GSFR_model);
SigIn = inverse_output.simout; %store the results for Monte Carlo

%Run Injection through GSFR to get the Base Case response with Injection
output = sim(GSFR_model);
GSFR_output = output.simout;
results{1,2} = GSFR_output; %store the results for Plotting later
results{2,2} = [R H K Fh Tr D];
%---------Part 2: Finding Injection for Estimated Values--------------




%----------------------Part 3: the Experiment-------------------------
%Insert the Power Injection from Part 2 into a GSFR with Random Values
%Determined by a Normal Distribution for each Variable or System condition
for i = 1 : Total_Num_Runs
    %Increment Run Number
    run_num = run_num + 1;
    %Randomize the Variables by Normal Distribution
    R = (0.05/3)*randn()+0.075;
    H = 1.17*randn()+5.5; %mean of 5.5 with dev 1.17
    K = 0.95;
    Fh = 0.034*randn()+0.25;
    Tr = 0.67*randn()+9;
    D = 1.0; 
    wn = sqrt((D*R+K)/(2*H*R*Tr));
    c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;
    
    %Run the GSFR model with new Variable values
    output = sim(GSFR_model);
    GSFR_output = output.simout;
    results{1,run_num} = GSFR_output; %store the results for Plotting later
    results{2,run_num} = [R H K Fh Tr D];
    
    
    %display progress
    waitbar((1/Total_Num_Runs)*(run_num/2),progress);
end
%----------------------Part 3: the Experiment-------------------------




%------------------Part #: Plotting the results-----------------------
%Update progress bar
waitbar(0.5,progress,'Simulations Complete -> Plotting results');
%Plot Results for Saving
hold on
for i = 3 : length(results) %Plot all the results in one Figure
    plot(results{1,i},':');
    waitbar(0.5+(1/Total_Num_Runs)*((i-2)/2),progress);
end
close(progress);
%Plot the Mean result last so it sits on top and is most obvious
plot(results{1,1},'k--'); 
plot(results{1,2},'k'); 
%Add appropriate title and axis labels
plot_title = sprintf('Envelope f(t) for B = %f with %d runs',B,Total_Num_Runs);
title(plot_title);
xlabel('Time (s) since LoG Event')
ylabel('Frequency (pu)')
hold off
%------------------Part #: Plotting the results-----------------------







%--------------------------Saving Results-----------------------------
folder_title = sprintf('results_B_%d_Runs_%d',Percent_Correction,Total_Num_Runs);
mkdir (folder_title);

%Save results to figure file, PNG, and raw Cell Array
file_title = sprintf('B_%d_Runs_%d',Percent_Correction,Total_Num_Runs);
saveas(gcf,file_title);
movefile([file_title '.fig'], folder_title);

file_title1 = sprintf('B_%d_Runs_%d.png',Percent_Correction,Total_Num_Runs);
saveas(gcf,file_title1);
movefile(file_title1, folder_title);

file_title2 = sprintf('B_%d_Runs_%d.mat',Percent_Correction,Total_Num_Runs);
save(file_title2, "results");
movefile(file_title2, folder_title);
%--------------------------Saving Results-----------------------------