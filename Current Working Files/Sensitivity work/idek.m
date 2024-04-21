
%Initialize the Expirement parameters
Pd = -0.1; 
Total_Num_Runs = 20;
Accuracy = 15; %Tolerance or Uncertainty

%Folder Creation for Result Storage
folder_title = sprintf('results');
mkdir (folder_title);


%Progress Bar for the Longer Monte Carlo Sims
progress = waitbar(0,'running simulation');

%Initialize the SFR Simulink Models
GSFR_model = "GSFR_Individual_Vars_sensitivity";
load_system(GSFR_model);   


%Initialize Storage Variables:
results = zeros(Total_Num_Runs,8);

%Run the Monte Carlo Test
for i = 1: Total_Num_Runs
    %Randomize the Variables by Normal Distribution
    D = (1*(Accuracy/100))*randn()+1; %1
    R = (0.05*(Accuracy/100))*randn()+0.05; %2
    H = (4*(Accuracy/100))*randn()+4; %3
    K = (0.95*(Accuracy/100))*randn()+0.95; %4
    Fh = (0.3*(Accuracy/100))*randn()+0.3; %5
    Tr = (8*(Accuracy/100))*randn()+8; %6

    wn = sqrt((D*R+K)/(2*H*R*Tr));
    c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*sqrt((D*R+K)/(2*H*R*Tr));

    %Run the GSFR model with new Variable values
    output = sim(GSFR_model);
    GSFR_output = output.simout;
    
    %calculate and store the minimum value
    results(i,7) = min(GSFR_output);
    results(i,8) = GSFR_output.Data(end, :);
    
    %store each variable
    results(i,1) = D;
    results(i,2) = R;
    results(i,3) = H;
    results(i,4) = K;
    results(i,5) = Fh;
    results(i,6) = Tr;

    %display progress
    waitbar((i/Total_Num_Runs),progress);
end
    %Update progress bar
    waitbar(1,progress,'Simulations Complete -> Plotting results');
    %Plot Results for Saving
    close(progress);


    %Save results to figure file, PNG, and raw Cell Array
    file_title2 = sprintf('results.mat');
    save(file_title2, "results");
    movefile(file_title2, folder_title);