%Initialize the Expirement parameters
Percent_Correction = 100;
Pd = -0.3; 
SigIn = [0 0]; %load an empty signal for Power injection (SigIn)
Total_Num_Runs = 20;
Accuracy = 5; %Tolerance or Uncertainty
Var_Mean = 2; %Expected Value of Inertia Within the Grid


%Initialize the SFR Simulink Models
GSFR_model = "GSFR_Individual_Vars";
load_system(GSFR_model);   
inverse_GSFR_model = "inverse_GSFR_Individual_Vars";
load_system(inverse_GSFR_model); 


%Initialize Storage Variables:
results = cell(2,Total_Num_Runs+2);
run_num = 2;


for j = 0 : 2.5: 50
    Accuracy = j;
    run_num = 2;
    
    %Folder Creation for Result Storage
    folder_title = sprintf('results_D_%f_Uncertainty_%f_Pd_%f',Var_Mean,Accuracy,Pd);
    mkdir (folder_title);


    %Progress Bar for the Longer Monte Carlo Sims
    clf; %Clear Plot
    progress = waitbar(0,'running simulation');


    %Part 1: Running The Assumed Values GSFR model
    %Initialize variables of the ORIGINAL SFR model (assumed network values)
    R = 0.075;
    H = 4;
    K = 0.95;
    Fh = 0.3;
    Tr = 9;
    D = Var_Mean;
    B = Percent_Correction/100; %Correction Factor Value:


    %Part 1.1: calculate GSFR(s) and the BASE CASE
    %calculate the extra variables
    wn = sqrt((D*R+K)/(2*H*R*Tr));
    c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;


    %run the Simulink Model for GSFR, store output
    %This is the base case, where all variables are
    %Equal to the mean variable
    output = sim(GSFR_model);
    GSFR_output = output.simout;
    results{1,1} = GSFR_output; %store the results for Plotting Later
    results{2,1} = D;


    %Run the Monte Carlo Test
    for i = 1 : Total_Num_Runs
        %Increment Run Number
        run_num = run_num + 1;
        %Randomize the Variables by Normal Distribution
        D = (Var_Mean*(Accuracy/100))*randn()+Var_Mean; %mean of Var_Mean with standard dev of Accuracy
        while D <= 0
            D = (Var_Mean*(Accuracy/100))*randn()+Var_Mean;
        end
        wn = sqrt((D*R+K)/(2*H*R*Tr));
        c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;

        %Run the GSFR model with new Variable values
        output = sim(GSFR_model);
        GSFR_output = output.simout;
        results{1,run_num} = GSFR_output; %store the results for Plotting later
        results{2,run_num} = D;


        %display progress
        waitbar((1/Total_Num_Runs)*(run_num/2),progress);
    end


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
    plot(results{1,1},'k--','DisplayName', 'Mean'); 
    yline(49.5/50, 'k', 'DisplayName', 'Statutory Limit');
    %Add appropriate title and axis labels
    plot_title = sprintf('Envelope f(t) for D mean = %f with %d percent uncertainty',Var_Mean,Accuracy);
    title(plot_title);
    xlabel('Time (s) since LoG Event')
    ylabel('Frequency (pu)')
    %legend
    hold off


    %Save results to figure file, PNG, and raw Cell Array
    file_title = sprintf('results');
    saveas(gcf,file_title);
    movefile([file_title '.fig'], folder_title);

    file_title1 = sprintf('results.png');
    saveas(gcf,file_title1);
    movefile(file_title1, folder_title);

    file_title2 = sprintf('results.mat');
    save(file_title2, "results");
    movefile(file_title2, folder_title);
end