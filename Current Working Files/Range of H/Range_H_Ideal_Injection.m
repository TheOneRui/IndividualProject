%---------------Initialize the Experiment and Models------------------
GSFR_model = "GSFR_Individual_Vars";
load_system(GSFR_model); %load the GSFR Model with Workspace Variables
inverse_GSFR_model = "inverse_GSFR_Individual_Vars";
load_system(inverse_GSFR_model); %load the inverse GSFR Model
%---------------Initialize the Experiment and Models------------------




%-----------------Initialize Experiment Variables---------------------
Percent_Correction = 100;
Pd = -0.1; %the LoG event
SigIn = [0 0]; %load an empty signal for Power injection (SigIn)

%Initialize Storage Variables:
results = cell(5,12);
%-----------------Initialize Experiment Variables---------------------





%-------------------Initialize Progress Tracker-----------------------
clf; %Clear Plot
progress = waitbar(0,'running simulation');
%-------------------Initialize Progress Tracker-----------------------


for j = 1: 12

    %----------------Part 1: Estimated GSFR Response----------------------
    %Initialize Estimated Variables (See Anderson's Paper)
    SigIn = [0 0];
    H = j;
    R = 0.04;
    K = 1;
    Fh = 0.4;
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
    results{1,j} = [H R K Fh Tr D];
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
    %---------Part 2: Finding Injection for Estimated Values--------------




    %----------------------Part 3: the Experiment-------------------------
    %Using Injection, check the lower bounds for H's innacuracy
    %store the last value that does not violate Statutory Limits
    for i = 1 : 99 
        H = j-((i/100)*j); 
        wn = sqrt((D*R+K)/(2*H*R*Tr));
        c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;

        %Run the GSFR model with new Variable values
        output = sim(GSFR_model);
        GSFR_output = output.simout;
        
        if min(GSFR_output) > 0.99
            results{2,j} = i;
        else
            results{4,j} = GSFR_output;
            break;
        end
    end
    
    
    %Using Injection, check the upper bounds for H's innacuracy
    %store the last value that does not violate Statutory Limits
    for i = 1 : 99 
        H = j+((i/100)*j); 
        wn = sqrt((D*R+K)/(2*H*R*Tr));
        c = (((2*H*R)+(((D*R)+(K*Fh))*Tr))/(2*(D*R+K)))*wn;

        %Run the GSFR model with new Variable values
        output = sim(GSFR_model);
        GSFR_output = output.simout;
        
        if max(GSFR_output) < 1.01
            results{3,j} = i;
        else
            results{5,j} = GSFR_output;
            break;
        end
    end
    %----------------------Part 3: the Experiment-------------------------
    %display progress
    waitbar((1/12)*(j),progress);
end

hold on
for i = 1 : length(results) %Plot all the results in one Figure
    plot(results{4,i});
    plot(results{5,i});
end
yline(0.99, 'k', 'DisplayName', 'lower Statutory Limit');
yline(1.01, 'k', 'DisplayName', 'upper Statutory Limit');
hold off

%--------------------------Saving Results-----------------------------
folder_title = sprintf('results_ideal_R_%f_K_%f_Fh_%f',R,K,Fh);
mkdir (folder_title);

%Save results to figure file, PNG, and raw Cell Array
file_title2 = sprintf('boundaries.mat');
save(file_title2, "results");
movefile(file_title2, folder_title);
%--------------------------Saving Results-----------------------------

close(progress);