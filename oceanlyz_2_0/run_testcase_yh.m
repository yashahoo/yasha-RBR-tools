    %Change current working directory to OCEANLYZ folder
    %Assume OCEANLYZ files are in 'C:\oceanlyz_matlab' folder
    cd('/Users/00068592/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools_v5/oceanlyz_2_0')

    %Create OCEANLYZ object
    %clear ocn %Optional
    ocn=oceanlyz;
    
    %Read data
    %Assume data file is named 'waterpressure_5burst.csv' and is stored in 'C:\oceanlyz_matlab\Sample_Data'
    current_folder=pwd;                  %Current (OCEANLYZ) path
    cd('Sample_Data/') %Change current path to Sample_Data folder
    water_pressure=importdata('waterpressure_5burst.csv'); %Load data
    cd(current_folder)                   %Change current path to OCEANLYZ folder
    
    %Input parameters
    ocn.data= water_pressure;
    ocn.InputType='pressure';
    ocn.OutputType='wave+waterlevel';
    ocn.AnalysisMethod='spectral';
    ocn.n_burst=5;
    ocn.burst_duration=1024;
    ocn.fs=10;
    ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.heightfrombed=0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.dispout='yes';
    ocn.Rho=1024;                     %Seawater density (Varies)

    %Run OCEANLYZ
    ocn.runoceanlyz()

    %Plot peak wave period (Tp)
    plot(ocn.wave.Tp(1,:))

    
    %% run with zero cross
    
        %Assume OCEANLYZ files are in 'C:\oceanlyz_matlab' folder
    cd('/Users/00068592/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools_v5/oceanlyz_2_0')

    %Create OCEANLYZ object
    %clear ocn %Optional
    ocn=oceanlyz;
    
    %Read data
    %Assume data file is named 'waterpressure_5burst.csv' and is stored in 'C:\oceanlyz_matlab\Sample_Data'
    current_folder=pwd;                  %Current (OCEANLYZ) path
    cd('Sample_Data/') %Change current path to Sample_Data folder
    water_pressure=importdata('waterpressure_5burst.csv'); %Load data
    cd(current_folder)                   %Change current path to OCEANLYZ folder
    
    %Input parameters
    ocn.data= water_pressure;
    ocn.InputType='pressure';
    ocn.OutputType='wave+waterlevel';
    ocn.AnalysisMethod='zerocross';
    ocn.n_burst=5;
    ocn.burst_duration=1024;
    ocn.fs=10;
    ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.heightfrombed=0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.dispout='yes';
    ocn.Rho=1024;                     %Seawater density (Varies)

    %Run OCEANLYZ
    ocn.runoceanlyz()

    %Plot peak wave period (Tp)
    figure;
    plot(ocn.wave.Ts(1,:))
    
    %% run with spectral separating sea swell
    
        %Assume OCEANLYZ files are in 'C:\oceanlyz_matlab' folder
    cd('/Users/00068592/Documents/MATLAB/m-files/yashatools/yasha_RBR_tools_v5/oceanlyz_2_0')

    %Create OCEANLYZ object
    %clear ocn %Optional
    ocn=oceanlyz;
    
    %Read data
    %Assume data file is named 'waterpressure_5burst.csv' and is stored in 'C:\oceanlyz_matlab\Sample_Data'
    current_folder=pwd;                  %Current (OCEANLYZ) path
    cd('Sample_Data/') %Change current path to Sample_Data folder
    water_pressure=importdata('waterpressure_5burst.csv'); %Load data
    cd(current_folder)                   %Change current path to OCEANLYZ folder
    
    %Input parameters
    ocn.data= water_pressure;
    ocn.InputType='pressure';
    ocn.OutputType='wave+waterlevel';
    ocn.AnalysisMethod='spectral';
    ocn.n_burst=5;
    ocn.burst_duration=1024;
    ocn.fs=10;
    ocn.fmin=0.05;                    %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmax=ocn.fs/2;                %Only required if ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorrCalcMethod='auto';   %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.Kpafterfmaxpcorr='constant';  %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fminpcorr=0.15;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.fmaxpcorr=0.55;               %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.heightfrombed=0.05;           %Only required if ocn.InputType='pressure' and ocn.AnalysisMethod='spectral'
    ocn.dispout='yes';
    ocn.Rho=1024;                     %Seawater density (Varies)
    ocn.SeparateSeaSwell='yes'

    %Run OCEANLYZ
    ocn.runoceanlyz()

    %Plot peak wave period (Tp)
    figure;
    plot(ocn.wave.Tp(1,:))