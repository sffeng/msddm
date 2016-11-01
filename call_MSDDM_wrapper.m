% Example script for how to call function MSDDM_wrapper
    % Main inputs to change are vectors of stage onsets (deadlines), drift (a) /noise (s) at each stage, and distribution of starting points at 1st stage
    % For convenience, also currently allows for looping over multiple different *uniform* thresholds across the MS-DDM 
        % (using different thresholds for each stage is theoretically possible, but requires a few modifications to the core functions first)
    % Also options for running simulations to compare results (this will obviously take much longer) and some basic plots
% Outputs:
    % aRT, aER = vectors of analytic expected decision time and error rate, one for each threshold being tested [using multi_stage_ddm_metrics]
    % aCDF_T, aCDF_Y = cell array of analytic CDFs (..._T = RT range, ..._Y = cumul prob), one for each threshold being tested [using multistage_ddm_fpt_dist]
    % simMeanRT, simMeanER, simCDF_T, simCDF_Y = same as above, but using MC simulations to generate each value estimate

% clc
% clear all
close all

for i=1:1
%% General wrapper parameters:

% Run MC simulations for the sake of comparison? (takes longer!)
    % NB: simulations can be made to be more efficient than currently (as with analytics!)
    runSimulations = 1;

% Generate plots? Would include following subplots: 
    % (1a-b) analytic vs. simulated CDF for *final* uniform threshold (if testing multiple)
    % (2a-c) analytic vs. simulated DT/ER across uniform thresholds tested
    doPlots = 1;

% Example model to run (set to 0 to manually set parameters below, otherwise specify example number from options below):
    useExample = 3.0;
%     useExample = 98.0;
%     useExample = 99.0;

    
%% Model parameters:

if ~useExample  % Set parameters manually
    
    % Support of initial condition (at first stage) and its density (i.e., can enter distribution or individual point)
    %     x0 = -0.2;   % starting point(s)
        x0 = 0;   % starting point(s)
        x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic

    % Thresholds to test (can enter vector or individual point)
        % (IMPORTANT: these are thresholds to test for different runs of the multi-stage DDM, each with a uniform threshold across stages
        %  Thresholds for different *stages* would be set in varthresh, and then the thresh variable is ignored
        thresh = [1.0];  % if using variable threshold, this can be ignored
    
    % The following vectors (deadlines,a,s,varthresh) need to be the same length:    
        % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
           deadlines = [0 1.5 2.5]; %% Must always start with zero
        % Vector of drift rates (i.e., a0, a1, a2....)
            a = [0.2 0.4 0.6];
        % Vector of diffusion rates (noise coefficients)
            s = [1.0 1.5 1.0];
        % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
            varthresh = nan;
    
else  % Use preset example parameter sets
    
    switch useExample
        
        case 1.0 %%%  EXAMPLE 1.0: 4 stages with variable drift/noise/threshold (testing one uniform threshold) 
            
            x0 = -0.2;   % starting point(s)
            x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
            thresh = [2.0];  % if using variable threshold, this can be ignored
            
            deadlines = [0 1 2 3]; % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
            a = [0.1 0.2 0.05 0.3]; % Vector of drift rates (i.e., a0, a1, a2....)
            s = [1 1.5 1.25 2]; % Vector of diffusion rates (noise coefficients)
            varthresh = nan; % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
            
            simContinuous = nan;
            contSim = [];

        case 1.1 %%%  EXAMPLE 1.1: Same as Example 1.0, but testing multiple uniform thresholds  
            
            x0 = -0.2;   % starting point(s)
            x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
            thresh = [0.5:0.5:4.0];  % if using variable threshold, this can be ignored
            
            deadlines = [0 1 2 3]; % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
            a = [0.1 0.2 0.05 0.3]; % Vector of drift rates (i.e., a0, a1, a2....)
            s = [1 1.5 1.25 2]; % Vector of diffusion rates (noise coefficients)
            varthresh = nan; % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
            
            simContinuous = nan;
            contSim = [];

        case 2.0 %%%  EXAMPLE 2.0: 30 stages with switching drift (e.g., shifting attention; Krajbich et al., 2010) 
            
            x0 = 0;   % starting point(s)
            x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
            thresh = [2.0];  % if using variable threshold, this can be ignored
            
            deadlines = sort([0 randsample(0.05:0.05:10,29,false)]); % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
            a = repmat([1.0,-0.75],1,length(deadlines)/2); % Vector of drift rates (i.e., a0, a1, a2....)
            s = ones(1,length(a)); % Vector of diffusion rates (noise coefficients)
            varthresh = nan; % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
 
            simContinuous = nan;
            contSim = [];
         
        case 3.0 %%%  EXAMPLE 3.0: gradually changing drift (e.g., shrinking spotlight; White et al., 2011) 
           
            x0 = 0;   % starting point(s)
            x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
            thresh = [2.0];  % if using variable threshold, this can be ignored
            
            numDiscreteSteps = 20
%             numDiscreteSteps = 100
            
            deadlines = linspace(0,5,numDiscreteSteps); % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
            a = linspace(-0.2,0.3,numDiscreteSteps); % Vector of drift rates (i.e., a0, a1, a2....)
            s = 1.0*ones(1,numDiscreteSteps); % Vector of diffusion rates (noise coefficients)
            varthresh = nan; % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
 
            simContinuous = 1; % Use continuous (rather than stage-wise) simulation for comparison
            step=0.005;
            contSim.a = [linspace(-0.2,0.3,5/step),0.3*ones(1,45/step)];
            contSim.s = [1.0*ones(1,50/step)];
            contSim.z = [2.0*ones(1,50/step)];

        case 4.0 %%%  EXAMPLE 4.0: collapsing to a lower bounds  
            
            x0 = 0;   % starting point(s)
            x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
            thresh = [2.0];  % if using variable threshold, this can be ignored
            
            deadlines = linspace(0,5,20); % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
            a = 0.15*ones(1,20); % Vector of drift rates (i.e., a0, a1, a2....)
            s = 1.0*ones(1,20); % Vector of diffusion rates (noise coefficients)
            varthresh = linspace(3.0,1e-6,20); % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
        
            simContinuous = 1; % Use continuous (rather than stage-wise) simulation for comparison
            step=0.005;
            contSim.a = [0.15*ones(1,10/step)];
            contSim.s = [1.0*ones(1,10/step)];
            contSim.z = [linspace(3.0,0,5/step),0*ones(1,5/step)];
        
        case 98.0 %%% Single stage DDM
            
            a = [0.2394    0.2311    0.2394    0.2311    0.2394]
            s = [2.3818    2.3818    2.3818    2.3818    2.3818]
            thresh = [1.0]
            x0 =  0
            x0dist =  1
            deadlines =  [0    0.3439    0.8959    1.8638    4.0257]
%             tfinal = 5
%             
%             
%             
%             x0 = 0.0;   % starting point(s)
%             x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
%             thresh = [2.0];  % if using variable threshold, this can be ignored
%             
%             deadlines = [0]; % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
%             a = [1.0]; % Vector of drift rates (i.e., a0, a1, a2....)
%             s = [1]; % Vector of diffusion rates (noise coefficients)
            varthresh = nan; % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
            
            simContinuous = nan;
            contSim = [];

        case 99.0 %%% Single stage DDM
            
            x0 = 0.0;   % starting point(s)
            x0dist = 1;  % probability distribution of starting point(s) -- 1 if deterministic
            thresh = [2.0];  % if using variable threshold, this can be ignored
            
            deadlines = [0]; % Times (secs) at which  each stage starts (i.e., t0, t1, t2....)
            a = [1.0]; % Vector of drift rates (i.e., a0, a1, a2....)
            s = [1]; % Vector of diffusion rates (noise coefficients)
            varthresh = nan; % Vector of thresholds (IF varying by stage) - set to NaN if using uniform thresh
            
            simContinuous = nan;
            contSim = [];
            
                        
        otherwise
            error([num2str(useExample),' is not a valid example number!']);
    end
end

%% Generate results:

% Current error-check to avoid threshold values of exactly zero:
    varthresh(abs(varthresh)<1e-15) = sign(varthresh(abs(varthresh)<1e-15))*1e-15;

% Forcing thresh to be ignored if varying thresholds across stages:    
if ~isnan(varthresh)
    thresh = []; % this extra step isn't necessary, just a reminder
end
    
% Run MSDDM wrapper, collect outputs:
    if runSimulations
        [aRT, aER, aRT_plus, aRT_minus, aCDF_T, aCDF_Y, aCDF_Y_plus, aCDF_Y_minus, ...
            simMeanRT, simMeanER,simMeanRT_plus,simMeanRT_minus, simCDF_T, simCDF_Y, simCDF_Y_plus, simCDF_Y_minus] = ...
            MSDDM_wrapper(a,s,varthresh,deadlines,thresh,x0,x0dist,runSimulations,doPlots,simContinuous,contSim);
    else
        [aRT, aER, aRT_plus, aRT_minus, aCDF_T, aCDF_Y, aCDF_Y_plus, aCDF_Y_minus] = ...
            MSDDM_wrapper(a,s,varthresh,deadlines,thresh,x0,x0dist,runSimulations,doPlots);
    end

%%  Convert CDF to PDF (if useful; ):
% NB: Take care when using this for MLE to properly bin to avoid discontinuity issues at stage transitions:
    for i = 1:length(aCDF_T)
        [tArrayPDF_plus{i},finalPDF_plus{i}] = cdf_to_pdf(aCDF_T{i},aCDF_Y_plus{i});
        [tArrayPDF_minus{i},finalPDF_minus{i}] = cdf_to_pdf(aCDF_T{i},aCDF_Y_minus{i});
    end

    
    dt = 0.01;
    interpCDF_A = interp1(aCDF_T{1},aCDF_Y{1},0:dt:25,'nearest','extrap');
    interpCDF_S = interp1(simCDF_T{1},simCDF_Y{1},0:dt:25,'nearest','extrap');
%     nansum(abs(simCDF_Y_minus{1}-aCDF_Y_minus{1}))*abs(diff(aCDF_T{1}(1:2)))
    nansum(abs(interpCDF_S-interpCDF_A))*dt
    
%     XinterpCDF_S{length(XinterpCDF_S)+1} = interpCDF_S;
    
end



%%
% errorSum = []
% for iii=1:length(XinterpCDF_S)
%     errorSum(iii,1) = nansum(abs(XinterpCDF_S{iii}-XinterpCDF_A_2))*dt;
%     errorSum(iii,2) = nansum(abs(XinterpCDF_S{iii}-XinterpCDF_A_5))*dt;
%     errorSum(iii,3) = nansum(abs(XinterpCDF_S{iii}-XinterpCDF_A_10))*dt;
%     errorSum(iii,4) = nansum(abs(XinterpCDF_S{iii}-XinterpCDF_A_20))*dt;
%     errorSum(iii,5) = nansum(abs(XinterpCDF_S{iii}-XinterpCDF_A_50))*dt;
%     errorSum(iii,6) = nansum(abs(XinterpCDF_S{iii}-XinterpCDF_A_100))*dt;
% end
% % 
% % 
% figure
% errorbar(1:6,mean(errorSum),[AS_nanSE(errorSum)],'k');
% set(gca,'TickDir','out');
% set(gca,'Box','off');
% ylabel('CDF error');
% xlabel('# stages');
% set(gca,'XTick',1:6);
% set(gca,'XTickLabel',{'2','5','10','20','50','100'});
% % 
% saveas(gcf,['tmpCDFerror_041216.ai'],'ai')
% saveas(gcf,['tmpCDFerror_041216.pdf'],'pdf')
