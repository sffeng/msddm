clc
clear all
close all



% Discretization step for simulating DDM
step=0.001;


% Support of initial condition and its density

x0=0; x0dist=1;

% No of Monte Carlo runs

realizations=10000;

%A=0.15; 

thresh=1;

S=1;

% Initialization

RT=zeros(1,realizations);
ER=zeros(1,realizations);

% Simulate the multistage DDM

numIters=1;
for iteri=1:numIters
    for N=1:realizations
        
        t=0;
        
        x(1)=x0;
        
        l=1;
        
        stop=0;
        
        
        while stop==0
            
            %stage=find(deadlines<=t,1,'last');
            
            A = min(0.3, -0.2 +0.1*t);
            
            x(l+1)= x(l) +A*step + S*randn*sqrt(step);
            t=t+step;
            l=l+1;
            
            %thresh= 3- 0.6*t;
            
            if (x(l)>=thresh || x(l)<=-thresh)
                stop=1;
                RT(N)=t;
                ER(N)=(x(l)<=-thresh);
            end
        end
    end
    
    RTplus=RT(ER==0);
    
    RTminus=RT(ER==1);
    
    
    
    [Nhist{iteri},T{iteri}] = hist(RT,100);
    
    [Np, Tp]=hist(RTplus,100);
    
    [Nm, Tm]=hist(RTminus,100);
    
    tfinal(iteri)= T{iteri}(end);
    
    T_inter{iteri}=0:0.01:tfinal(iteri);
    Y_sim{iteri}=cumsum(Nhist{iteri})/realizations;    
    sim_inter{iteri}=interp1(T{iteri}, Y_sim{iteri}, T_inter{iteri},'nearest','extrap');
end


% stage_vec=2:20; %5:5:30;%10:10:100;
stage_vec=2:15; %5:5:30;%10:10:100;

for ii=1:length(stage_vec)
    
    stages= stage_vec(ii);
    %stages=50;
    
    % Vector of drift rates
    
    a= linspace(-0.2, 0.3, stages);
    %0.15*ones(1,stages);
    
    % Vector of diffusion rates
    
    s=1.0*ones(1,stages);
    
    % Thresholds
    
    threshold= 1*ones(1,stages);
    %linspace(3.0,0.01,stages);
    
    % Times at which  each stage starts
    
    deadlines=linspace(0,5,stages);
    
    
    [T_anal,Y_anal,Y_plus_anal,Y_minus_anal]=multistage_ddm_fpt_dist(a,s,threshold,x0,x0dist,deadlines,tfinal);
        
    
    %cdf_er(ii,:)= ana_inter-sim_inter;
    for iteri=1:numIters
        ana_inter=interp1(T_anal, Y_anal, T_inter{iteri},'nearest','extrap');
        cdf_er(iteri,ii)= sqrt(sum((ana_inter-sim_inter{iteri}).^2)*0.01);
    end
end



close all;

figure
plot(stage_vec, abs(cdf_er),'k','linewidth',2)
% ms = mean(cdf_er);
% ses = AS_nanSE(cdf_er);
% % errorbar(ms,ses,'k');
set(gca,'TickDir','out');
set(gca,'Box','off');
ylabel('CDF Error (MSE)');
xlabel('# Stages');
set(gca,'XTick',1:length(stage_vec));
xlim([0.5 14.5]);
set(gca,'XTickLabel',{'2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
% % % set(gca,'XTickLabel',{'2','5','10','20','50','100'});
% saveas(gcf,['tmpCDFerror_s15_MSE_041316.ai'],'ai')
% saveas(gcf,['tmpCDFerror_s15_MSE_041316.pdf'],'pdf')

%axis([0 10 0 0.05 ])
% 
% figure
% 
% 
% plot(T, Y_sim,'k','linewidth',2);
% 
% hold on
% 
% plot(T_anal,Y_anal, 'r--','linewidth',2)
% 
% 
% xlabel('Decision Time'); ylabel('CDF')
% 
% 
% legend('Simulation', 'Analytic','location','best')
% 
% figure
% 
% plot(Tp, cumsum(Np)/realizations,'k','linewidth',2);
% 
% hold on
% 
% plot(T_anal,Y_plus_anal, 'r--','linewidth',2)
% 
% 
% xlabel('Correct Decision Time'); ylabel('CDF')
% 
% 
% legend('Simulation', 'Analytic','location','best')
% 
% 
% figure
% 
% plot(Tm, cumsum(Nm)/realizations,'k','linewidth',2);
% 
% hold on
% 
% plot(T_anal,Y_minus_anal, 'r--','linewidth',2)
% 
% 
% xlabel('Erroneous Decision Time'); ylabel('CDF')
% 
% 
% legend('Simulation', 'Analytic','location','best')

