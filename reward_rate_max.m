
function opt_thresh =reward_rate_max(a ,s, deadlines, x0, x0dist, z1, z2)


% This function implements golden section rule to compute the threshold
% that maximizes the reward rate.  
%Inputs:
% a =vector of drift rates
% s = vector of diffusion rates 
% deadlines = vector of start times of stages
% x0 = initial condition , x0dist = distribution of initial condition
% z1  =lower bound on threshold 
% z2  = upper bound on threshold

% Output
% opt_thresh = optimal threshold



epsilon=0.0001;               % accuracy value
iter= 100;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618

k=0;                            % number of iterations


 % computing x values
 
x1=z1+(1-tau)*(z2-z1);            
x2=z1+tau*(z2-z1);


% computing values in x points

f_x1=-reward_rate(a ,s, deadlines, x1*ones(size(a)), x0, x0dist);                     
f_x2=-reward_rate(a ,s, deadlines, x2*ones(size(a)), x0, x0dist); 



while ((abs(z2-z1)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        z2=x2;
        x2=x1;
        x1=z1+(1-tau)*(z2-z1);
        f_x2=f_x1;
        f_x1=-reward_rate(a ,s, deadlines, x1*ones(size(a)), x0, x0dist); 
        
    else
        z1=x1;
        x1=x2;
        x2=z1+tau*(z2-z1);
        f_x1=f_x2;
        f_x2=-reward_rate(a ,s, deadlines, x2*ones(size(a)), x0, x0dist); 
        
    end
    
end


% chooses minimum point
if(f_x1<f_x2)
    
    opt_thresh=x1;
    
    %max_val=-f_x1;
    
else
    opt_thresh=x2;
    %max_val=-f_x2;
end


function rr =reward_rate(a ,s, deadlines, threshold, x0, x0dist)

[aRT, aER, ~, ~]= multi_stage_ddm_metrics(a ,s, deadlines, threshold, x0, x0dist);

% T_mot = motor time  
% Dp = additional response time after erroneous decision

Tmot =0.3; 
Dp=0; 

rr = (1-aER)./(aRT + Tmot + Dp*aER);




