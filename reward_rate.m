function rr =reward_rate(a ,s, deadlines, threshold, x0, x0dist)

[aRT, aER, ~, ~]= multi_stage_ddm_metrics(a ,s, deadlines, threshold, x0, x0dist);

Tmot =0.3; 
Dp=0; 

rr = (1-aER)./(aRT + Tmot + Dp*aER);


