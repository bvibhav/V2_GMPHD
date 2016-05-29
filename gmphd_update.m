function [ Hyp ] = gmphd_update( Hyp, model, sensorScan)
%GMPHD_UPDATE Summary of this function goes here
%   Detailed explanation goes here

nHyp = numel(Hyp);

% construction of PHD update components
for j = 1:nHyp
    Hyp(j).neta = model.H * Hyp(j).mk;
    Hyp(j).Sk = model.R + model.H * Hyp(j).Pk * model.H';
    Hyp(j).Kk = Hyp(j).Pk * model.H' * pinv(Hyp(j).Sk);
    Hyp(j).Pk = (eye(4) - Hyp(j).Kk * model.H) * Hyp(j).Pk;
end
    
% update
l_count = 0;
for i_obs = 1:numel(sensorScan.xMeas)
    l_count=l_count+1;
    z = [sensorScan.xMeas(i_obs);...
         sensorScan.yMeas(i_obs)];
    w_sum = 0;
    for j = 1:nHyp
        Hyp(l_count*nHyp+j).wk = model.pD*Hyp(j).wk;
        Hyp(l_count*nHyp+j).mk = Hyp(j).mk + Hyp(j).Kk*(z-Hyp(j).neta);
        Hyp(l_count*nHyp+j).Pk = Hyp(j).Pk;
        w_sum = w_sum + Hyp(l_count*nHyp+j).wk;
    end
    for j = 1:nHyp
      Hyp(l_count*nHyp+j).wk = Hyp(l_count*nHyp+j).wk/w_sum;
    end
end
% nHyp = l_count*nHyp + nHyp;
for j = 1:nHyp
    Hyp(j).wk = (1 - model.pD) * Hyp(j).wk;
end

end

