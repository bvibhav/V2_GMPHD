%% Intialize 
clc;
clear;
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultFigureWindowStyle','docked');
figure(101); clf(101); axis([-500 500 -500 500]);
figure(102); clf(102); axis([-500 500 -500 500]);

%% Generate simulated data
try
  load('measurements.mat');
catch
  generate_data();
  load('measurements.mat');
end

%% Plot simulated measurements
figure(101); box on; grid on; hold on;
title('Measurements');
for j = 1:numel(sensorMeasurements)
    plot(sensorMeasurements{1,j}.xMeas,sensorMeasurements{1,j}.yMeas,'.r');
end

%% structure for hypotheses and tracks
duration = 100;
structHyp = struct(...
    'wk',0,...                      % Probability for the hypothesis to exist
    'mk',zeros(4,1),...               % Mean of the hypothesis
    'Pk',zeros(4),...                    % Covariance of the hypothesis
    'Sk', zeros(4),...
    'Kk', 0,...
    'neta', 0);

HypP = structHyp;
HypP = structHyp;
nHyp = 0;

%% Filter
pD = .9;
pS = .99;

dT = 1;
F = [1 dT;...
     0 1];
F = [F zeros(size(F)); zeros(size(F)) F];

H = [1 0 0 0;...
     0 0 1 0];

R = 100*eye(2);
Q = 10*eye(4);

for k = 1:100
    % Prediction of existing targets
    for j = 1:numel(nHyp)
        if(nHyp==0) break; end
        HypP(j).wk = pS * HypP(j).wk;
        HypP(j).mk = F*HypP(j).mk;
        HypP(j).Pk = Q + F * HypP(j).Pk * F';
    end
  
    % Prediction of new births
    [x_pos,y_pos] = meshgrid([-250 250]);
    for j = 1:4
        HypP(nHyp+j).wk = .25 + abs(randn)/100;
        HypP(nHyp+j).mk = [x_pos(j);...
                           0.1;...
                           y_pos(j);...
                           0.1];
        HypP(nHyp+j).Pk = 25000*eye(4).*diag([1 .001 1 .001]);
%         figure(102); hold on;
%         h_ellips(j) = ellips(x_pos(j),y_pos(j),...
%                               diag([HypP(nHyp+j).Pk(1,1) HypP(nHyp+j).Pk(3,3)]),'r');
    end
    nHyp = nHyp + 4;
    
    % construction of PHD update components
    for j = 1:nHyp
        HypP(j).neta = H * HypP(j).mk;
        HypP(j).Sk = R + H * HypP(j).Pk * H';
        HypP(j).Kk = HypP(j).Pk * H' * pinv(HypP(j).Sk);
        HypP(j).Pk = (eye(4) - HypP(j).Kk * H) * HypP(j).Pk;
    end
    
    % update
    for j = 1:nHyp
        HypP(j).wk = (1 - pD) * HypP(j).wk;
    end
    L_val=0;
    for i_obs = 1:numel(sensorMeasurements{1,k}.xMeas)
        L_val=L_val+1;
        z = [sensorMeasurements{1,k}.xMeas(i_obs);...
             sensorMeasurements{1,k}.xMeas(i_obs)];
        w_sum = 0;
        for j = 1:nHyp
            HypP(L_val*nHyp+j).wk = pD*HypP(j).wk;
            HypP(L_val*nHyp+j).mk = HypP(j).mk + HypP(j).Kk*(z-HypP(j).neta);
            HypP(L_val*nHyp+j).Pk = HypP(j).Pk;
            w_sum = w_sum + HypP(L_val*nHyp+j).wk;
        end
        for j = 1:nHyp
          HypP(L_val*nHyp+j).wk = HypP(L_val*nHyp+j).wk/w_sum;
        end
    end
    nHyp = L_val*nHyp + nHyp;
    
    HypN = structHyp;
    % Merging/Pruning
    wk = extractfield(HypP,'wk');
    mk = extractfield(HypP,'mk');
    mk = reshape(mk,[numel(HypP(1).mk), numel(HypP)]);
    l = 0;
    I = find(wk >= 0.05);  % Pruning 
    I_full = I;
    while(~isempty(I))
      l=l+1;
      [~,j] = max(wk(I));   % index of maximum wt in pruned targets
      j = I(j);             % index of maximum wt in actual hypotheses 
      % Compute L(equality of gaussian components) with component j
      L_val = [];
      for i_merge = 1:numel(I)
        L_tmp = (HypP(I(i_merge)).mk - HypP(j).mk)' * pinv(HypP(I(i_merge)).Pk) * (HypP(I(i_merge)).mk - HypP(j).mk);
        L_val = [L_val L_tmp];
      end
      L = find(L_val<=1);
      I(L);
      HypN(l).wk = sum(wk(I(L)));
      HypN(l).mk = sum(repmat(wk(I(L)),4,1).*mk(:,I(L)),2)/HypN(l).wk;
      for i_sum = 1:numel(L)
        HypN(l).Pk = HypP(I(L(i_sum))).Pk + ...
          (HypN(l).mk - HypP(I(L(i_sum))).mk)*...
          (HypN(l).mk - HypP(I(L(i_sum))).mk)';
        HypN(l).Pk = wk(I(L(i_sum))) * HypN(l).Pk;
      end
      HypN(l).Pk = HypN(l).Pk/HypN(l).wk;
%       disp(HypN(l).Pk)
      I(L) = [];
    end
    HypP = HypN;
    nHyp = numel(HypN);
    for i = 1:numel(HypP)
      figure(102); hold on;
      plot(HypP(i).mk(1),HypP(i).mk(3),'.b');
    end
end