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
HypU = structHyp;
nHyp = 0;

%% Filter
pDetection = .9;
H = [1 0 0 0;...
     0 0 1 0];
R = 100*eye(2);

for k = 1:1
    % Prediction of new births
    [x_pos,y_pos] = meshgrid(-325:325:325);
    for j = 1:9
        nHyp = nHyp + 1;
        HypP(j).wk = .1;
        HypP(j).mk = [x_pos(j);...
                      0;...
                      y_pos(j);...
                      0];
        HypP(j).Pk = 10000*eye(4).*diag([1 .001 1 .001]);

%         figure(102); hold on;
%         h_ellips(j) = ellips(x_pos(j),y_pos(j),Hyp(j).cov,'r');
    end
    
    % Prediction of existing targets
    if nHyp > 0;
        NONE = -1;
    end
    
    % construction of PHD update components
    for j = 1:nHyp
        HypP(j).neta = H * HypP(j).mk;
        HypP(j).Sk = R + H * HypP(j).Pk * H';
        HypP(j).Kk = HypP(j).Pk * H' * pinv(HypP(j).Sk);
        HypP(j).Pk = (eye(4) - HypP(j).Kk * H) * HypP(j).Pk;
    end
    
    % update
    for j = 1:nHyp
        HypU(j).wk = (1 - pDetection) * HypP(j).wk;
        HypU(j).mk = HypP(j).mk;
        HypU(j).Pk = HypP(j).Pk;
    end
    L=0;
    for i_obs = 1:numel(sensorMeasurements{1,k}.xMeas)
        L=L+1;
        z = [sensorMeasurements{1,k}.xMeas(i_obs);...
             sensorMeasurements{1,k}.xMeas(i_obs)];
        w_sum = 0;
        for j = 1:nHyp
            HypU(L*nHyp+j).wk = pDetection*HypU(j).wk;
            HypU(L*nHyp+j).mk = HypU(j).mk + HypP(j).Kk*(z-HypP(j).neta);
            HypU(L*nHyp+j).Pk = HypP(j).Pk;
            w_sum = w_sum + HypU(L*nHyp+j).wk;
        end
        for j = 1:nHyp
          HypU(L*nHyp+j).wk = HypU(L*nHyp+j).wk/w_sum;
        end
    end
    nHyp = L*nHyp + nHyp;
end