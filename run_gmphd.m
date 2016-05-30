%% Clear and csstomize MATLAB 
clc;
clear;
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultFigureWindowStyle','docked');
figure(101); clf(101); axis([-500 500 -500 500]);
figure(102); clf(102); axis([-500 500 -500 500]);

%% GMPHD modelling params
model.prune_T = .05;
model.merge_U = 10;
model.pD = 1;
model.pS = .9;
model.dT = 1;
model.noise_process = .1;
model.noise_sensor = 300;
model.F = [1 model.dT;...
           0 1];
model.F = [model.F zeros(size(model.F));...
           zeros(size(model.F)) model.F];
model.Q = model.noise_process * eye(4);
model.H = [1 0 0 0;...
           0 0 1 0];
model.R = model.noise_sensor * eye(2);

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

%% plot groundtruth
figure(102); box on; grid on; hold on;
title('Tracks');
for j = 1:numel(groundTruth)
    h_102(1) = plot(groundTruth(j).track.x, groundTruth(j).track.y,'-b');
end
h_102(2) = plot(-500,500,'.k');
legend(h_102,'Groundtruth','GMPHD');

%% structure for hypotheses and tracks
duration = 100;
structHyp = struct(...
    'wk',-1,...                % Probability for the hypothesis to exist, keep -1 to lets functions know first iteration
    'mk',zeros(4,1),...        % Mean of the hypothesis
    'Pk',zeros(4),...          % Covariance of the hypothesis
    'Sk', zeros(4),...
    'Kk', 0,...
    'neta', 0);

HypP = structHyp;

% values for test
% HypP.wk = 1;
% HypP.mk = [-100, 0, -400, 7]';
% Hyp.Pk = eye(4);

%% Filter
for k = 1:numel(sensorMeasurements)
    % Prediction
    HypP = gmphd_predict(HypP, model, sensorMeasurements{1,k},k);
    
    % Update
    HypP = gmphd_update( HypP, model, sensorMeasurements{1,k});
    
    % Prune and Merge
    HypP = gmphd_merge( HypP, model.prune_T, model.merge_U );
    
    wk = extractfield(HypP,'wk');
    disp(['sum of wk:' num2str(sum(wk))])
    figure(102);hold on; box on; grid on;
    for i = 1:round(sum(wk))
        if(i>numel(wk))
            break;
        end
        plot(HypP(i).mk(1),HypP(i).mk(3),'.k');
    end
    % State extraction
%     Xk = [];
%     for i = 1:numel(HypP)
%         if(HypP(i).wk > .5)
%             for j = 1:round(HypP(i).wk)
%                 Xk = [Xk, HypP(i).mk];
%             end
%         end
%     end
%     
%     if(~isempty(Xk))
%         figure(102); hold on;
%         plot(Xk(1,:),Xk(3,:),'.k');
%     end
    
    
    pause(.01);
    disp(['Iteration:' num2str(k)]);
end