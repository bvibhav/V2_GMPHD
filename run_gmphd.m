%% Intialize 
clc;
clear;

set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultFigureWindowStyle','docked');

figure(101); clf(101); axis([-500 500 -500 500]);
figure(102); clf(102); axis([-500 500 -500 500]);

%% Generate simulated data
% generate_data();
load('measurements.mat');

%% Plot simulated measurements
figure(101); box on; grid on; hold on;
title('Measurements');
for i = 1:numel(sensorMeasurements)
    plot(sensorMeasurements{1,i}.xMeas,sensorMeasurements{1,i}.yMeas,'.r');
end
    
%%