function generate_data()
%GENERATE_DATA Summary of this function goes here
%   Generate simulated data for GMPHD filter

figure(101); clf(101); hold on; axis([-500 500 -500 500]); box on; grid on;
figure(102); clf(102); hold on; axis([-500 500 -500 500]); box on; grid on;

% Initialize targets
nTargets = 3;
simTime = 100;

for i = 1:nTargets
    target(i).x = 0;
    target(i).y = 0;
    target(i).vx = 1;
    target(i).vy = 1;
    target(i).spawn = 'n';
    target(i).spawnTime = -1;
    target(i).spawnTarget.vx = 0;
    target(i).spawnTarget.vy = 0;
end

% Initialize individual target properties
target(1).x = -100;
target(1).y = -400;
target(1).vx = 0;
target(1).vy = 7;

target(2).x = 100;
target(2).y = 400;
target(2).vx = 0;
target(2).vy = -7;
% target(2).spawn = 'y';
% target(2).spawnTime = 30;
% target(2).spawnTarget.vx = 2;
% target(2).spawnTarget.vy = -4;
target(2).x = 300;
target(2).y = 0;
target(2).vx = -7;
target(2).vy = 0;

target(3).x = -300;
target(3).y = -300;
target(3).vx = 6;
target(3).vy = 6;
target(3).spawn = 'y';
target(3).spawnTime = simTime/2;
target(3).spawnTarget.vx = -2;
target(3).spawnTarget.vy = 6;

% Initialize target states
for i = 1:nTargets
    target(i).state = [ target(i).x;...
                        target(i).vx;...
                        target(i).y;...
                        target(i).vy];
    groundTruth(i).track.x = [];
    groundTruth(i).track.y = [];
    groundTruth(i).track.t = [];
end

% Sampling time
dT = 1;

% State transition model
F = [1 dT;...
     0 1];
   
% Input control matrix or matrix for accelerative white gaussian noise
G = [dT^2; ...
    dT];
  
% For 2D systems
F = [F zeros(size(F)); zeros(size(F)) F];
G = [G zeros(size(G)); zeros(size(G)) G];

% noise factor for generating trajectories
pNoise = .05;

sensorMeasurements = {};
for i = 1:simTime
  sensor.xMeas = [];
  sensor.yMeas = [];
  
  for i_targetCounter = 1:nTargets
    figure(101); hold on;
    plot(target(i_targetCounter).state(1), target(i_targetCounter).state(3),'.b');
    groundTruth(i_targetCounter).track.x = [groundTruth(i_targetCounter).track.x    target(i_targetCounter).state(1)];
    groundTruth(i_targetCounter).track.y = [groundTruth(i_targetCounter).track.y    target(i_targetCounter).state(3)];
    groundTruth(i_targetCounter).track.t = [groundTruth(i_targetCounter).track.t    i];
    
    if(target(i_targetCounter).spawn=='y' & target(i_targetCounter).spawnTime==i)
      nTargets = nTargets + 1;
      target(nTargets).x = target(i_targetCounter).state(1);
      target(nTargets).y = target(i_targetCounter).state(3);
      target(nTargets).vx = target(i_targetCounter).spawnTarget.vx;
      target(nTargets).vy = target(i_targetCounter).spawnTarget.vy;
      target(nTargets).spawn = 'n';
      target(nTargets).spawnTime = -1;
      target(nTargets).spawnTarget.vx = 0;
      target(nTargets).spawnTarget.vy = 0;
      target(nTargets).state = [ target(nTargets).x;...
                          target(nTargets).vx;...
                          target(nTargets).y;...
                          target(nTargets).vy];
        groundTruth(nTargets).track.x = [];
        groundTruth(nTargets).track.y = [];
        groundTruth(nTargets).track.t = [];
    end
    
    SQRT = sqrt(target(i_targetCounter).state(2)^2 + target(i_targetCounter).state(4)^2);
    pNoiseX = pNoise * target(i_targetCounter).state(2)/SQRT + .05;
    pNoiseY = pNoise * target(i_targetCounter).state(4)/SQRT + .05;

    wk = [pNoiseX;pNoiseY] .* wgn(2,1,10);
    target(i_targetCounter).state = F * target(i_targetCounter).state+ G * wk;
    
    if rand < .95 % Pdetection
      sensor.xMeas = [sensor.xMeas target(i_targetCounter).state(1)+5*randn];
      sensor.yMeas = [sensor.yMeas target(i_targetCounter).state(3)+5*randn];
    end
  end
  
  PRN = poissrnd(12.5);
  disp(['PRN is ' num2str(PRN)])
  for i_randCounter = 1:PRN
    sensor.xMeas = [sensor.xMeas -500+1000*rand];
    sensor.yMeas = [sensor.yMeas -500+1000*rand];
  end
  % Add random noise points
%   for i_randCounter = 1:floor(rand*5)
%     sensor.xMeas = [sensor.xMeas -500+1000*rand];
%     sensor.yMeas = [sensor.yMeas -500+1000*rand];
%   end
  pause(.01);
  sensorMeasurements{i} = sensor;
  figure(102); hold on;
  plot(sensor.xMeas,sensor.yMeas,'.r');
end

save('measurements.mat', 'sensorMeasurements', 'groundTruth');

% disp('Press any key to continue'); pause;
% figure(101); clf(101); box on; grid on; hold on;
% figure(102); clf(102); box on; grid on; hold on;

end

