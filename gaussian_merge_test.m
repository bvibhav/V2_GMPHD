clc;
clear;

set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultFigureWindowStyle','docked');

figure(101); clf(101); axis([-500 500 -500 500]);

m1 = [0 0]';
P1 = 1000*eye(2).*diag([1 1]);
figure(101); hold on;
h_ellips(1) = ellips(m1(1),m1(2),diag([P1(1,1) P1(2,2)]),'r')

m2 = [100 0]';
P2 = 1000*eye(2).*diag([1 1]);
figure(101); hold on;
h_ellips(1) = ellips(m2(1),m2(2),diag([P2(1,1) P2(2,2)]),'b')

L = (m1-m2)' * (pinv(P1).*pinv(P2)) * (m1-m2)