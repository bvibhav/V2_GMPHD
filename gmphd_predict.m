function Hyp = gmphd_predict( Hyp, model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nHyp = numel(Hyp);
if Hyp(1).wk==-1
    nHyp = 0;
end

% Prediction of existing targets
for j = 1:nHyp
    Hyp(j).wk = model.pS * Hyp(j).wk;
    Hyp(j).mk = model.F * Hyp(j).mk;
    Hyp(j).Pk = model.Q + model.F * Hyp(j).Pk * model.F';
end

% Prediction of new births
[x_pos,y_pos] = meshgrid([-250 250]);
for j = 1:4
    Hyp(nHyp+j).wk = 1 + abs(randn)/1000;
    Hyp(nHyp+j).mk = [x_pos(j);...
                       0.1;...
                       y_pos(j);...
                       0.1];
    Hyp(nHyp+j).Pk = 2000000*eye(4).*diag([1 .001 1 .001]);
end

end

