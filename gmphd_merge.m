function [ HypN ] = gmphd_merge( Hyp, prune_T, merge_U )
%Function to merge the hypothesis made my gmphd update
% INPUT VARIABLES
%   Hyp     - hypotheses
%   prune_T - pruning threshold
%   merge_U - merging threshold
%
% OUTPUT VARIABLES
%   HypN  - merged hypotheses 
%
% TODO: add Jmax for maximum number of gaussian terms

% Copy the structure format
HypN = Hyp;
HypN(2:end) = []; 

% Extract these values to avoid loop usage
wk = extractfield(Hyp,'wk');
mk = extractfield(Hyp,'mk');
mk = reshape(mk,[numel(Hyp(1).mk), numel(Hyp)]);

l_count = 0;
I = find(wk >= prune_T);  % Pruning 

while(~isempty(I))
    l_count = l_count+1;
    [~,j] = max(wk(I));   % index of maximum wt in pruned targets
    j = I(j);             % index of maximum wt in actual hypotheses 
    
    % Compute L(equality of gaussian components) with component j
    L_val = [];
    for i_merge = 1:numel(I)
        L_tmp = (Hyp(I(i_merge)).mk - Hyp(j).mk)' * pinv(Hyp(I(i_merge)).Pk) * (Hyp(I(i_merge)).mk - Hyp(j).mk);
        L_val = [L_val L_tmp];
    end
    
    % Find all L that are less than the merging threshold
    L = find(L_val <= merge_U);
    
    % Merge the hypotheses that were found close in L
    HypN(l_count).wk = sum(wk(I(L)));
    HypN(l_count).mk = sum(repmat(wk(I(L)),4,1).*mk(:,I(L)),2)/HypN(l_count).wk;
    for i_sum = 1:numel(L)
        HypN(l_count).Pk = Hyp(I(L(i_sum))).Pk + ...
          (HypN(l_count).mk - Hyp(I(L(i_sum))).mk)*...
          (HypN(l_count).mk - Hyp(I(L(i_sum))).mk)';
        HypN(l_count).Pk = wk(I(L(i_sum))) * HypN(l_count).Pk;
    end
    HypN(l_count).Pk = HypN(l_count).Pk/HypN(l_count).wk;
    
    % Remove hypotheses that are already merged and added as new hypotheses
    I(L) = [];
end

end