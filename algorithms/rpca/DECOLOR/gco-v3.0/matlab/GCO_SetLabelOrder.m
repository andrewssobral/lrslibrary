function GCO_SetLabelOrder(Handle,Order)
% GCO_SetLabelOrder   Set label order for Expansion/Swap moves.
%   GCO_SetLabelOrder(Handle,Order) tells Expansion/Swap to select labels
%   in a specific order when Order contains integers from 1..NumLabels.
%   By default, Expansion/Swap use a consistent, prescribed order of labels
%   in until convergence. 
%   Example:
%      GCO_SetLabelOrder(Handle,5:10); % only operate on labels 5..10
%      GCO_SetLabelOrder(Handle,randperm(NumLabels)); % random label order
%      

GCO_LoadLib();
gco_matlab('gco_setlabelorder',Handle,int32(Order));
end
