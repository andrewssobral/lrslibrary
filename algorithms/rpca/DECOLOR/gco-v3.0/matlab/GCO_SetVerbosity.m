function GCO_SetVerbosity(Handle,Level)
% GCO_SetVerbosity   Print status messages during Expansion/Swap.
%   Level 0 prints no output (full speed).
%   Level 1 prints cycle-level status messages.
%   Level 2 prints expansion/swap-level status messages.
%   The current energy is printed as
%        E=Total (E=DataCost+SmoothCost+LabelCost)
%   At level 2, the size of each binary graph cut problem is also
%   printed (# vars).
%
%   Note that printing may have an effect on overall run time (tic/toc)
%   though the internal computation times are printed in milliseconds
%   and exclude the time to print.
%
%   Example:
%     >> GCO_SetVerbosity(Handle,2); % Level 2 output
%     >> GCO_Expansion(Handle);
%     gco>> starting alpha-expansion w/ adaptive cycles
%     gco>> initial energy: 	E=18 (E=17+0+1)
%     gco>>   after expansion(3): 	E=17 (E=14+1+2);	 4 vars; 	(1 of 8);	 0.003 ms
%     gco>>   after expansion(7): 	E=15 (E=11+1+3);	 2 vars; 	(2 of 8);	 0.002 ms
%             ...
%     gco>> after cycle  1: 	E=12 (E=6+2+4); 	8 expansions(s);

GCO_LoadLib();
gco_matlab('gco_setverbosity',Handle,int32(Level));
end
