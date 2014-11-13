function [er,output2,output3] = errorFunction( errFcn, varargin )
% err = errorFunction( errFcn, arguments )
%   just returns the evaluation of errFcn( arguments ).
%   This seems pointless, but the point is that this function
%   has memory of all the errors, so you look at this later.
%
% errHist = errorFunction()
%   will return the history of errors and also clear the memory
%
% ... = errorFunction( ticReference )
%   will allow the function to enable time logging
%   It will subtract off the time used for the error function
%
% [errHist, timeLog] = errorFunction( )
%   will output the time logging information.
%   The time log shows at what time, relative to ticReference,
%   the error was recorded
%
% The function uses persistent variables so don't forget to clear it
%   and set it with a new ticReference!
%
% Stephen Becker, March 6 2014 stephen.beckr@gmail.com
% March 11 2014, updated to subtract off the time involved
%   for computing the error. Also allows error function
%   to return a row vector instead of a scalar
% March 14 2014, if ticReference is really a cell array with
%   {ticReference,maxTime}, then this will throw an error
%   after reaching "maxTime" seconds

persistent erHist timeLog ticReference timeCorrection maxTime
if nargin <= 1
    if nargin == 1 && ~isempty(errFcn)
        if iscell( errFcn )
            ticReference = errFcn{1};
            maxTime      = errFcn{2};
        else
            maxTime     = Inf;
            ticReference = errFcn;
        end
    else
        ticReference = [];
    end
    er = erHist;
    if nargout >= 2, output2 = timeLog; end
    if nargout >= 3, output3 = timeCorrection; end
    erHist = [];
    timeLog= [];
    timeCorrection = [];
    return
end
if isempty( maxTime ), maxTime = Inf; end

tc1 = tic;
er = errFcn( varargin{:} );
erTime = toc(tc1);
erHist(end+1,:) = er;

% And mark at what time we took this error measurement
if ~isempty( ticReference )
    % And we also want to add "erTime" to ticReference
    % but this isn't straightforward since ticReference is some huge
    % uint64 and not in seconds. So instead we update a correction
    if isempty( timeCorrection), timeCorrection = 0; end
    timeCorrection = timeCorrection + erTime;
    timeLog(end+1) = toc( ticReference ) - timeCorrection;
    
    if timeLog(end) > maxTime
        error('errorFunction:timedOut','Reached maximum allowed time');
    end
end