function [] = StopLoadingSpinner_f(varargin)
%STOPLOADINGSPINNER_F Summary of this function goes here
%   Detailed explanation goes here
%mandatory inputs
%   'hdls'               -> Java object handles to the spinner and container
%possible inputs (case insensitive):
%   'successText'      -> String to display
%   'stopfreeze'       -> boolean whether icon shall freeze(true) or vanish (false) after being stopped

% defaul values
In.stopfreeze = false;
% parse property name - value pairs into input struct
i = 1;
while i <= length(varargin)
    if strcmpi(varargin{i},'hdls')
        In.hdls = varargin{i+1};
    elseif strcmp(varargin{i},'successtext')
        In.successtext = varargin{i+1};
    elseif strcmp(varargin{i},'stopfreeze')
        In.stopfreeze = varargin{i+1};
    else
        error(['Unknown property name: ', varargin{i}]);
    end
    i = i+2;
end
% stop spinner and...
    if isfield(In, 'successtext')
    %... show text
        In.hdls.jObj.setPaintsWhenStopped(true);
        In.hdls.jObj.setBusyText(In.successtext);
    end
    if In.stopfreeze == true
    %... freeze it
        In.hdls.jObj.setPaintsWhenStopped(true);
    end
    In.hdls.jObj.stop;
end
