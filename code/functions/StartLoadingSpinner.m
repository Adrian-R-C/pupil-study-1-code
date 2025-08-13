function [hdls, In] = StartLoadingSpinner_f(varargin)
%STARTLOADINGSPINNER_F Summary of this function goes here
%   Detailed explanation goes here
%possible inputs (case insensitive):
%   'position'          -> [top bottom]
%   'backgroundColor'   -> [R G B]
%   'size'              -> either '16px' or '32px'
%   'framesize'         -> [width height] in px. Size of the surrounding frame. Default: Icon size
%   'parent'            -> handle to either a Uipanel, a figure or a Uitoolbar
%   'text'              -> text to be shown under the spinner as String
%   'dotcolor'          -> 'darkblue' or 'white'

% set default values
In.text = '';
In.position = [0 0];
In.size = '16px';
In.dotcolor = 'darkblue';
In.framesize = [18 18];
% parse property name - value pairs into input struct
i = 1;
while i <= length(varargin)
    if strcmpi(varargin{i},'position')
        In.position = varargin{i+1};
    elseif strcmp(varargin{i},'backgroundcolor')
        In.bgcolor = varargin{i+1};
    elseif strcmpi(varargin{i},'size')
        In.size = varargin{i+1};
    elseif strcmpi(varargin{i},'parent')
        In.parent = varargin{i+1};
    elseif strcmpi(varargin{i},'text')
        In.text = varargin{i+1};
    elseif strcmpi(varargin{i},'dotcolor')
        In.dotcolor = varargin{i+1};
    elseif strcmpi(varargin{i},'framesize')
        In.framesize = varargin{i+1};
    else
        error(['Unknown property name: ', varargin{i}]);
    end
    i = i+2;
end
% get possible sizes
iconsSizeEnums = javaMethod('values','com.mathworks.widgets.BusyAffordance$AffordanceSize');
% set size
switch In.size
    case '16px'
        In.enumsize = iconsSizeEnums(1);
        In.pxsize = 16;
    case '32px'
        In.enumsize = iconsSizeEnums(3);
        In.pxsize = 32;
        if In.framesize(1) == 18 && In.framesize(2) == 18
            In.framesize = [34 34];
        end
    otherwise
        error('Wrong size option!')
end
% set text
if ~isfield(In, 'text')
    In.text = '';
end
% create Java Object
jObj = com.mathworks.widgets.BusyAffordance(In.enumsize, In.text);  % icon, label
% set dotcolor (default: false)
switch In.dotcolor
    case 'darkblue'
        jObj.useWhiteDots(false);  
    case 'white'
        jObj.useWhiteDots(true);  
    otherwise
        error('Wrong dotcolor option!');
end
% set position and parent
if ~(strcmp('figure', get(In.parent, 'type')) || strcmp('uipanel', get(In.parent, 'type')) || strcmp('uitoolbar', get(In.parent, 'type')))
    error('Wrong type of parent handle');
else
    [hdls.hdlJObj, hdls.hdlContainer] = javacomponent(jObj.getComponent, [In.position(1),In.position(2),In.framesize(1),In.framesize(2)], In.parent);
end
% set background color
if isfield(In, 'bgcolor')
    jObj.getComponent.setBackground(java.awt.Color(In.bgcolor(1), In.bgcolor(2), In.bgcolor(3)));
end
% finally show spinner
jObj.start;
hdls.jObj = jObj;
end
