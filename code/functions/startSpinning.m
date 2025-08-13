function  handle=startSpinning(string)

    f1 = figure; % creating figure
    set(f1, 'Toolbar', 'none', 'Menubar', 'none');
    f1pos = get(f1, 'Position'); % determine position
    framesize = [400 400]; % set framesize
    [handle] = StartLoadingSpinner('parent', f1, 'size', '32px', 'position', [(f1pos(3)/2)-framesize(1)/2 (f1pos(4)/2)-framesize(1)/2], 'framesize', framesize, 'text', string); % start spinner
      
    pause(0.01)

end