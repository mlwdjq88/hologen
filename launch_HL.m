addpath('../mpm');

% mic library
mpm addpath

% Add this root directory to path
[cDirThis, cName, cExt] = fileparts(mfilename('fullpath'));
addpath(genpath(cDirThis));



% Build HL UI
hFigure = figure(...
    'name', 'Hololens GDS generator GUI v1.200924',...
    'Units', 'pixels',...
    'Position', [500, 300,  900, 600],...
    'handlevisibility','off',... %out of reach gcf
    'numberTitle','off',...
    'Toolbar','none',...
    'Menubar','none');

hl = hologen.ui.HL_generator;
hl.build(hFigure,-500,-300);


%% generating holo lens for 1d LSI (visible light)
% hologen.simulation.genHL_LSI_onAxis_Visible();

%% generating holo lens for 2d LSI (visible light)
% hologen.simulation.genHL_QWLSI_onAxis_Visible();

%% generating holo lens for 2d LSI (EUV light)
%% hologen.simulation.genHL_QWLSI_onAxis_EUV();
