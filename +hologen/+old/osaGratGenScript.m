% Script to generate zts gratings and OSA
addpath('./mgds');
cXLSTemplate = 'cellRef.xlsx';
mgds_SUPER = MGDS();
mgds_SUPER.init('ZTS_grating_OSA_package_mask', ...
                        'autogen structures', true);

mgds = MGDS();
mgds.clearAll();

% Render to super MGDS:
bProduction = true;
bGratingLabelsOff = false; % speeds up render

% Fields:
% MAIN_APERTURE, ALIGN_15_30, LETTER_L,
% LETTER_R, ALIGN_15_35, MAIN_GRATING

%% Open:

cField = 'OPEN';

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
mgds.makeRect(cField, 80, 80, 'center', 0);

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end

%% Open:

cField = 'SMALL_OPEN';

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
mgds.makeRect(cField, 30, 30, 'center', 0);

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end

%% Open:

cField = 'BIGGER_OPEN';

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
mgds.makeRect(cField, 125, 125, 'center', 0);

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end
%% Open:

cField = 'BIG_OPEN';

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
mgds.makeRect(cField, 200, 200, 'center', 0);

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end

%% Open:

cField = 'BIG_OPEN_CIRCLE';

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
mgds.makeCircle(cField, 87.5, [0,0], 100, 0);


if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end

%% Grats:
cField = 'GRAT_HV';
dPitch = 0.6;
dCD = dPitch/2;

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);


% make prim:
gPrimitiveName = sprintf('%s_pr', cField);
gPrimitive = mgds.makeRect(gPrimitiveName, dCD, dCD, [0, 0], 0);


gAtomName = sprintf('%s_at', cField);
gAtom = mgds.makeRef(gAtomName, gPrimitive, [0, 0; dCD, dCD], 0);

% Make grating:
mgds.makeARef(cField, gAtom, 133, 133, dPitch, dPitch, -40, -40);
if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end




cField = 'GRAT_45';
mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
gPrimitiveName = sprintf('%s_pr', cField);
gPrimitive = mgds.makeRect(gPrimitiveName, dCD, dCD, [0, 0], 0);
gAtomName = sprintf('%s_at', cField);
gAtom = mgds.makeRef(gAtomName, gPrimitive, [0, 0], 45*pi/180);
mgds.makeARef(cField, gAtom, 188, 188, dPitch/sqrt(2), dPitch/sqrt(2), -40, -40);


if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end


%% Letters
cField = 'LETTER_L';
mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
sfLetterL = 'sfLetterL';
mgds.makePolygonText(sfLetterL, 0,0, 0, 'L', 60, 'center', false, 0);
mgds.makeRef(cField, sfLetterL, [-2.5, 6], 0);


if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end

%% Letters
cField = 'LETTER_R';
mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
sfLetterR = 'sfLetterR';
mgds.makePolygonText(sfLetterR, 0,0, 0, 'R', 60, 'center', false, 0);
mgds.makeRef(cField, sfLetterR, [0, 6], 0);

cFieldBF = 'LETTER_R_BF';
% Let's make R brightField:
mgds.scheduleBinaryOperation(cFieldBF, cField, 'BIGGER_OPEN', 'XOR');

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end

%% Windmill:

cField = 'WINDMILL';
sfField = 'WINDMILL_SUB';
mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);

gPrimitive = sprintf('%s_prim', cField);
mgds.makeShape(gPrimitive, [0, 20, 40-.25, 0, 0], [0, 0, 40-.25, 20, 0], 0);
mgds.makeRef(sfField, gPrimitive, [0, 0], 0);
mgds.makeRef(sfField, gPrimitive, [80, 0], pi/2);
mgds.makeRef(sfField, gPrimitive, [80, 80], pi);
mgds.makeRef(sfField, gPrimitive, [0, 80], 3*pi/2);

gRect = sprintf('%s_rect', cField);
mgds.makeRect(gRect, .5, .5, 'center', 0);
mgds.makeRef(sfField, gRect, [40, 40], 0);


gLines = sprintf('%s_lines', cField);
% make lines:
% for k = 1:2:38
%     mgds.makeRect(gLines, k, 1, [40 - k/2, 38-k], 0);
% end

for k = 1:39
    mgds.makeRect(gLines, k, .5, [40 - k/2, 39.5-k], 0);
end

mgds.makeRef(sfField, gLines, [0, 0], 0);
mgds.makeRef(sfField, gLines, [80, 0], pi/2);
mgds.makeRef(sfField, gLines, [80, 80], pi);
mgds.makeRef(sfField, gLines, [0, 80], 3*pi/2);

mgds.makeRef(cField, sfField, [-40, -40], 0);

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end
%% Main grating:

cField = 'MAIN_GRATING';
dPitch = 1.644*sqrt(2);
dCD = dPitch/2;

mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);
gPrimitiveName = sprintf('%s_pr', cField);
gPrimitive = mgds.makeRect(gPrimitiveName, dCD, dCD, 'center', 0);

idx = 0:(dPitch*sqrt(2))/2:43*dPitch*sqrt(2);
idx = idx - mean(idx);
[dX, dY] = meshgrid(idx);
nIdx = sqrt(dX(:).^2 + dY(:).^2) < 20/2;

dX(nIdx) = [];
dY(nIdx) = [];

mgds.makeRef(cField, gPrimitive, [dX(:), dY(:)], 45*pi/180);

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end


%% Main Obs:

cField = 'MAIN_APERTURE';
mgds.init(cField, 'autogen structures', true, 'hide grating labels', bGratingLabelsOff);

R = 6;
th = linspace(0, 2*pi, 51);
dX = R*cos(th);
dY = R*sin(th);

mgds.makeShape(cField, dX, dY, 0);
if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end
%% Align gratings

cField = 'ALIGN_TOP';
mgds.init(cField, 'autogen structures', true);

T = 24;
DC = 12;
N = 2;
cVLinesName = sprintf('%s_vl', cField);
for k = -N:N
    if k == 0
        height = DC*3;
    else
        height = DC*2.5;
    end
    mgds.makeRect(cVLinesName, DC,height, [k*T - DC/2, -DC*2.5], 0);
end

mgds.makeRef(cField, cVLinesName, [-1.5*T, -DC*2], 0);
mgds.makeRef(cField, cVLinesName, [DC*2, 1.5*T], pi/2);

cField = 'ALIGN_BOTTOM';

T2 = 28;
DC = 12;
N = 2;
cVLinesName2 = sprintf('%s_vl2', cField);
for k = -N:N
    mgds.makeRect(cVLinesName2, DC, DC*2.5, [k*T2 - DC/2, -DC*2.5], 0);
end

mgds.makeRef(cField, cVLinesName2, [-1.5*T, -DC*5], 0);
mgds.makeRef(cField, cVLinesName2, [DC*5, 1.5*T], pi/2);


mgds.scheduleBinaryOperation('ALIGN_BOTTOM_BF', 'ALIGN_BOTTOM', 'BIG_OPEN', 'XOR');
mgds.scheduleBinaryOperation('ALIGN_TOP_BF', 'ALIGN_TOP', 'BIG_OPEN', 'XOR');

if bProduction
    mgds_SUPER.import(mgds);
else
    mgds.makeGDS();
end



%% Create refereces to subfields
mgds_SUPER.makeRefFromXLS(cXLSTemplate, 'OSA');
mgds_SUPER.makeRefFromXLS(cXLSTemplate, 'GRATING');
% mgds_SUPER.makeRefFromXLS(cXLSTemplate, 'OSA_XORMASK');
% mgds_SUPER.makeRefFromXLS(cXLSTemplate, 'GRATING_XORMASK');

%% Schedule binary operations:
% mgds_SUPER.scheduleBinaryOperation('OSA_MAIN', 'OSA_XORMASK', 'OSA', 'XOR');
% mgds_SUPER.scheduleBinaryOperation('GRATING_MAIN', 'GRATING_XORMASK', 'GRATING', 'XOR');




mgds_SUPER.makeGDS();


    

