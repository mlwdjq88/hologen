%% test script
if ~exist('hl')||~ishandle(hl.hFigure)
    launch_hl;
end

%% parameter setting
hl.uieLambda.set(135);
hl.uieT.set(200);
hl.uieFocalLength.set(3);
hl.uieNA.set(0.0875);
hl.uieOffsetAngle.set(6);
hl.uieObscuration.set(0);
hl.uieMinFeature.set(1);
hl.uieFileName.set('test');
hl.uieOffset.set('[0, 0]');
hl.uieShift.set('[0, 0]');
hl.uipHologram.setSelectedIndex(uint8(1));
hl.cb(hl.uipHologram);

%% file path, select gds file here
hl.uieFilePath.set('C:\Users\mlwdjq\Documents\Code\hologen\Data\gds');


%% simulate pattern
hl.cb(hl.uibGenPattern);
drawnow;

%% simulate pattern
hl.cb(hl.uibGenGDS);
