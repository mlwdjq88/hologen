classdef HL_generator < mic.Base
    
    
    properties (Constant)
        dWidth  = 870;
        dHeight =  600;
        
        % Axes tab IDS
        U8SIM          = 1
    end
    
    properties
        
        cAppPath = fileparts(mfilename('fullpath'))
        
        % Graphical elements
        hFigure     % Main figure (not overwritable)
        
        % valuables
        stShapes
        dTh
        dI
        dHL
        x_um
        y_um
        
        % axis tab
        uitgAxesDisplay     % displays axes:simulation
        haFieldAmp
        
        
        % parameters
        hpPara
        uieLambda
        uieT
        uieFocalLength
        uieNA
        uieOffsetAngle
        uieObscuration
        uieMinFeature
        uieFileName
        uieOffset
        uieShift
        uipHologram
        
        % control
        hpControl
        uieFilePath
        uibSetPath
        uilFileList
        uibGenGDS
        uibOpenGDS
        uibGenPattern
        
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        function this = HL_generator()
            this.init()
        end
        
        
        
        function init(this)
            
            % axis tab
            this.uitgAxesDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Simulation'});
            
            % parameters
            this.uieLambda       = mic.ui.common.Edit('cLabel', 'Wavelength(nm)', 'cType', 'd');
            this.uieT     = mic.ui.common.Edit('cLabel', 'T(um)', 'cType', 'd');
            this.uieFocalLength = mic.ui.common.Edit('cLabel', 'Focal length(mm)', 'cType', 'd');
            this.uieNA  = mic.ui.common.Edit('cLabel', 'NA', 'cType', 'd');
            this.uieOffsetAngle  = mic.ui.common.Edit('cLabel', 'Offset angle(deg)', 'cType', 'd');
            this.uieObscuration      = mic.ui.common.Edit('cLabel', 'Obs. radius(mm)', 'cType', 'd');
            this.uieMinFeature      = mic.ui.common.Edit('cLabel', 'Min feature(nm)', 'cType', 'd');
            this.uieFileName      = mic.ui.common.Edit('cLabel', 'File name', 'cType', 'c');
            this.uieOffset         = mic.ui.common.Edit('cLabel', 'Center offset(mm)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieShift         = mic.ui.common.Edit('cLabel', 'Center shift(mm)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipHologram    = mic.ui.common.Popup('cLabel', 'Hologram', 'ceOptions',...
                {'QWLSI','LSI'}, 'lShowLabel', true,'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieLambda.set(13.5);
            this.uieT.set(540);
            this.uieFocalLength.set(3);
            this.uieNA.set(0.0875);
            this.uieOffsetAngle.set(6);
            this.uieObscuration.set(0);
            this.uieMinFeature.set(1);
            this.uieFileName.set('New');
            this.uieOffset.set('[0, 0]');
            this.uieShift.set('[0, 0]');
            this.uipHologram.setSelectedIndex(uint8(1));
            this.dTh = 32;
            
            % control
            this.uieFilePath    = mic.ui.common.Edit('cLabel', 'GDS file path', 'cType', 'c','fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSetPath    = mic.ui.common.Button('cText', 'Set path', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uilFileList         = mic.ui.common.List('cLabel', 'GDS file list', ...
                'lShowDelete', false, 'lShowMove', false, 'lShowRefresh', false,'fhOnChange', @(src, evt)this.cb(src));
            this.uibGenGDS   = mic.ui.common.Button('cText', 'Gen GDS', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibGenPattern   = mic.ui.common.Button('cText', 'Gen pattern', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibOpenGDS   = mic.ui.common.Button('cText', 'Open GDS', 'fhDirectCallback', @(src, evt)this.cb(src));
            cDataDir = fullfile(this.cAppPath,'..','..','Data','gds');
            this.uieFilePath.set(cDataDir);
            
            
        end
        
        % Callback handler
        function cb(this, src,evt)
            switch src
                case {this.uieOffset, this.uieShift}
                    this.validateCouplesEditBox(src, '[0, 0]');
                    
                case this.uieFilePath
                    cDataDir = this.uieFilePath.get();
                    % get data file list
                    dirOutput=dir(fullfile(cDataDir,'*.gds'));
                    fileNames={dirOutput.name}';
                    this.uilFileList.setOptions(fileNames);
                    
                case this.uibSetPath
                    cDataDir = fullfile(this.cAppPath,'..','..','Data','gds');
                    dataPath = uigetdir(cDataDir);
                    this.uieFilePath.set(dataPath);
                    
                case this.uilFileList
                    datalist = this.uilFileList.getOptions();
                    value = this.uilFileList.getSelectedIndexes();
                    this.uieFileName.set(datalist{value}(1:end-4));
                    
                case this.uipHologram
                    switch this.uipHologram.getSelectedIndex()
                        case 1
                            this.dTh = 32; % QWLSI
                        case 2
                            this.dTh = 8; % LSI
                    end
                    
                case this.uibGenGDS
                    this.genGDS();
                    this.cb(this.uieFilePath);
                case this.uibGenPattern
                    this.genPattern();
                    this.replot(this.U8SIM);
                    
                case this.uibOpenGDS
                    cDataDir = this.uieFilePath.get();
                    datalist = this.uilFileList.getOptions();
                    value = this.uilFileList.getSelectedIndexes();
                    filename= fullfile(cDataDir,datalist{value});
                    winopen(filename);
            end
        end
        
        % simulate hologram
        function genPattern(this)
            T_um = this.uieT.get(); % grating pitch
            f_um = this.uieFocalLength.get()*1e3; % lens focal length
            lambda_um = this.uieLambda.get()*1e-3; % wavelength
            NA = this.uieNA.get();
            offset = eval(this.uieOffset.get())*1e3;
            xOffset = offset(1);
            yOffset = offset(2);
            N = 2000;
            subR_um=f_um*tan(asin(NA)); % lens radius
            this.x_um = linspace(-subR_um,subR_um,N);
            this.y_um = linspace(-subR_um,subR_um,N);
            [xp,yp] = meshgrid(this.x_um,this.y_um);
            pupil = zeros(N);
            pupil(xp.^2+yp.^2<= subR_um.^2) = 1;
            xp = xp + xOffset;
            yp = yp + yOffset;
            this.x_um = this.x_um + xOffset;
            this.y_um = this.y_um + yOffset;
            theta = this.uieOffsetAngle.get(); % offset angle;
            delta_um = 2*f_um*tan(asin(lambda_um/T_um));
            switch this.uipHologram.getSelectedIndex()
                case 1
                    % QWLSI
                    this.dI=hologen.utils.getIntensity_QWLSI(xp,yp,delta_um,f_um,lambda_um,theta/180*pi);
                case 2
                    % LSI
                    this.dI=hologen.utils.getIntensity_LSI(xp,yp,delta_um,f_um,lambda_um,theta/180*pi);
            end
            this.dHL = this.dI;
            this.dHL(this.dHL<this.dTh)=0;
            this.dHL(this.dHL>=this.dTh)=1;
            this.dHL = this.dHL.*pupil;
        end
        
        % gds generate function
        function genGDS(this)
            % parameter setting
            f = this.uieFocalLength.get(); % lens focal length
            incidentAngle= this.uieOffsetAngle.get(); % offset angle
            T = this.uieT.get()/1e3; % grating pitch
            db = 10000000; %set unit to anstrom
            lambda = this.uieLambda.get()*1e-6; % wavelength
            delta=2*f*tan(asin(lambda/T));
            NA = this.uieNA.get();
            offset = eval(this.uieOffset.get());
            xOffset = offset(1);
            yOffset = offset(2);
            subR=f*tan(asin(NA)); % lens radius
            R=subR+sqrt(xOffset.^2+yOffset.^2); % parent lens radius
            Nx=500;
            Ny=500;
            obscuration = this.uieObscuration.get(); % central obscuration
            printResolution = this.uieMinFeature.get()*1e-6; % minimum printable feature
            ringSampling=50;
            ringSamplingNum=Nx/ringSampling;
            
            DownSamplingAccuCtrl=0.001; % ratio
            CoordsAccuCtrl=0.0001; % intensity
            filename=fullfile(this.uieFilePath.get(),[this.uieFileName.get(),'.gds']);
            filenamestr=['F',num2str(f),'_T',num2str(T),'_wl',num2str(lambda*1e6)];
            outputFile=hologen.utils.OpenGDS(filename,filenamestr);
            incidentAngle=incidentAngle/180*pi;
            
            RingNum=(sqrt(f^2+R^2)-f)/lambda;
            divnum=2*floor(RingNum/ringSamplingNum)+1;
            dx=2*R/divnum;
            dy=2*R/divnum;
            
            Mx=divnum;
            My=divnum;
            fliped=divnum/2+1;
            shift = eval(this.uieShift.get())*db;
            xShift = shift(1); % coordinates shift in gds
            yShift = shift(2);
            hologram = this.uipHologram.getSelectedIndex();
            
            th = this.dTh; % threshold for generating hololens
            reGen = 0;
            if reGen ~= 1
            tic
            fprintf('Initialing coordinates...\n');
            ceout=cell(Mx,My);
            parfor_progress(Mx);
            parfor p=1:Mx
                parfor_progress;
                for q=1:My
                    ceout{p,q} = hologen.utils.getCoords(lambda,delta,f,R,...
                        subR,xOffset,yOffset,db,Nx,Ny,p,q,Mx,My,dx,dy,fliped,...
                        incidentAngle,th,hologram);
                end
            end
            parfor_progress(0);
            fprintf('Initialing coordinates took %0.1f\n',toc);
            fprintf('Processing coordinates... \n');
            tic,
            Ixy=[];
            uxy=[];
            this.stShapes=[];
            sns=0;
            for p=1:Mx
                for q=1:My
                    if ~isempty(ceout{p,q})
                        this.stShapes=[this.stShapes,ceout{p,q}.Ixy0];
                        Ixy=[Ixy,ceout{p,q}.Ixy];
                        if ~isempty(ceout{p,q}.uxy)
                            ceout{p,q}.uxy(:,4)=ceout{p,q}.uxy(:,4)+sns;
                        end
                        uxy=[uxy;ceout{p,q}.uxy];
                        sns=length(Ixy);
                    end
                end
            end
            % stitching splitted shapes and seperate each shape
            while ~isempty(uxy)
                t=uxy(1,4);
                temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)~=t);
                if isempty(temp)
                    temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)==t);
                    if isempty(temp)
                        break;
                    end
                    %         if length(temp)==1
                    %            uxy(:,:)=uxy([end,1:end-1],:);
                    %            continue;
                    %         end
                    uxy(temp,:)=[];
                    temp3=find(uxy(:,4)==t);
                    if isempty(temp3)
                        try
                            this.stShapes(end+1)=Ixy(t);
                        catch
                            this.stShapes=Ixy(t);
                        end
                    end
                    continue;
                end
                k=uxy(temp,4);
                Ixy(t).xr=[Ixy(t).xr;Ixy(k).xr];
                Ixy(t).yr=[Ixy(t).yr;Ixy(k).yr];
                uxy(temp,:)=[];
                uxy(1,:)=[];
                temp2=find(uxy(:,4)==k);
                if ~isempty(temp2)
                    uxy(temp2,4)=t;
                end
                temp3=find(uxy(:,4)==t);
                if isempty(temp3)
                    this.stShapes(end+1)=Ixy(t);
                end
            end
                      
%             for t=length(this.stShapes):-1:1
%                 if ~all((this.stShapes(t).xr-xOffset).^2+(this.stShapes(t).yr-yOffset).^2>Rb.^2&...
%                         (this.stShapes(t).xr-xOffset).^2+(this.stShapes(t).yr-yOffset).^2<=(subR).^2)
%                     this.stShapes(t)=[];
%                 end
%             end
            fprintf('Processing coordinates took %0.1f\n',toc);
            end
            fprintf('Generating ploygons... \n');
            tic,
            lambda=lambda*db;
            printResolution = printResolution*db;
            delta=delta*db;
            f=f*db;
            xOffset = xOffset*db;
            yOffset = yOffset*db;
            subR=subR*db;
            Rb=obscuration*db;
            len=length(this.stShapes);
            parfor_progress(len);
            for q=1:len
                parfor_progress;
%                 if q~=200
%                     continue;
%                 end
                cxy=[this.stShapes(q).xr,this.stShapes(q).yr];
%                 if q==300
%                 figure(2),plot(cxy(:,1),cxy(:,2),'.'),hold on;drawnow;
%                 end
                cxy=unique(cxy,'rows');
                cn=length(cxy);
%                 k=k+1;
                cx=cxy(:,1);
                cy=cxy(:,2);
                %     if q~=70
                %         continue;
                %     end
                switch hologram
                    case 1
                        B=hologen.utils.CgetBoundariesFromLabelQWLSI(cx,cy,cn,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
                    case 2
                        B=hologen.utils.CgetBoundariesFromLabelLSI(cx,cy,cn,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
                end
%                 Bs = B;
% for si=1:length(B)
%     figure(2),plot(B{si}(:,1),B{si}(:,2),'-');hold on;pause(0.5);
%     
% end
                 Bs=hologen.utils.CdownSamplingUsingRealCoords(B,lambda,delta,f,DownSamplingAccuCtrl);
                
                for si=length(Bs):-1:1
                    
                    while 1 % remove glitch
                        num=length(Bs{si});
                        angle=zeros(num,1);
                        for js=1:num
                            if js==1
                                angle(1)=hologen.utils.calAngle(Bs{si}(end,1),Bs{si}(end,2),Bs{si}(1,1),Bs{si}(1,2),Bs{si}(2,1),Bs{si}(2,2));
                            elseif js==num
                                angle(num)=hologen.utils.calAngle(Bs{si}(num-1,1),Bs{si}(num-1,2),Bs{si}(num,1),Bs{si}(num,2),Bs{si}(1,1),Bs{si}(1,2));
                            else
                                angle(js)=hologen.utils.calAngle(Bs{si}(js-1,1),Bs{si}(js-1,2),Bs{si}(js,1),Bs{si}(js,2),Bs{si}(js+1,1),Bs{si}(js+1,2));
                            end
                        end
                        da=diff([angle;angle(1)]);
                        index=abs(da)>5; % remove glitch, if this value is too small, it may remove correct shapes
                        if sum(abs(angle)>3.1)>1
                            index=index|abs(angle)>3.1;
                        end
                        Bs{si}(index,:)=[];
                        [~,rss]=cart2pol(Bs{si}(:,1),Bs{si}(:,2));
                        if max(rss)-min(rss)<printResolution
                            Bs{si}=[];
                        end
                        if sum(index)==0||length(Bs{si})<3
                            break;
                        end
                    end
                end
                
                
                
                %% generate boundary
                for ns=1:length(Bs)
                    xy=Bs{ns};
                    if isempty(xy)||size(xy,1)<3
                        continue;
                    end
                    if ((xy(1,1)-xOffset)^2+(xy(1,2)-yOffset)^2)>subR^2||((xy(1,1)-xOffset)^2+(xy(1,2)-yOffset)^2)<Rb^2
                        continue;
                    end
                    xy(end+1,:)=xy(1,:);
                    xy(:,1)=xy(:,1)+xShift;
                    xy(:,2)=xy(:,2)+yShift;
                    Np=size(xy,1)-1;
                    hologen.utils.CreateBoundary(outputFile,xy',Np);
                end
            end
            parfor_progress(0);
            fprintf('Generating ploygons took %0.1f\n',toc);
            %% finish GDS
            hologen.utils.CloseGDS(outputFile);
            fprintf('GDS file is generated!\n');
        end
        
        
        
        % validates whether a char edit box evaluates to a Nx2 matrix,
        % colors accordingly.  Empty value is changed to []
        function [lOut, vals] = validateCouplesEditBox(~, src, cDefaultVal)
            lOut = true;
            vals = [];
            if isempty(src.get())
                src.styleDefault();
                src.set(cDefaultVal);
                return
            end
            try
                vals = eval(src.get());
                [~, sc] = size(vals);
                if (sc == 2 || sc == 0)
                    src.styleDefault();
                    lOut = true;
                else
                    src.styleBad();
                    lOut = false;
                end
            catch
                % can't read this edit box
                src.styleBad();
                lOut = false;
            end
        end
        
        % Main redraw function. Pass tab indices to refresh axes
        function replot(this, dTabIdx)
            
            switch dTabIdx
                
                case this.U8SIM
                    imagesc(this.haFieldAmp,this.x_um,this.y_um,this.dHL);
                    axis(this.haFieldAmp,'xy','equal');
                    xlabel(this.haFieldAmp,'x/um'),
                    ylabel(this.haFieldAmp,'y/um');
            end
            
        end
        
        
        
        
        function build(this, hFigure, dOffsetX, dOffsetY)
            if nargin <3
                dOffsetX = 0;
                dOffsetY = 0;
            elseif nargin == 3
                dOffsetY = dOffsetX;
                dOffsetX =  hFigure;
            end
            
            % build the main window
            if nargin == 2||nargin == 4
                this.hFigure = hFigure;
            else
                this.hFigure = figure(...
                    'name', 'GDS propagation GUI v1.200701',...
                    'Units', 'pixels',...
                    'Position', [5 - dOffsetX, 5 - dOffsetY,  this.dWidth, this.dHeight],...
                    'handlevisibility','off',... %out of reach gcf
                    'numberTitle','off',...
                    'Toolbar','none',...
                    'Menubar','none');
            end
            
            
            % Build all containers first:
            drawnow
            
            % Axes
            dTgPx = 20;
            dTgPy = 20;
            this.uitgAxesDisplay.build(this.hFigure, dTgPx, dTgPy, 550, 550);
            
            % Axes:Far field
            uitField = this.uitgAxesDisplay.getTabByName('Simulation');
            
            
            this.haFieldAmp = axes('Parent', uitField, ...
                'Units', 'pixels', ...
                'Position', [80, 75, 400, 400], ...
                'XTick', [], 'YTick', []);
            
            
            
            
            % parameters
            this.hpPara = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Parameters',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 310 260 270] ...
                );
            this.uieLambda.build (this.hpPara,20,20,100,20);
            this.uieT.build (this.hpPara,140,20,100,20);
            this.uieFocalLength.build (this.hpPara,20,60,100,20);
            this.uieNA.build (this.hpPara,140,60,100,20);
            this.uieOffsetAngle.build (this.hpPara,20,100,100,20);
            this.uieObscuration.build (this.hpPara,140,100,100,20);
            this.uieMinFeature.build (this.hpPara,20,140,100,20);
            this.uieFileName.build (this.hpPara,140,140,100,20);
            this.uieOffset.build (this.hpPara,20,180,100,20);
            this.uieShift.build (this.hpPara,140,180,100,20);
            this.uipHologram.build (this.hpPara,20,220,100,20);
            
            % control
            this.hpControl = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Control',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 30 260 270] ...
                );
            this.uieFilePath.build (this.hpControl,20,20,220,20);
            this.uibSetPath.build (this.hpControl,140,60,100,20);
            this.uilFileList.build (this.hpControl,20,83,220,105);
            this.uibGenGDS.build (this.hpControl,140,240,100,20);
            this.uibGenPattern.build (this.hpControl,20,240,100,20);
            this.uibOpenGDS.build (this.hpControl,140,210,100,20);
            
            drawnow;
        end
        
    end
    
end

