function drawlattice_wsw(Offset, Scaling, limit)
%DRAWLATTICE - Draws the AT lattice to a figure
%  drawlattice(Offset, Scaling)
%
%  Written by Greg Portmann
%  Modified by S.W.Wang on Nov.15,2019

global THERING
if nargin < 3
    limit = length(THERING);
end

% Minimum icon width in meters
MinIconWidth = .09;

if nargin < 1
    Offset = 0;
end
Offset = Offset(1);
if nargin < 2
    Scaling = 1;
end
Scaling = Scaling(1);

SPositions = findspos(THERING, 1:(limit+1));
L = SPositions(end);
plot([0 L], [0 0]+Offset, 'k');

% Remember the hold state then turn hold on
HoldState = ishold;
hold on;

% Make default icons for elements of different physical types
for i = 1:limit
    SPos = SPositions(i);
    if isfield(THERING{i},'BendingAngle') & THERING{i}.BendingAngle
        % make icons for bending magnets
        IconHeight = .3;
        IconColor = [1 1 0];
        IconWidth = THERING{i}.Length;
        if IconWidth < MinIconWidth % meters
            IconWidth = MinIconWidth;
            SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
        end
        vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
        vy = [IconHeight IconHeight -IconHeight -IconHeight];
        h = patch(vx, Scaling*vy+Offset, IconColor,'LineStyle','-');
        set(h, 'UserData', i);

    elseif isfield(THERING{i},'K') & THERING{i}.K
        % Quadrupole
        if THERING{i}.K > 0
            % Focusing quadrupole
            IconHeight = .8;
            IconColor = [1 0 0];
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
            vy = [IconHeight IconHeight 0 0];
        else
            % Defocusing quadrupole
            IconHeight = .8;
            IconColor = [0 0 1];
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                % Center about starting point
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
            vy = [0 0 -IconHeight -IconHeight];
        end
        h = patch(vx, Scaling*vy+Offset, IconColor,'LineStyle','-');
        set(h, 'UserData', i);
       
    elseif isfield(THERING{i},'PolynomB') & length(THERING{i}.PolynomB)>2 & (THERING{i}.PolynomB(3) | any(strcmpi(THERING{i}.FamName,{'SF','SFF','SD','SDD'})))
        
        % Sextupole
        if THERING{i}.PolynomB(3)>0 | any(strcmpi(THERING{i}.FamName,{'SF','SFF'}))
            % Focusing sextupole
            IconHeight = .5;
            IconColor = [1 0 1]; % pink
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos          SPos+.33*IconWidth  SPos+.66*IconWidth  SPos+IconWidth   SPos+IconWidth   SPos+.66*IconWidth  SPos+.33*IconWidth      SPos          SPos];
            vy = [IconHeight/3      IconHeight          IconHeight        IconHeight/3    -IconHeight/3      -IconHeight          -IconHeight     -IconHeight/3  IconHeight/3];
%             vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
%             vy = [IconHeight IconHeight -IconHeight -IconHeight];
        elseif THERING{i}.PolynomB(3)<0
            % Defocusing sextupole
            IconHeight = .5;
            IconColor = [0 1 0]; % green
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos          SPos+.33*IconWidth  SPos+.66*IconWidth  SPos+IconWidth   SPos+IconWidth   SPos+.66*IconWidth  SPos+.33*IconWidth      SPos          SPos];
            vy = [IconHeight/3      IconHeight          IconHeight        IconHeight/3    -IconHeight/3      -IconHeight          -IconHeight     -IconHeight/3  IconHeight/3];
%             vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
%             vy = [IconHeight IconHeight -IconHeight -IconHeight];
        end
        h = patch(vx, Scaling*vy+Offset, IconColor,'LineStyle','-');
        set(h, 'UserData', i);

    elseif isfield(THERING{i},'Frequency') & isfield(THERING{i},'Voltage')
        % RF cavity
        IconColor = [1 0.5 0];
        h = plot(SPos, 0+Offset, 'o', 'MarkerFaceColor', IconColor, 'Color', IconColor, 'MarkerSize', 4);
        set(h, 'UserData', i);
        
    elseif strcmpi(THERING{i}.FamName,'BPM')
        % BPM
        IconColor = 'k';
        h = plot(SPos, 0+Offset, '.', 'Color', IconColor);
        set(h, 'UserData', i);
        
    elseif strcmpi(THERING{i}.FamName,'TV')
        % TV screen
        IconHeight = .7;
        IconColor = [.5 0 0];  %'k';
        h = plot(SPos, IconHeight, 'Marker','Square', 'MarkerFaceColor', IconColor, 'Color', IconColor, 'MarkerSize', 3.5);
        set(h, 'UserData', i);
        
    elseif any(strcmpi(THERING{i}.FamName,{'COR','XCOR','YCOR','HCOR','VCOR'})) | isfield(THERING{i},'KickAngle')
        % Corrector
        IconWidth = THERING{i}.Length;
        IconHeight = 1.0;  % was .8
        vx = [SPos   SPos];
        
        % Draw a line above for a HCM and below for a VCM
        % If it's not in the ML, then draw a line above and below
        CMFound = 1;

        if CMFound
            IconColor = [0 0 0];
            vy = [-IconHeight IconHeight];
            if IconWidth < MinIconWidth
                h = plot(vx, Scaling*vy+Offset, 'Color', IconColor);
            else
                IconWidth = THERING{i}.Length;
                vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
                vy = [IconHeight IconHeight -IconHeight -IconHeight];
                h = patch(vx, Scaling*vy+Offset, IconColor, 'LineStyle', '-');
                if IconWidth < MinIconWidth*2 % meters
                    set(h, 'EdgeColor', IconColor);
                end
            end
            set(h, 'UserData', i);
            CMFound = 0;
        end
        
    end
end
% Leave the hold state as it was at the start
if ~HoldState
    hold off
end
