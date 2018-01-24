function [nDetected, Ranges, Bearings, IDs, measuredID] = processLaserData(laserData, Map)
    global Xest;
    global plotHandleObject;
    
    angles = [0:360]'*pi/360;
    cosAngles = cos(angles);
    sinAngles = sin(angles);
    mask13 = uint16(2^13-1);
    d = 0.5; %Offset of the laser scanner
    nDetected = 0;
    IDs = [];
    Ranges = [];
    adjRanges = zeros(length(angles), 1);
    adjAngles = angles;
    Bearings = [];
    measuredID = [];
    
    
    if laserData.n > 0
        allRanges = double(bitand(laserData.Ranges(:,laserData.n),mask13)) * 0.01; % Extracting the ranges and putting them into meters
        Intensity = bitshift(laserData.Ranges(:,laserData.n),-13);
        ii = find(Intensity>0);
        
        for i = 1:length(angles),
            adjRanges(i) = sqrt(d^2 + allRanges(i)^2 - 2*d*allRanges(i)*cos(angles(i) + pi/2));
            temp = acos((d^2+adjRanges(i)^2-allRanges(i)^2)/(2*d*adjRanges(i)));
            if i < 180,
                adjAngles(i) = pi/2 - temp;
            elseif i == 180
                adjAngles(i) = pi/2;
            else
                adjAngles(i) = pi/2 + temp;
            end
            adjAngles(i) = real(adjAngles(i));
        end
        
        cosAngles = cos(adjAngles);
        sinAngles = sin(adjAngles);
%         
%         xx = allRanges .* cosAngles;
%         yy = allRanges .* sinAngles;
         xx = adjRanges .* cosAngles;
         yy = adjRanges .* sinAngles;

%         OOIs = DetectOOIs(allRanges,ii, xx, yy);
        OOIs = DetectOOIs(adjRanges,ii, xx, yy);
        
        high = find(OOIs.Color);
        if ~isempty(high),
            Ranges = sqrt((OOIs.Centers(high, 1)).^2 + (OOIs.Centers(high, 2)).^2);
            %Bearings = atan2(OOIs.Centers(high, 2), OOIs.Centers(high, 1)) - Xest(3) + pi/2;
            Bearings = atan2(OOIs.Centers(high, 2), OOIs.Centers(high, 1));

            theta = Xest(3) - pi/2;
            R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            p = R * [OOIs.Centers(high, 1)'; OOIs.Centers(high, 2)'];
            xxO = p(1, :) + Xest(1);
            yyO = p(2, :) + Xest(2);

            %ii is index of measured ranges
            %uu is index of global co-ords
            [nDetected, ~, ind, mapInd] = dataAssociation(xxO, yyO, Map.x, Map.y);
            IDs = mapInd;
            measuredID = ind;
            set(plotHandleObject.OOIsHandle, 'xdata', xxO, 'ydata', yyO);
%             set(plotHandleObject.OOIsHandle, 'xdata', Map.x(IDs), 'ydata', Map.y(IDs));
        end
    end
end