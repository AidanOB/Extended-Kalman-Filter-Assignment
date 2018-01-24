function r = DetectOOIs(ranges,ii, XX, YY)
    r.N = 0;
    r.Centers = [];
    r.Sizes   = [];
    r.Color = [];
    if length(ii) < 1,
        return;
    end
    
    rangesLength = length(ranges);
    rangesEndIndexes = zeros(rangesLength, 1);
    rangesStartIndexes = zeros(rangesLength, 1);
    rangesStartIndexes(1) = 1;
    %d_thresh= 0.07;
    
    for i = 1:rangesLength - 1,
        d_ri_riP1 = sqrt(ranges(i)^2 + ranges(i+1)^2 - 2*ranges(i)*ranges(i+1)*cosd(0.5));
        d_thresh = 0.05 + abs((ranges(i) - ranges(i+1)) / (ranges(i) + ranges(i+1)))*5;
        if d_ri_riP1 > d_thresh,
            rangesEndIndexes(i) = 1;
            rangesStartIndexes(i+1) = 1;
        end
        
    end
    %rangesEndIndexes
    rEndInd = find(rangesEndIndexes > 0);
    rStartInd = find(rangesStartIndexes > 0);
    rStartInd(end) = [];
    j = 1;
    for i = 1:length(rEndInd),
        xDia = XX(rStartInd(i)) - XX(rEndInd(i));
        yDia = YY(rStartInd(i)) - YY(rEndInd(i));
        diaDist = sqrt(xDia^2 + yDia^2);
        if (diaDist < 0.2) && (diaDist > 0.1),
            r.Centers(j, 1) = mean(XX(rStartInd(i):rEndInd(i)));
            r.Centers(j, 2) = mean(YY(rStartInd(i):rEndInd(i)));
            r.Sizes(j) = diaDist;
            r.N = r.N + 1;
            r.Color(j) = ismember(1, ismember(rStartInd(i):rEndInd(i), ii));
            j = j+1;
        end
        
    end

    %Similar to the DetectPoles.p method
    %Separating the segments for the OOI's
    a = diff(ii);
    inf = zeros(length(a), 1);
    endIndex = find([a inf]>1);
    endIndex = [endIndex; length(ii)];
    k = [true;diff(ii(:))~=1];
    
    s = cumsum(k);
    runLengths =  histc(s,1:s(end));
    startIndex = find(k);
    N = max(s);
    
    %Arranging the positions for the segments
    xList = XX(ii);
    yList = YY(ii);
    Centers = zeros(N, 2);
    Sizes = zeros(N, 1);
    Color =zeros(N,1);

    for i = 1:N,
        %Forms the centers for x and y co-ordinates
        Centers(i, 1) = mean(xList(startIndex(i):endIndex(i)));
        Centers(i, 2) = mean(yList(startIndex(i):endIndex(i)));

        %Finding the diameters
        xDist = xList(startIndex(i)) - xList(endIndex(i));
        yDist = yList(startIndex(i)) - yList(endIndex(i));
        Sizes(i) = sqrt(xDist^2 + yDist^2);
        %Sizes(i) = 0;
        Color(i) = 1;
    end
    if r.N > 1,
        [C, I]  = min(pdist2(r.Centers, Centers));
        for i = length(C):-1:1,
           if C(i) < 0.2
               Centers(i, :) = [];
               Sizes(i) = [];
               Color(i) = [];
           end
        end
    end
    Sizes = Sizes';
    r.Centers = [r.Centers; Centers];
    r.Sizes = [r.Sizes, Sizes];
    r.Color = [r.Color, Color'];

end