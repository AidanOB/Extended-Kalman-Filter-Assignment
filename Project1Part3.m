function Project1Part3()
% Part3.m - version modified of DemoEKF.m originally provided by lecturer, MTRN4010 - 2014 
% by AIDAN O'BRIEN – z3324494
% Project 01. Part3
% 27 April 2014
    clear all;
    addpath('./PossumEssentials/');
    addpath('./SubFunctions/');
    global globalMap;
    global Xest;
    global CurrAttitude;
    global bias;
    global plotHandleObject;
    setupPlots();
    
    % Setting up the constants for use in the online implementation
    angles = [0:360]'*pi/360;
    mask13 = uint16(2^13-1);
    d = 0.5;
    %Run the program every 100ms
    DT = 0.1;
    % Defining the stdDeviations of the noises present
    stdDevGyro = 1*pi/180;  % 2 degrees/second , standard deviation of the gyros' noise
    stdDevSpeed = 0.1;      %0.3m/s is a large error
    sdev_rangeMeasurement = 0.25;	% std. of noise in range measurements. 0.25m
    sdev_bearing = 0.5 * pi/180; %In degrees, 0.5 degrees is conservative
    
    API = API_PossumUGV();
    if API.Ok < 1,
        disp('API Initialisation Failed');
        return;     %Breaks out of the program.
    end;

    adjRanges = zeros(length(angles), 1);
    adjAngles = angles;
    m = API.ReadLMS200(2, 1);
    if m.n > 0
    Ranges = double(bitand(m.Ranges(:,m.n),mask13)) * 0.01; % Extracting the ranges and putting them into meters
    for k = 1:length(angles),
        adjRanges(k) = sqrt(d^2 + Ranges(k)^2 - 2*d*Ranges(k)*cos(angles(k) + pi/2));
        temp = acos((d^2+adjRanges(k)^2-Ranges(k)^2)/(2*d*adjRanges(k)));
        if k < 180,
            adjAngles(k) = pi/2 - temp;
        elseif k == 180
            adjAngles(k) = pi/2;
        else
            adjAngles(k) = pi/2 + temp;
        end
        adjAngles(k) = real(adjAngles(k));
    end
    cosAngles = cos(adjAngles);
    sinAngles = sin(adjAngles);
    for j = 1:m.n
    xx(:, j) = adjRanges(:, j) .* cosAngles;
    yy(:, j) = adjRanges(:, j) .* sinAngles;
    end
        else
        disp('Laser data not loaded correctly, please try again...');
        return;
    end
    
    %Rather than bothering with saving and loading, these are the map
    %values
    globalMap.x = [3.7139; 1.3864; 2.4646; -0.1427; -1.3685];
    globalMap.y = [2.2515; 3.0803; 5.7840;  4.5880;  3.2070];
    
    %Plot the map and first scan
    set(plotHandleObject.globalMapHandle, 'xdata', xx, 'ydata', yy);
    set(plotHandleObject.globalOOIsHandle, 'xdata', globalMap.x, 'ydata', globalMap.y);

    %Initial estimated position
    Xest = [ 0; 0; pi/2];
    XestKin = Xest;
    CurrAttitude = [0; 0; pi/2];
    bias = [0; 0; -0.0155];
    
    % Initial co-variance is zero, because known position is perfect.
    P = zeros(3, 3);

    %Covariance of the noise
    Pn = diag([stdDevSpeed * stdDevSpeed, stdDevGyro * stdDevGyro]);
    
    %Uncertainty in the model
    Qi = diag([(0.15)^2 , (0.15)^2, (1*pi/180)^2]);

    cont = 1;
    count = 1;

    while(cont)
        tic;
        readL = API.ReadLMS200(2, 1);
        readI = API.ReadIMU(1, 100);
        readS = API.ReadSpeed(1, 100);
        
        [phi, dist] = processKinMod(readI, readS);
        
        Xest = ProcessKinematicModel(Xest, dist, phi);
        
        %This estimates the covariance of the model after prediction 
        J = [[1, 0, -dist*sin(Xest(3))]; [0, 1, dist*cos(Xest(3))]; [0, 0, 1]]; 
        %Jacobian of noise
        Ju = [[DT*cos(Xest(3)) , 0] ; [DT*sin(Xest(3)) , 0] ; [0, DT]];
        
        %Covariance due to the inputs
        Qu= Ju*Pn*Ju';
        
        %Combined covariance of the model
        Q = Qi + Qu;
        
        %The new covariance
        P = J*P*J'+Q;
        
        % Kalman Filter Prediction step done
        
        %Start Kalman filter from Observations Step
        
        if readL.n > 0
            [nDetectedLandmarks, MeasuredRanges, MeasuredBearings, IDs, mID] = processLaserData(readL, globalMap);
        else
            nDetectedLandmarks = 0;
        end

        for u=1:nDetectedLandmarks,
            
            ID = IDs(u);
            rID = mID(u);

            eDX = (globalMap.x(ID,1)-Xest(1));       % (xu-x)
            eDY = (globalMap.y(ID,1)-Xest(2));       % (yu-y)
            eDD = sqrt((eDX).*(eDX) + (eDY).*(eDY)); %   so : sqrt( (xu-x)^2+(yu-y)^2 ) 
            estDPhi = atan2(eDY, eDX) - Xest(3) + pi/2;
        
            H1 = [-eDX./eDD, -eDY./eDD, 0];   % Jacobian of h(X)
            H2 = [eDY./(eDD.^2), -eDX./(eDD.^2), -1]; %Jacobian of Bearing
            H = [H1; H2];
            
            ExpectedRange = eDD;   % just a coincidence: we already calculated them for the Jacobian, so I reuse it. 
            ExpectedBearing = estDPhi;
        
            % Evaluate residual (innovation)  "Y-h(Xe)"
            z1  = MeasuredRanges(rID) - ExpectedRange;      
            z2  = MeasuredBearings(rID) - ExpectedBearing;
            z2 = mod(z2+pi, 2*pi) - pi;
            z = [z1; z2];

            
            R1 = (sdev_rangeMeasurement * sdev_rangeMeasurement) * 4;
            R2 = (sdev_bearing * sdev_bearing) * 4;%[0, sdev_bearing * sdev_bearing];
            R = [R1, 0; 0, R2];
        
            % Some intermediate steps for the EKF (as in the lecture notes)
            S = R + H*P*H';
            iS = inv(S);
            K = P*H'*iS;

            % ----- do it...I am getting  X(k+1|k+1) and P(k+1|k+1)
            Xest = Xest+K*z;

            P = P-P*H'*iS*H*P;     % update the Covariance

        % -----  individual EKF update done ...
        end; 
        
        set(plotHandleObject.RobotPositionHandle, 'xdata', Xest(1), 'ydata', Xest(2));
        set(plotHandleObject.RobotHeadingHandle, 'xdata', Xest(1), 'ydata', Xest(2), 'udata', cos(Xest(3)), 'vdata', sin(Xest(3)));
        Xhistory(1, count) = Xest(1);
        Xhistory(2, count) = Xest(2);
        count = count +1;
        set(plotHandleObject.PositionHistory, 'xdata', Xhistory(1, :), 'ydata', Xhistory(2, :));
                
        pause(DT-toc);
%        pause();
        cont = API.Ok;
        
    end


end