function [phi, dist] = processKinMod(readIMU, readSpeed)
    global CurrAttitude;
    global bias;
    phi = 0;
    dist = 0;
    
    if readIMU.n > 1,
        dtI = zeros(1, readIMU.n);
        time = double(readIMU.tcx(1,:))*0.0001;
        time = time - time(1);
        for i = 1:readIMU.n,
            if i > 1,
                    dtI(i) = time(i) - time(i - 1);
            end
            Gyros = double(readIMU.data(4:6, i));
            Gyros(1) = 0;
            Gyros(2) = 0;
            deltaAngle = IntegrateAttitudeStep(Gyros - bias, dtI(i), CurrAttitude);
            CurrAttitude = CurrAttitude + deltaAngle;
            CurrAttitude(1) = 0;
            CurrAttitude(2) = 0;
            CurrAttitude(3) = mod(CurrAttitude(3)+pi, 2*pi) - pi;
            phi = CurrAttitude(3);
        end
    end
    
    if readSpeed.n > 1,
        dtS = zeros(1, readSpeed.n);
        time = double(readSpeed.tcx(1,:))*0.0001;
        time = time - time(1);
        speed = readSpeed.data(1, :);
        
        for j = 1:readSpeed.n,
            if j > 1,
                dtS(j) = time(j) - time(j - 1);
            end
            deltaDist = IntegrateDistStep(speed(j), dtS(j));
            dist = dist + deltaDist;
        end
    end

end

function [delta] = IntegrateDistStep(speed, dt)
    delta = speed * dt;
end