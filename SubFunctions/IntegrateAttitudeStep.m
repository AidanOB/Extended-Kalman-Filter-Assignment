function attitudeAdjustment = IntegrateAttitudeStep(gyros, dt, currAttitude)
    %gyros  : vector with the gyros measurement scaled to rad/s, for a
    %single change
    %dt     : time step between vectors
    %currAttitude : the initial attitude that is being calculated from

    % Sorting the inputs into easy to read variables
    ang = currAttitude;
    wx  = gyros(1);
    wy  = gyros(2);
    wz  = gyros(3);
    
    % Performing a single trig function for each variable that is used
    % multiple times
    cosang1 = cos(ang(1));
    cosang2 = cos(ang(2));
    sinang1 = sin(ang(1));

    % Integration for roll pitch and yaw
    roll    = dt * (wx + (wy * sinang1 +wz * cosang1) * tan(ang(2)));
    pitch   = dt * (wy * cosang1 - wz * sinang1);
    yaw     = dt * ((wy * sinang1 + wz * cosang1) / cosang2);

    % Returning the newAttitude
    attitudeAdjustment = [roll; pitch; yaw];
end