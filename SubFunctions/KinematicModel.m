function Xnext = KinematicModel(X, dist, phi)
    Xnext = X + [dist*cos(X(3)); dist*sin(X(3)); phi];
end