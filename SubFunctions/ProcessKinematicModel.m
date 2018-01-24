function Xnext = ProcessKinematicModel(X, dist, phi)
    Xnext = X + [dist*cos(X(3)); dist*sin(X(3)); 0];
    Xnext(3) = phi;
end