function setupPlots()
    global plotHandleObject;
    figure(1); clf; hold on;
    plotHandleObject.globalMapHandle = plot(0, 0, 'k.');
    plotHandleObject.currentMapHandle = plot(0, 0, 'b.');
    plotHandleObject.globalOOIsHandle = plot(0, 0, 'c*');
    plotHandleObject.OOIsHandle = plot(0, 0, 'r*');
    plotHandleObject.RobotPositionHandle = plot(0, 0, 'r^');
    plotHandleObject.RobotHeadingHandle = quiver(0, 0, 0, 0, 'r');
    plotHandleObject.PositionHistory = plot(0, 0, 'g-');

    axis([-4, 6, 0, 10]); grid on;
    hold off;
end