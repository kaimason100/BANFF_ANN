function dist = endEffectorDistance(xState, target, l1, l2)
    th1 = xState(1); th2 = xState(2);
    endEffector = [l1*cos(th1) + l2*cos(th1 + th2); l1*sin(th1) + l2*sin(th1 + th2)];
    dist = hypot(endEffector(1) - target(1), endEffector(2) - target(2));
end
