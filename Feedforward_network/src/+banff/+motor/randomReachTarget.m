function target = randomReachTarget(l1, l2)
    r = 0.5 + (l1 + l2 - 0.5) * rand;
    angle = 2*pi*rand;
    target = [r*cos(angle); r*sin(angle)];
end
