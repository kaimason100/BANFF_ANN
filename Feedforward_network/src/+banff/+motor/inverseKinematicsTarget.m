function xTargetJoint = inverseKinematicsTarget(target, l1, l2)
    x = target(1); y = target(2);
    c2 = (x^2 + y^2 - l1^2 - l2^2) / (2*l1*l2);
    c2 = max(min(c2, 1), -1);
    theta2 = acos(c2);
    k1 = l1 + l2*cos(theta2); k2 = l2*sin(theta2);
    theta1 = atan2(y, x) - atan2(k2, k1);
    xTargetJoint = [theta1; theta2; 0; 0];
end
