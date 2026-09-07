function phi = encodeJointError(eJoint)
    phi = [sin(eJoint(1)); cos(eJoint(1)); sin(eJoint(2)); cos(eJoint(2)); eJoint(3); eJoint(4)];
end
