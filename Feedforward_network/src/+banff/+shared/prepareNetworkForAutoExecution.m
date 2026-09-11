function [net, useGPU, executionDevice] = prepareNetworkForAutoExecution(net)
%PREPARENETWORKFORAUTOEXECUTION Move a dlnetwork to an available GPU or CPU.
assert(isa(net, 'dlnetwork'), 'Expected a dlnetwork.');

useGPU = false;
executionDevice = "CPU";
try
    useGPU = canUseGPU();
catch
    try
        gpuDevice();
        useGPU = true;
    catch
        useGPU = false;
    end
end

if useGPU
    net = moveNetwork(net, @gpuArray);
    try
        device = gpuDevice();
        executionDevice = "GPU: " + string(device.Name);
    catch
        executionDevice = "GPU";
    end
end
end

function net = moveNetwork(net, moveFcn)
learnables = net.Learnables;
for k = 1:numel(learnables.Value)
    value = learnables.Value{k};
    if isa(value, 'dlarray')
        learnables.Value{k} = dlarray(moveFcn(extractdata(value)));
    else
        learnables.Value{k} = moveFcn(value);
    end
end
net.Learnables = learnables;
end
