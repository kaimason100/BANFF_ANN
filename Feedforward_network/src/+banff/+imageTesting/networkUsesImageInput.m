function tf = networkUsesImageInput(net)
    tf = false;
    if isprop(net, 'Layers') && ~isempty(net.Layers)
        layerClass = lower(class(net.Layers(1)));
        tf = contains(layerClass, 'imageinput');
    end
end
