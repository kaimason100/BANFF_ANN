function X = ensure4d(X)
    sz = size(X);
    if numel(sz) == 3 && sz(3) ~= 1
        X = reshape(X, 28, 28, 1, sz(3));
    elseif numel(sz) == 2 && sz(1) == 28*28
        X = reshape(X, 28, 28, 1, sz(2));
    elseif numel(sz) == 4
    else
        error('Unsupported image array shape.');
    end
end
