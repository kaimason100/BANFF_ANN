function dx = dynamicsDerivative(dynamicsName, x)
inputSize = size(x);
if isvector(x)
    x = x(:);
    inputSize = size(x);
end
dx = int_dyn(x, string(dynamicsName), 0, 0, 0, 'Simulate');
assert(isequal(size(dx), inputSize), '%s derivative returned size [%s] for input size [%s].', ...
    dynamicsName, num2str(size(dx)), num2str(inputSize));
end
