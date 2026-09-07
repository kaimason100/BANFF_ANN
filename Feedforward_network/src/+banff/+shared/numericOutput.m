function A = numericOutput(Y)
    if isa(Y, 'dlarray'), A = extractdata(Y); else, A = Y; end
    try, A = gather(A); catch, end
    A = double(A);
end
