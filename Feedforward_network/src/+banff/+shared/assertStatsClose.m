function assertStatsClose(actual, expected, label, statName)
    actual = double(actual);
    expected = double(expected);
    assert(isequal(size(actual), size(expected)), '%s %s size mismatch.', label, statName);
    tol = 1e-10 * max(1, max(abs(expected(:))));
    assert(all(abs(actual(:) - expected(:)) <= tol), '%s %s was not fitted from the training split.', label, statName);
end
