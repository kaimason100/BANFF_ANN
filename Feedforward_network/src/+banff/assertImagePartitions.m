function assertImagePartitions(X, XTest, stats, trainIdx, valIdx)
%ASSERTIMAGEPARTITIONS Check the exact representation presented to the network.
A = banff.shared.applyImageStandardization(X,stats);
B = banff.shared.applyImageStandardization(XTest,stats);
A = reshape(single(A),[],size(X,4));
B = reshape(single(B),[],size(XTest,4));
banff.assertNoSharedInputs(A(:,trainIdx).',A(:,valIdx).','image train/validation');
banff.assertNoSharedInputs(A.',B.','image development/official test');
end
