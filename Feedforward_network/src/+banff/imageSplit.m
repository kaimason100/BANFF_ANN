function [trainIdx,valIdx] = imageSplit(images,labels,seed)
%IMAGESPLIT Preserve the original randperm split when images are distinct.
% Repeat images stay in one group and exact labelled repeats occur only once.
X = reshape(images,[],size(images,4));
[~,keep] = unique([double(X.'), double(labels(:))],'rows','stable');
keep = sort(keep);
[~,~,group] = unique(X(:,keep).','rows','stable');
rng(double(seed),'twister');
order = randperm(max(group));
nTrain = round(0.8*numel(order));
assert(nTrain > 0 && nTrain < numel(order),'Not enough distinct images.');
% Preserve the original shuffled row order (including when there are no duplicates).
groups = accumarray(group,(1:numel(group)).',[],@(v){v});
trainIdx = keep(vertcat(groups{order(1:nTrain)}));
valIdx = keep(vertcat(groups{order(nTrain+1:end)}));
trainIdx = trainIdx(:).'; valIdx = valIdx(:).';
banff.assertPredictorPartitions(X,labels,{trainIdx,valIdx});
end
