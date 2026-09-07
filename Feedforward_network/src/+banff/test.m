function test(task,seeds)
%TEST Shared tabular/image evaluation with saved-split and duplicate checks.
if nargin < 2, seeds = 0:9; end
validateattributes(seeds,{'numeric'},{'vector','integer','>=',0,'<=',9,'nonempty'});
task = char(task);
if any(strcmp(task,{'iris','breast_cancer','car_quality','mushroom'}))
    banff.testClassification(task,seeds);
elseif any(strcmp(task,{'abalone','toyota'}))
    banff.testRegression(task,seeds);
elseif any(strcmp(task,{'mnist','mnist_fashion','kmnist','afro_mnist_ethiopic','afro_mnist_nko','afro_mnist_osmanya','afro_mnist_vai'}))
    banff.testImages(task,seeds);
else
    error('Unknown dataset task: %s. Use the dedicated closed-loop test Live Script.',task);
end
end
