function schedule = publicationSchedule()
%PUBLICATIONSCHEDULE Exact feedforward task/seed/epoch/learning-rate schedule.
% Grouped splitting is the deliberate methodological difference in this branch.
tasks = ["iris","breast_cancer","car_quality","mushroom","abalone","toyota", ...
    "afro_mnist_ethiopic","afro_mnist_nko","afro_mnist_osmanya","afro_mnist_vai", ...
    "kmnist","mnist_fashion","mnist","Pong","LQR_two_link_arm", ...
    "Lorenz","MO0","MO5","MO7","MO13","Rikitake","SprottB","SprottC","SprottS"];
schedule = table('Size',[0,6], ...
    'VariableTypes',{'string','double','string','double','double','logical'}, ...
    'VariableNames',{'Task','Seed','Stage','Epochs','LearnRate','DerivativeField'});
for task = tasks
    isDS = any(task == tasks(16:end));
    for seed = 0:9
        if isDS
            epochs = 100000; lr = 0.01;
            if (any(task == ["MO5","MO13"]) && seed == 3) || ...
                    (task == "MO7" && seed == 8) || ...
                    (task == "Rikitake" && any(seed == [0,4]))
                lr = 0.005;
            end
        else
            cfg = banff.settings(task);
            epochs = cfg.MaxEpochs; lr = cfg.LearnRate;
        end
        schedule(end+1,:) = {task,seed,"initial",epochs,lr,isDS};
        if task == "MO0"
            schedule(end+1,:) = {task,seed,"continuation",50000,0.001,true};
        end
    end
end
end
