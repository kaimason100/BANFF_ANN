function seed = splitSeed(task)
%SPLITSEED Preserve the published dataset seeds, independent of weight seed.
if strcmp(char(task),'abalone'), seed = 1; else, seed = 42; end
end
