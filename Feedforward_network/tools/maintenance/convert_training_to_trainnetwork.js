const fs = require("fs");
const path = require("path");
const os = require("os");
const childProcess = require("child_process");

const repo = process.cwd();

function readMlxCode(relPath) {
  const xml = childProcess.execFileSync("unzip", ["-p", relPath, "matlab/document.xml"], { encoding: "utf8" });
  const match = xml.match(/<!\[CDATA\[([\s\S]*?)\]\]><\/w:t>/);
  if (!match) throw new Error(`No MATLAB code block found in ${relPath}`);
  return match[1].replace(/\r\n/g, "\n");
}

function writeMlxCode(relPath, code) {
  const abs = path.join(repo, relPath);
  const work = fs.mkdtempSync(path.join(os.tmpdir(), "mlx-trainnetwork-"));
  childProcess.execFileSync("unzip", ["-q", abs, "-d", work]);
  const docPath = path.join(work, "matlab", "document.xml");
  let xml = fs.readFileSync(docPath, "utf8");
  const escaped = code.replace(/\r\n/g, "\n").replaceAll("]]>", "]]]]><![CDATA[>");
  xml = xml.replace(/<!\[CDATA\[[\s\S]*?\]\]><\/w:t>/, () => `<![CDATA[${escaped}]]></w:t>`);
  fs.writeFileSync(docPath, xml);
  const backup = `${abs}.bak_trainnetwork`;
  fs.copyFileSync(abs, backup);
  fs.rmSync(abs, { force: true });
  childProcess.execFileSync("zip", ["-q", "-r", abs, ...fs.readdirSync(work)], { cwd: work });
  fs.rmSync(work, { recursive: true, force: true });
  fs.rmSync(backup, { force: true });
}

function removeMetrics(code) {
  return code
    .replace(/,\s*\.\.\.\n\s*'Metrics'\s*,\s*\[[^\]]*\]/g, "")
    .replace(/,\s*\.\.\.\n\s*'Metrics'\s*,\s*["'][^"']+["']/g, "")
    .replace(/,\s*\.\.\.\n\s*Metrics\s*=\s*["'][^"']+["']/g, "")
    .replace(/,\s*'Metrics'\s*,\s*(\[[^\]]*\]|"[^"]*"|'[^']*')/g, "")
    .replace(/,\s*Metrics\s*=\s*("[^"]*"|'[^']*')/g, "");
}

function dedupeBestValidationOutput(code) {
  return code
    .replace(
      /('ValidationFrequency'\s*,\s*10,\s*'OutputNetwork'\s*,\s*'best-validation-loss',\s*\.\.\.\n\s*)'ValidationFrequency'\s*,\s*10,\s*'OutputNetwork'\s*,\s*'best-validation-loss',\s*/g,
      "$1"
    )
    .replace(
      /(ValidationFrequency\s*=\s*10,\s*OutputNetwork\s*=\s*'best-validation-loss',\s*\.\.\.\n\s*)ValidationFrequency\s*=\s*10,\s*OutputNetwork\s*=\s*'best-validation-loss',\s*/g,
      "$1"
    );
}

function addClassificationLayer(code) {
  code = code.replace(/softmaxLayer\];/g, "softmaxLayer\n    classificationLayer];");
  code = code.replace(/softmaxLayer\]/g, "softmaxLayer\n    classificationLayer]");
  code = code.replace(/(\s*softmaxLayer\([^\n]+\))\];/g, "$1\n        classificationLayer];");
  return code;
}

function addRegressionLayer(code) {
  if (code.includes("regressionLayer")) return code;
  return code.replace(/(\n\s*fullyConnectedLayer\([^\n]+\))\];/g, "$1\n    regressionLayer];");
}

function ensureBestValidationOutput(code) {
  code = dedupeBestValidationOutput(code);
  if (code.includes("OutputNetwork")) return code;
  code = code.replace(/('ValidationFrequency'\s*,\s*10)/g, "$1, 'OutputNetwork', 'best-validation-loss'");
  code = code.replace(/(ValidationFrequency\s*=\s*10)/g, "$1, OutputNetwork='best-validation-loss'");
  return code;
}

function commonTrainNetworkConversion(code) {
  code = removeMetrics(code);
  code = code.replace(/net\s*=\s*dlnetwork\(layers\);\n/g, "");
  code = code.replace(/net\s*=\s*dlnetwork\(layerGraph\(layers\)\);\n/g, "");
  code = code.replace(/'Verbose'\s*,\s*false/g, "'Verbose', 0");
  code = ensureBestValidationOutput(code);
  code = dedupeBestValidationOutput(code);
  return code;
}

function fixTabularClassificationScoring(code) {
  const before = code;
  code = code.replace(
    /\[~,\s*idx\]\s*=\s*max\(scores,\s*\[\],\s*1\);\s*\npredictedLabels\s*=\s*categorical\(idx\(:\),\s*1:K,\s*categories\(YTrain\)\);/g,
    "idx = scoresToClassIndex(scores, K, numel(YTest));\npredictedLabels = categorical(idx(:), 1:K, categories(YTrain));"
  );

  if (code !== before && !code.includes("function idx = scoresToClassIndex(")) {
    code = `${code.trimEnd()}\n\nfunction idx = scoresToClassIndex(scores, nClasses, nObs)\n%SCORESTOCLASSINDEX Convert trainNetwork/predict scores to class indices.\nif isdlarray(scores), scores = extractdata(scores); end\nscores = gather(scores);\nif iscell(scores), scores = scores{1}; end\nscores = squeeze(scores);\nif ismatrix(scores) && size(scores, 2) == nClasses && size(scores, 1) == nObs\n    [~, idx] = max(scores, [], 2);\nelseif ismatrix(scores) && size(scores, 1) == nClasses && size(scores, 2) == nObs\n    [~, idx] = max(scores, [], 1);\n    idx = idx(:);\nelse\n    error('Unexpected score size %s for %d classes and %d observations.', mat2str(size(scores)), nClasses, nObs);\nend\nidx = idx(:);\nend\n`;
  }
  return code;
}

function convertTabularClassification(code) {
  code = commonTrainNetworkConversion(addClassificationLayer(code));
  code = code.replace(/\{dlarray\(single\(XVal\),\s*'CB'\),\s*YVal\}/g, "{XVal.', YVal}");
  code = code.replace(/\[net,\s*info\]\s*=\s*trainnet\(dlarray\(single\(XTrain\),\s*'CB'\),\s*YTrain,\s*net,\s*'crossentropy',\s*options\);/g,
    "[net, info] = trainNetwork(XTrain.', YTrain, layers, options);");
  code = fixTabularClassificationScoring(code);
  return code;
}

function convertTabularRegression(code) {
  code = commonTrainNetworkConversion(addRegressionLayer(code));
  code = code.replace(/\{dlarray\(single\(XVal\),\s*'CB'\),\s*YValNorm\}/g, "{XVal.', YValNorm}");
  code = code.replace(/\{dlarray\(single\(XVal\),\s*'CB'\),\s*YTrain\}/g, "{XVal.', YTrain}");
  code = code.replace(/\{dlarray\(single\(XVal\),\s*'CB'\),\s*YVal\}/g, "{XVal.', YVal}");
  code = code.replace(/\[net,\s*info\]\s*=\s*trainnet\(dlarray\(single\(XTrain\),\s*'CB'\),\s*YTrainNorm,\s*net,\s*'mse',\s*options\);/g,
    "[net, info] = trainNetwork(XTrain.', YTrainNorm, layers, options);");
  code = code.replace(/\[net,\s*info\]\s*=\s*trainnet\(dlarray\(single\(XTrain\),\s*'CB'\),\s*YTrain,\s*net,\s*'mse',\s*options\);/g,
    "[net, info] = trainNetwork(XTrain.', YTrain, layers, options);");
  return code;
}

function convertMnist(code) {
  code = commonTrainNetworkConversion(addClassificationLayer(code));
  code = code.replace(/\[net,\s*info\]\s*=\s*trainnet\(XTrain,\s*YTrain,\s*net,\s*'crossentropy',\s*options\);/g,
    "[net, info] = trainNetwork(XTrain, YTrain, layers, options);");
  return code;
}

function convertDynamics(code) {
  code = commonTrainNetworkConversion(addRegressionLayer(code));
  code = code.replace(/N\s*=\s*size\(inpRaw,\s*3\);/g, "N = size(inpRaw, 2);");
  code = code.replace(/MiniBatchSize=size\(XTrain,\s*3\)/g, "MiniBatchSize=size(XTrain, 1)");
  code = code.replace(/sequenceInputLayer\(/g, "featureInputLayer(");
  code = code.replace(/featureInputLayer\(size\(XTrain,\s*1\)\)/g, "featureInputLayer(size(XTrain, 2))");
  code = code.replace(
    /X\s*=\s*dlarray\(inp,\s*'CTB'\);\s*Y\s*=\s*dlarray\(tgt,\s*'CTB'\);\s*XTrain\s*=\s*X\(:,\s*:,\s*1:idx\);\s*YTrain\s*=\s*Y\(:,\s*:,\s*1:idx\);\s*XVal\s*=\s*X\(:,\s*:,\s*idx\+1:end\);\s*YVal\s*=\s*Y\(:,\s*:,\s*idx\+1:end\);/g,
    "XTrain = inp(:, 1:idx).'; YTrain = tgt(:, 1:idx).';\nXVal = inp(:, idx+1:end).'; YVal = tgt(:, idx+1:end).';"
  );
  code = code.replace(/\[net,\s*trainingInfo\]\s*=\s*trainnet\(XTrain,\s*YTrain,\s*net,\s*'mse',\s*options\);/g,
    "[net, trainingInfo] = trainNetwork(XTrain, YTrain, layers, options);");
  return code;
}

function convertPong(code) {
  code = commonTrainNetworkConversion(addClassificationLayer(code));
  code = code.replace(/\[net,\s*tr\]\s*=\s*trainnet\(trainFeatures,\s*trainLabels,\s*net,\s*'crossentropy',\s*options\);/g,
    "[net, tr] = trainNetwork(trainFeatures, trainLabels, layers, options);");
  return code;
}

function convertMotor(code) {
  code = commonTrainNetworkConversion(addRegressionLayer(code));
  code = code.replace(/\[net,\s*tr\]\s*=\s*trainnet\(XTrain,\s*YTrain,\s*net,\s*'mse',\s*options\);/g,
    "[net, tr] = trainNetwork(XTrain, YTrain, layers, options);");
  return code;
}

const files = [
  ...fs.readdirSync("examples/classification").filter((f) => f.endsWith(".mlx")).map((f) => `examples/classification/${f}`),
  ...fs.readdirSync("examples/regression").filter((f) => f.endsWith(".mlx")).map((f) => `examples/regression/${f}`),
  ...fs.readdirSync("examples/dynamical_systems").filter((f) => f.endsWith(".mlx")).map((f) => `examples/dynamical_systems/${f}`),
  "examples/control/motor_control.mlx",
  "examples/pong/Pong.mlx",
];

for (const relPath of files) {
  if (relPath.endsWith("playPong.mlx")) continue;
  let code = readMlxCode(relPath);
  const before = code;
  if (relPath.includes("examples/dynamical_systems/")) {
    code = convertDynamics(code);
  } else if (relPath.includes("examples/pong/")) {
    code = convertPong(code);
  } else if (relPath.includes("examples/control/")) {
    code = convertMotor(code);
  } else if (relPath.includes("examples/regression/")) {
    code = convertTabularRegression(code);
  } else if (relPath.includes("MNIST_") || relPath.includes("mnist_")) {
    code = convertMnist(code);
  } else if (relPath.includes("examples/classification/")) {
    code = convertTabularClassification(code);
  }
  if (code !== before) {
    writeMlxCode(relPath, code);
    console.log(`converted ${relPath}`);
  }
}
