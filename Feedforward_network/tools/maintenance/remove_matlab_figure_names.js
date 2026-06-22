const fs = require("fs");
const path = require("path");
const os = require("os");
const childProcess = require("child_process");

const repo = process.cwd();

function listMlx(dir) {
  const out = [];
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.isDirectory()) out.push(...listMlx(p));
    if (entry.isFile() && entry.name.endsWith(".mlx")) out.push(p);
  }
  return out;
}

function stripFigureName(xml) {
  let next = xml;
  next = next.replace(/figure\(\s*'Name'\s*,\s*sprintf\([^)]*\)\s*,\s*/g, "figure(");
  next = next.replace(/figure\(\s*'Name'\s*,\s*sprintf\([^)]*\)\s*\)/g, "figure");
  next = next.replace(/figure\(\s*'Name'\s*,\s*'[^']*'\s*,\s*/g, "figure(");
  next = next.replace(/figure\(\s*'Name'\s*,\s*'[^']*'\s*\)/g, "figure");
  return next;
}

function updateMlx(file) {
  const work = fs.mkdtempSync(path.join(os.tmpdir(), "mlx-edit-"));
  childProcess.execFileSync("unzip", ["-q", file, "-d", work]);
  const docPath = path.join(work, "matlab", "document.xml");
  const original = fs.readFileSync(docPath, "utf8");
  const updated = stripFigureName(original);
  if (updated === original) {
    fs.rmSync(work, { recursive: true, force: true });
    return false;
  }
  fs.writeFileSync(docPath, updated);
  const tempOut = path.join(os.tmpdir(), `${path.basename(file)}.${Date.now()}.tmp`);
  childProcess.execFileSync("zip", ["-qr", tempOut, "."], { cwd: work });
  fs.copyFileSync(tempOut, file);
  fs.rmSync(tempOut, { force: true });
  fs.rmSync(work, { recursive: true, force: true });
  return true;
}

let count = 0;
for (const root of ["examples", "src", "tests"]) {
  const rootPath = path.join(repo, root);
  if (!fs.existsSync(rootPath)) continue;
  for (const file of listMlx(rootPath)) {
    if (updateMlx(file)) count += 1;
  }
}

console.log(`Updated ${count} live scripts.`);
