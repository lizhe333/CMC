import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

function argsOf(argv) {
  const out = {};
  for (let i = 2; i < argv.length; i += 1) {
    const token = argv[i];
    if (token === "--inspect-only") out.inspectOnly = true;
    else if (token.startsWith("--")) out[token.slice(2)] = argv[++i];
  }
  return out;
}

function columnName(oneBased) {
  let n = oneBased;
  let s = "";
  while (n > 0) {
    n -= 1;
    s = String.fromCharCode(65 + (n % 26)) + s;
    n = Math.floor(n / 26);
  }
  return s;
}

function nearly(a, b, tol = 1e-9) {
  return Math.abs(a - b) <= tol;
}

async function readFormalRows(csvPath) {
  const text = await fs.readFile(csvPath, "utf8");
  const lines = text.trim().split(/\r?\n/);
  const header = lines[0].split(",");
  const ix = Object.fromEntries(header.map((name, i) => [name, i]));
  const required = [
    "time_s", "position_cm", "domain_in_material",
    "moisture_kg_per_kg", "surface_moisture_kg_per_kg",
  ];
  for (const name of required) {
    if (!(name in ix)) throw new Error(`missing source column: ${name}`);
  }

  const positions = Array.from({ length: 21 }, (_, i) => i / 10);
  const grouped = new Map();
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const cells = line.split(",");
    const t = Number(cells[ix.time_s]);
    const p = Number(cells[ix.position_cm]);
    const key = t.toString();
    if (!grouped.has(key)) grouped.set(key, { time: t, values: new Map(), surfaces: [] });
    const record = grouped.get(key);
    const pIndex = positions.findIndex((x) => nearly(x, p, 1e-8));
    if (pIndex < 0) throw new Error(`unexpected position ${p} at t=${t}`);
    const inside = cells[ix.domain_in_material].toLowerCase() === "true";
    const moisture = inside ? Number(cells[ix.moisture_kg_per_kg]) : null;
    if (inside && !Number.isFinite(moisture)) throw new Error(`invalid moisture at t=${t}, p=${p}`);
    record.values.set(pIndex, moisture);
    record.surfaces.push(Number(cells[ix.surface_moisture_kg_per_kg]));
  }

  const records = [...grouped.values()].sort((a, b) => a.time - b.time);
  if (!records.length || !nearly(records[0].time, 60)) throw new Error("output must begin at 60 s");
  for (let i = 0; i < records.length; i += 1) {
    const r = records[i];
    if (r.values.size !== positions.length) throw new Error(`incomplete positions at t=${r.time}`);
    const surface0 = r.surfaces[0];
    if (!Number.isFinite(surface0) || r.surfaces.some((x) => !nearly(x, surface0, 1e-12))) {
      throw new Error(`inconsistent surface values at t=${r.time}`);
    }
    r.surface = surface0;
    if (i < records.length - 1 && !nearly(r.time, 60 * (i + 1), 1e-8)) {
      throw new Error(`non-minute row before exact endpoint: ${r.time}`);
    }
  }
  const exact = records.at(-1).time;
  if (nearly(exact / 60, Math.round(exact / 60), 1e-10)) {
    throw new Error("exact endpoint unexpectedly duplicates an integer-minute row");
  }
  return {
    positions,
    records,
    rows: records.map((r) => [r.time, ...positions.map((_, i) => r.values.get(i)), r.surface]),
  };
}

const args = argsOf(process.argv);
if (args["verify-only"]) {
  if (!args.source || !args.report) throw new Error("--source and --report are required for verification");
  const candidateBlob = await FileBlob.load(args["verify-only"]);
  const candidate = await SpreadsheetFile.importXlsx(candidateBlob);
  const candidateSheet = candidate.worksheets.getItemAt(0);
  const { positions, records, rows } = await readFormalRows(args.source);
  const lastColumn = columnName(positions.length + 2);
  const lastRow = records.length + 1;
  const actual = candidateSheet.getRange(`A2:${lastColumn}${lastRow}`).values;
  if (actual.length !== rows.length) throw new Error(`row count mismatch: ${actual.length} vs ${rows.length}`);
  let numericMaxAbs = 0;
  let blankMismatch = 0;
  let valueMismatch = 0;
  for (let i = 0; i < rows.length; i += 1) {
    if (actual[i].length !== rows[i].length) throw new Error(`column count mismatch at row ${i + 2}`);
    for (let j = 0; j < rows[i].length; j += 1) {
      const expected = rows[i][j];
      const got = actual[i][j];
      if (expected === null) {
        if (!(got === null || got === undefined || got === "")) blankMismatch += 1;
      } else if (typeof got !== "number" || !Number.isFinite(got)) {
        valueMismatch += 1;
      } else {
        const diff = Math.abs(got - expected);
        numericMaxAbs = Math.max(numericMaxAbs, diff);
        if (diff > 1e-12 * Math.max(1, Math.abs(expected))) valueMismatch += 1;
      }
    }
  }
  const used = candidateSheet.getUsedRange();
  const errors = await candidate.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
    options: { useRegex: true, maxResults: 300 },
    summary: "saved result4 formula error scan",
  });
  const status = blankMismatch === 0 && valueMismatch === 0 && used.address === `A1:${lastColumn}${lastRow}` ? "PASS" : "FAIL";
  if (args.preview) {
    const startRow = Math.max(1, lastRow - 5);
    const endpointPreview = await candidate.render({
      sheetName: candidateSheet.name,
      range: `A${startRow}:${lastColumn}${lastRow}`,
      scale: 1.5,
      format: "png",
    });
    await fs.mkdir(path.dirname(args.preview), { recursive: true });
    await fs.writeFile(args.preview, new Uint8Array(await endpointPreview.arrayBuffer()));
  }
  const report = {
    status,
    workbook: path.resolve(args["verify-only"]),
    source: path.resolve(args.source),
    used_range: used.address,
    expected_used_range: `A1:${lastColumn}${lastRow}`,
    row_count: rows.length,
    column_count: rows[0].length,
    numeric_max_abs_difference: numericMaxAbs,
    blank_mismatch_count: blankMismatch,
    value_mismatch_count: valueMismatch,
    formula_error_scan_ndjson: errors.ndjson,
  };
  await fs.writeFile(args.report, `${JSON.stringify(report, null, 2)}\n`, "utf8");
  process.stdout.write(`${JSON.stringify(report, null, 2)}\n`);
  if (status !== "PASS") process.exitCode = 1;
  process.exit();
}

if (!args.template) throw new Error("--template is required");

const templateBlob = await FileBlob.load(args.template);
const workbook = await SpreadsheetFile.importXlsx(templateBlob);
const sheet = workbook.worksheets.getItemAt(0);

if (args.inspectOnly) {
  const inspection = await workbook.inspect({
    kind: "workbook,sheet,region,computedStyle",
    sheetId: sheet.name,
    range: "A1:W4",
    include: "id,name,values,formulas",
    maxChars: 12000,
  });
  process.stdout.write(`${inspection.ndjson}\n`);
  process.exit(0);
}

if (!args.source || !args.output || !args.preview || !args.report) {
  throw new Error("--source, --output, --preview, and --report are required for authoring");
}

const { positions, records, rows } = await readFormalRows(args.source);
const header = ["时间\\到药材中心的距离", ...positions.map((x) => x.toFixed(1)), "药材表面"];
const lastColumn = columnName(header.length);
const lastRow = rows.length + 1;

const oldUsed = sheet.getUsedRange();
if (oldUsed) oldUsed.clear({ applyTo: "contents" });
sheet.getRange(`A1:${lastColumn}${lastRow}`).values = [header, ...rows];

const all = sheet.getRange(`A1:${lastColumn}${lastRow}`);
all.format.font = { name: "宋体", size: 10, color: "#000000" };
all.format.verticalAlignment = "center";
all.format.horizontalAlignment = "center";
all.format.wrapText = false;
all.format.rowHeight = 18;

const head = sheet.getRange(`A1:${lastColumn}1`);
head.format.font = { name: "宋体", size: 10, bold: false, color: "#000000" };
head.format.rowHeight = 20;

if (lastRow > 2) sheet.getRange(`A2:A${lastRow - 1}`).format.numberFormat = "0";
sheet.getRange(`A${lastRow}`).format.numberFormat = "0.000000000";
sheet.getRange(`B2:${lastColumn}${lastRow}`).format.numberFormat = "0.0000";
sheet.getRange("A:A").format.columnWidth = 24;
sheet.getRange(`B:${lastColumn}`).format.columnWidth = 9;
sheet.getRange(`${lastColumn}:${lastColumn}`).format.columnWidth = 12;
sheet.freezePanes.freezeRows(1);

workbook.recalculate();
const keyCheck = await workbook.inspect({
  kind: "table",
  sheetId: sheet.name,
  range: `A1:${lastColumn}8`,
  include: "values,formulas",
  tableMaxRows: 8,
  tableMaxCols: 23,
  maxChars: 12000,
});
const endpointCheck = await workbook.inspect({
  kind: "table",
  sheetId: sheet.name,
  range: `A${lastRow}:${lastColumn}${lastRow}`,
  include: "values,formulas",
  tableMaxRows: 2,
  tableMaxCols: 23,
  maxChars: 6000,
});
const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 300 },
  summary: "result4 final formula error scan",
});

const preview = await workbook.render({
  sheetName: sheet.name,
  range: `A1:${lastColumn}18`,
  scale: 1.5,
  format: "png",
});
await fs.mkdir(path.dirname(args.preview), { recursive: true });
await fs.writeFile(args.preview, new Uint8Array(await preview.arrayBuffer()));

await fs.mkdir(path.dirname(args.output), { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(args.output);

const report = {
  status: "PASS",
  template: path.resolve(args.template),
  source: path.resolve(args.source),
  output: path.resolve(args.output),
  sheet: sheet.name,
  row_count: records.length,
  column_count: header.length,
  first_time_s: records[0].time,
  exact_endpoint_s: records.at(-1).time,
  minute_rows: records.length - 1,
  exact_endpoint_rows: 1,
  out_of_domain_blank_cells: rows.reduce((n, row) => n + row.slice(1, -1).filter((x) => x === null).length, 0),
  key_check_ndjson: keyCheck.ndjson,
  endpoint_check_ndjson: endpointCheck.ndjson,
  formula_error_scan_ndjson: formulaErrors.ndjson,
};
await fs.writeFile(args.report, `${JSON.stringify(report, null, 2)}\n`, "utf8");
process.stdout.write(`${JSON.stringify(report, null, 2)}\n`);
