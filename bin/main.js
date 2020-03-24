#!/usr/bin/env node
const program = require('commander');
const fs = require('fs');
const path = require('path');
const JSZip = require('jszip');
const Docxtemplater = require('docxtemplater');


function fatalIf(test, message, exitCode) {

  if (test) {

    console.log(message);
    process.exit(exitCode);
  }
}

function validateFile(filepath, format, inputFile) {

  fatalIf(
    typeof filepath === 'undefined',
    "FATAL: " + " File has not been defined",
    1
  );
  fatalIf(
    inputFile && ! fs.existsSync(filepath),
    "FATAL: Input " + filepath + " does not exists",
    2
  );
  fatalIf(
    ! inputFile && fs.existsSync(filepath),
    "FATAL: File " + filepath + " already exists",
    3
  );
  fatalIf(
    ! filepath.endsWith(format),
    "FATAL: Input file does not have the expected file extension: " + format,
    4
  );
}


program
  .option('-d, --data <data.json>', 'JSON file containing the data')
  .option('-t, --template <template.docx>', 'The template file the report should be generated from')
  .option('-o, --output <report.docx.>', 'Where the clinical report should be written to')
  .parse(process.argv);

// Validate the command line arguments
validateFile(program.data, "json", true);
validateFile(program.template, "docx", true);
validateFile(program.output, "docx", false);

// Load the JSON data filepath
const jsonData = JSON.parse(fs.readFileSync(program.data));

//Load the docx file as a binary
const zip = new JSZip(
  fs.readFileSync(path.resolve(program.template), 'binary')
);

const doc = new Docxtemplater();
doc.loadZip(zip);
doc.setData(jsonData);

try {
    // render the document (replace all variables)
    doc.render()
}
catch (error) {
    var e = {
        message: error.message,
        name: error.name,
        stack: error.stack,
        properties: error.properties,
    }
    console.log(JSON.stringify({error: e}));
    // The error thrown here contains additional information when logged with JSON.stringify (it contains a property object).
    throw error;
}

const buf = doc.getZip().generate({type: 'nodebuffer'});

// buf is a nodejs buffer, you can either write it to a file or do anything else with it.
fs.writeFileSync(path.resolve(program.output), buf);
