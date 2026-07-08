const logBox = document.getElementById("logBox");

function appendLog(message, type = "log") {
    if (!logBox) return;

    const prefix =
        type === "error" ? "[ERROR] " :
        type === "warn"  ? "[WARN]  " :
                           "[INFO]  ";

    logBox.textContent += prefix + message + "\n";
    logBox.scrollTop = logBox.scrollHeight;
}

(function mirrorConsoleToLogBox() {
    const oldLog = console.log;
    const oldWarn = console.warn;
    const oldError = console.error;

    console.log = function(...args) {
        oldLog.apply(console, args);
        appendLog(args.map(String).join(" "), "log");
    };

    console.warn = function(...args) {
        oldWarn.apply(console, args);
        appendLog(args.map(String).join(" "), "warn");
    };

    console.error = function(...args) {
        oldError.apply(console, args);
        appendLog(args.map(String).join(" "), "error");
    };
})();

const G_KEYS = [
    "we","xwe","ywe","zwe","awe","bwe","cwe","dwe",
    "ewe","fwe","gwe","hwe","iwe","jwe","kwe","lwe"
];

const B_KEYS = [
    "Be","ae","ye",
    "_1e","_2e","_3e","_4e","_5e","_6e",
    "_7e","_8e","_9e","_10e","_11e","_12e","_13e"
];

const GENERAL_KEYS = [
    "name","state","kaiser","intmethod",
    "Vmax","space","ladderspace","errortol",
    "m1","m2","netcharge",
    "Te","De","Ke","re"
];
const LABELS = {

    // -----------------------------
    // General
    // -----------------------------
    name: "\\(Molecule\\)",
    state: "\\(State\\)",

    kaiser: "\\(Kaiser\\)",
    intmethod: "\\(Method\\)",

    Vmax: "\\(v_{\\max}\\;(cm^{-1})\\)",
    space: "\\(\\Delta v\\)",
    ladderspace: "\\(ladder\\;step\\)",
    errortol:"\\(Tolerance\\)",
    
    m1: "\\(m_1\\;(amu)\\)",
    m2: "\\(m_2\\;(amu)\\)",
    netcharge: "\\(Charge\\;(e)\\)",
    
    Te: "\\(T_e\\;(cm^{-1})\\)",
    De: "\\(D_e\\;(cm^{-1})\\)",
    Ke: "\\(k_e\\;(cm^{-1}Å^{-2})\\)",
    re: "\\(r_e\\;(Å)\\)",

    // -----------------------------
    // Vibrational coefficients
    // -----------------------------
    we: "\\(\\omega_e\\)",
    xwe: "\\(\\omega_e x_e\\)",
    ywe: "\\(\\omega_e y_e\\)",
    zwe: "\\(\\omega_e z_e\\)",

    awe: "\\(\\omega_e a_e\\)",
    bwe: "\\(\\omega_e b_e\\)",
    cwe: "\\(\\omega_e c_e\\)",
    dwe: "\\(\\omega_e d_e\\)",
    ewe: "\\(\\omega_e e_e\\)",
    fwe: "\\(\\omega_e f_e\\)",
    gwe: "\\(\\omega_e g_e\\)",
    hwe: "\\(\\omega_e h_e\\)",
    iwe: "\\(\\omega_e i_e\\)",
    jwe: "\\(\\omega_e j_e\\)",
    kwe: "\\(\\omega_e k_e\\)",
    lwe: "\\(\\omega_e l_e\\)",

    // -----------------------------
    // Rotational coefficients
    // -----------------------------
    Be: "\\(B_e\\)",
    ae: "\\(\\alpha_e\\)",
    ye: "\\(\\gamma_e\\)",

    _1e: "\\(\\delta_1\\)",
    _2e: "\\(\\delta_2\\)",
    _3e: "\\(\\delta_3\\)",
    _4e: "\\(\\delta_4\\)",
    _5e: "\\(\\delta_5\\)",
    _6e: "\\(\\delta_6\\)",
    _7e: "\\(\\delta_7\\)",
    _8e: "\\(\\delta_8\\)",
    _9e: "\\(\\delta_9\\)",
    _10e: "\\(\\delta_{10}\\)",
    _11e: "\\(\\delta_{11}\\)",
    _12e: "\\(\\delta_{12}\\)",
    _13e: "\\(\\delta_{13}\\)"
};

const DEFAULTS = {

    // General
    name: "Unknown",
    state: "X",

    kaiser: "0",
    intmethod: "Fleming",

    Vmax: "100",
    space: "0.01",
    ladderspace: "3",
    errortol: "1e-9",

    m1: "1.0",
    m2: "1.0",
    netcharge: "0",

    Te: "0.0",
    De: "0.0",
    Ke: "0.0",
    re: "0.0",

    // G(v)
    we:  "0.0",
    xwe: "0.0",
    ywe: "0.0",
    zwe: "0.0",
    awe: "0.0",
    bwe: "0.0",
    cwe: "0.0",
    dwe: "0.0",
    ewe: "0.0",
    fwe: "0.0",
    gwe: "0.0",
    hwe: "0.0",
    iwe: "0.0",
    jwe: "0.0",
    kwe: "0.0",
    lwe: "0.0",

    // B(v)
    Be:  "0.0",
    ae:  "0.0",
    ye:  "0.0",

    _1e:  "0.0",
    _2e:  "0.0",
    _3e:  "0.0",
    _4e:  "0.0",
    _5e:  "0.0",
    _6e:  "0.0",
    _7e:  "0.0",
    _8e:  "0.0",
    _9e:  "0.0",
    _10e: "0.0",
    _11e: "0.0",
    _12e: "0.0",
    _13e: "0.0"
};

[...G_KEYS, ...B_KEYS].forEach(k => DEFAULTS[k] = "0");


function createFields(keys, containerId) {
    const box = document.getElementById(containerId);

    keys.forEach(key => {
        const label = document.createElement("label");
        label.innerHTML = LABELS[key] || key;

        const input = document.createElement("input");
        input.id = "field_" + key;
        input.dataset.key = key;

        box.appendChild(label);
        box.appendChild(input);
    });
}

createFields(G_KEYS, "gvFields");
createFields(B_KEYS, "bvFields");
createFields(GENERAL_KEYS, "generalFields");
MathJax.typesetPromise();


let RKRModule = null;
let generatedDAT = "";

const fileInput      = document.getElementById("fileInput");
const inputText      = document.getElementById("inputText");
const loadButton     = document.getElementById("loadButton");
const generateButton = document.getElementById("generateButton");
const downloadButton = document.getElementById("downloadButton");

createRKRModule().then(Module => {
    RKRModule = Module;
    console.log("RKR WebAssembly loaded.");
});

loadButton.addEventListener("click", () => {
    if (fileInput.files.length === 0) {
        alert("Please choose an input file.");
        return;
    }

    const reader = new FileReader();

    reader.onload = function(e) {
        const text = e.target.result;

        inputText.value = text;              // hidden raw copy
        const params = parseKeyValueText(text);
        fillFieldsFromParams(params);
    };

    reader.readAsText(fileInput.files[0]);
});

generateButton.addEventListener("click", async () => {
    if (RKRModule == null) {
        alert("WebAssembly module is still loading.");
        return;
    }

    const input = buildInputTextFromFields();
    inputText.value = input;

    showLoading();
    generateButton.disabled = true;

    console.log("--------------------------------------");
    console.log("Running RKR...");
    console.log("--------------------------------------");

    await nextFrame();

    try {
        generatedDAT = RKRModule.RunRKRFromText(input);

        if (generatedDAT.startsWith("ERROR")) {
            alert(generatedDAT);
            return;
        }

        plotDAT(generatedDAT);
        console.log("Calculation completed successfully.");
        console.log("Download your Potential Energy Curve (PEC.dat).");
    }
    catch (err) {
        console.error(err);
        alert(err);
    }
    finally {
        hideLoading();
        generateButton.disabled = false;
    }
});

downloadButton.addEventListener("click", () => {
    if (!generatedDAT) {
        alert("Nothing has been generated yet.");
        return;
    }

    const blob = new Blob([generatedDAT], { type: "text/plain" });
    const url = URL.createObjectURL(blob);

    const a = document.createElement("a");
    a.href = url;
    a.download = "PEC.dat";
    a.click();

    URL.revokeObjectURL(url);
});

function parseDAT(datText) {
    const r = [];
    const V = [];

    const lines = datText.split(/\r?\n/);

    for (const line of lines) {
        const clean = line.trim();

        if (clean.length === 0) continue;
        if (clean.startsWith("#")) continue;

        const parts = clean.split(/\s+/);

        if (parts.length < 2) continue;

        const rr = Number(parts[0]);
        const vv = Number(parts[1]);

        if (Number.isFinite(rr) && Number.isFinite(vv)) {
            r.push(rr);
            V.push(vv);
        }
    }

    return { r, V };
}

function plotDAT(datText) {
    const data = parseDAT(datText);

    if (data.r.length === 0) {
        alert("No valid r,V data found to plot.");
        return;
    }

    const trace = {
        x: data.r,
        y: data.V,
        mode: "lines",
        type: "scatter",
        name: "RKR potential"
    };

    const layout = {
        title: {
            text: "RKR Potential Energy Curve",
            font: { size: 14 }
        },
        xaxis: {
            title: "r (Å)",
            zeroline: false
        },
        yaxis: {
            title: "V(r) (cm⁻¹)",
            zeroline: false
        },
        margin: {
            l: 60,
            r: 20,
            t: 45,
            b: 55
        }
    };

    Plotly.newPlot("plot", [trace], layout, {
        responsive: true,
        scrollZoom: true
    });
}


function parseKeyValueText(text) {
    const params = {};

    text.split(/\r?\n/).forEach(line => {
        line = line.split("#")[0].split("%")[0].trim();
        if (!line || !line.includes("=")) return;

        const [key, ...rest] = line.split("=");
        params[key.trim()] = rest.join("=").trim();
    });

    return params;
}

function fillFieldsFromParams(params) {
    const allKeys = [...GENERAL_KEYS, ...G_KEYS, ...B_KEYS];

    // first reset everything
    allKeys.forEach(key => {
        const field = document.getElementById("field_" + key);
        if (field) field.value = DEFAULTS[key] ?? "0";
    });

    // then fill only what the file contains
    Object.entries(params).forEach(([key, value]) => {
        const field = document.getElementById("field_" + key);
        if (field) field.value = value;
    });
}

function buildInputTextFromFields() {
    const allKeys = [...GENERAL_KEYS, ...G_KEYS, ...B_KEYS];

    return allKeys.map(key => {
        const field = document.getElementById("field_" + key);
        let value = field ? field.value.trim() : "";

        if (value === "")
            value = DEFAULTS[key] ?? "0";

        return `${key} = ${value}`;
    }).join("\n");
}


function fillDefaultFields() {
    Object.entries(DEFAULTS).forEach(([key, value]) => {
        const field = document.getElementById("field_" + key);
        if (field && field.value.trim() === "")
            field.value = value;
    });
}

fillDefaultFields();


const loadingOverlay = document.getElementById("loadingOverlay");

function showLoading() {
    loadingOverlay.classList.remove("hidden");
}

function hideLoading() {
    loadingOverlay.classList.add("hidden");
}

function nextFrame() {
    return new Promise(resolve => requestAnimationFrame(resolve));
}


