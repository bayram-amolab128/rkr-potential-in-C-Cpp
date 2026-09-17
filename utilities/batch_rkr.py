#!/usr/bin/env python3
"""
batch_rkr.py

Batch runner for molecular-data files.

For every non-empty, non-comment line in a file list, the script:
  1. Creates an indexed case directory.
  2. Runs Le Roy's executable in its own working directory.
  3. Runs your executable in its own working directory.
  4. Preserves every file produced by each executable.
  5. Captures stdout/stderr and return codes.
  6. Writes a master summary.csv.

This is intentionally designed so programs that always write files with the
same names (for example fort.10, output.dat, RKR.out, etc.) cannot overwrite
outputs from another molecule.

Example
-------
python batch_rkr.py ^
    --filelist molecular_filelist.txt ^
    --leroy "C:\\RKR\\LeRoy-RKR.exe" ^
    --ours "C:\\RKR\\RKR_CPP.exe" ^
    --outdir batch_results

By default BOTH executables receive the molecular data through stdin, i.e.
equivalent to:

    LeRoy-RKR.exe < molecule.dat
    RKR_CPP.exe   < molecule.dat

If one program instead expects the input filename as a command-line argument:

    --ours-mode arg

Then it is run approximately as:

    RKR_CPP.exe input.dat

File-list format
----------------
One input path per line:

I2B.dat
Li2_X.dat
Li2_A.dat

Blank lines and lines beginning with # are ignored.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class RunResult:
    program: str
    executable: str
    mode: str
    returncode: int | None
    elapsed_s: float
    status: str
    workdir: Path
    message: str = ""


def safe_name(text: str) -> str:
    """Return a filesystem-friendly name."""
    out = []
    for ch in text:
        if ch.isalnum() or ch in ("-", "_", "."):
            out.append(ch)
        else:
            out.append("_")
    cleaned = "".join(out).strip("._")
    return cleaned or "case"


def read_filelist(filelist: Path) -> list[Path]:
    """Read molecular input paths from a text file."""
    files: list[Path] = []
    base = filelist.parent.resolve()

    for lineno, raw in enumerate(filelist.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()

        if not line or line.startswith("#"):
            continue

        # Allow quoted paths in the file list.
        if len(line) >= 2 and line[0] == line[-1] and line[0] in ('"', "'"):
            line = line[1:-1]

        p = Path(line).expanduser()

        # Relative paths are interpreted relative to the file-list location.
        if not p.is_absolute():
            p = base / p

        p = p.resolve()

        if not p.is_file():
            raise FileNotFoundError(
                f"File-list line {lineno}: input file does not exist:\n  {p}"
            )

        files.append(p)

    if not files:
        raise ValueError(f"No input files found in {filelist}")

    return files


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8", errors="replace")


def run_program(
    *,
    label: str,
    executable: Path,
    mode: str,
    source_input: Path,
    workdir: Path,
    timeout: float | None,
) -> RunResult:
    """
    Run one executable inside its own isolated directory.

    mode='stdin':
        executable < input.dat

    mode='arg':
        executable input.dat
    """
    workdir.mkdir(parents=True, exist_ok=True)

    # Give each program a local copy so it can freely create/modify nearby files.
    local_input = workdir / source_input.name
    shutil.copy2(source_input, local_input)

    exe = executable.resolve()

    if mode == "stdin":
        command = [str(exe)]
        stdin_handle = local_input.open("rb")
        command_description = f'"{exe}" < "{local_input.name}"'
    elif mode == "arg":
        command = [str(exe), local_input.name]
        stdin_handle = None
        command_description = f'"{exe}" "{local_input.name}"'
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    write_text(workdir / "command.txt", command_description + "\n")
    write_text(workdir / "source_input_path.txt", str(source_input) + "\n")

    stdout_path = workdir / "stdout.txt"
    stderr_path = workdir / "stderr.txt"

    start = time.perf_counter()

    try:
        with stdout_path.open("wb") as stdout_file, stderr_path.open("wb") as stderr_file:
            completed = subprocess.run(
                command,
                cwd=workdir,
                stdin=stdin_handle,
                stdout=stdout_file,
                stderr=stderr_file,
                timeout=timeout,
                check=False,
            )

        elapsed = time.perf_counter() - start
        rc = completed.returncode
        status = "OK" if rc == 0 else "FAILED"
        message = ""

    except subprocess.TimeoutExpired:
        elapsed = time.perf_counter() - start
        rc = None
        status = "TIMEOUT"
        message = f"Process exceeded timeout of {timeout} s"

    except Exception as exc:
        elapsed = time.perf_counter() - start
        rc = None
        status = "ERROR"
        message = f"{type(exc).__name__}: {exc}"

    finally:
        if stdin_handle is not None:
            stdin_handle.close()

    write_text(workdir / "status.txt", status + "\n")
    write_text(workdir / "returncode.txt", "" if rc is None else f"{rc}\n")
    write_text(workdir / "elapsed_seconds.txt", f"{elapsed:.9f}\n")

    if message:
        write_text(workdir / "runner_error.txt", message + "\n")

    return RunResult(
        program=label,
        executable=str(exe),
        mode=mode,
        returncode=rc,
        elapsed_s=elapsed,
        status=status,
        workdir=workdir,
        message=message,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-run Le Roy RKR and a second executable over molecular data files."
    )

    parser.add_argument(
        "--filelist",
        type=Path,
        required=True,
        help="Text file containing one molecular-data filepath per line.",
    )
    parser.add_argument(
        "--leroy",
        type=Path,
        required=True,
        help="Path to Le Roy executable.",
    )
    parser.add_argument(
        "--ours",
        type=Path,
        required=True,
        help="Path to your C++ RKR executable.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("batch_results"),
        help="Directory in which indexed results are stored (default: batch_results).",
    )
    parser.add_argument(
        "--leroy-mode",
        choices=("stdin", "arg"),
        default="stdin",
        help="How Le Roy receives input: stdin or filename argument (default: stdin).",
    )
    parser.add_argument(
        "--ours-mode",
        choices=("stdin", "arg"),
        default="stdin",
        help="How your executable receives input: stdin or filename argument (default: stdin).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Optional timeout in seconds for each executable run.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue processing later molecules if one run fails.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    filelist = args.filelist.expanduser().resolve()
    leroy = args.leroy.expanduser().resolve()
    ours = args.ours.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()

    if not filelist.is_file():
        print(f"ERROR: file list not found: {filelist}", file=sys.stderr)
        return 2

    if not leroy.is_file():
        print(f"ERROR: Le Roy executable not found: {leroy}", file=sys.stderr)
        return 2

    if not ours.is_file():
        print(f"ERROR: your executable not found: {ours}", file=sys.stderr)
        return 2

    try:
        inputs = read_filelist(filelist)
    except Exception as exc:
        print(f"ERROR reading file list: {exc}", file=sys.stderr)
        return 2

    outdir.mkdir(parents=True, exist_ok=True)

    summary_path = outdir / "summary.csv"
    manifest_path = outdir / "manifest.txt"

    manifest_lines = [
        f"filelist={filelist}",
        f"leroy={leroy}",
        f"ours={ours}",
        f"leroy_mode={args.leroy_mode}",
        f"ours_mode={args.ours_mode}",
        f"number_of_inputs={len(inputs)}",
        "",
    ]

    rows: list[dict[str, object]] = []
    overall_failure = False

    print(f"Inputs : {len(inputs)}")
    print(f"Output : {outdir}")
    print()

    for index, input_path in enumerate(inputs, start=1):
        stem = safe_name(input_path.stem)
        case_id = f"{index:04d}_{stem}"
        case_dir = outdir / case_id

        leroy_dir = case_dir / "LeRoy"
        ours_dir = case_dir / "OurProgram"

        case_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, case_dir / input_path.name)
        write_text(case_dir / "source_input_path.txt", str(input_path) + "\n")

        print(f"[{index:04d}/{len(inputs):04d}] {input_path.name}")

        leroy_result = run_program(
            label="LeRoy",
            executable=leroy,
            mode=args.leroy_mode,
            source_input=input_path,
            workdir=leroy_dir,
            timeout=args.timeout,
        )
        print(
            f"    LeRoy      : {leroy_result.status:7s} "
            f"({leroy_result.elapsed_s:.4f} s)"
        )

        ours_result = run_program(
            label="OurProgram",
            executable=ours,
            mode=args.ours_mode,
            source_input=input_path,
            workdir=ours_dir,
            timeout=args.timeout,
        )
        print(
            f"    OurProgram : {ours_result.status:7s} "
            f"({ours_result.elapsed_s:.4f} s)"
        )

        manifest_lines.append(
            f"{case_id}\t{input_path}\t"
            f"LeRoy={leroy_result.status}\tOurProgram={ours_result.status}"
        )

        for result in (leroy_result, ours_result):
            rows.append(
                {
                    "index": index,
                    "case_id": case_id,
                    "input_name": input_path.name,
                    "input_path": str(input_path),
                    "program": result.program,
                    "executable": result.executable,
                    "input_mode": result.mode,
                    "status": result.status,
                    "returncode": "" if result.returncode is None else result.returncode,
                    "elapsed_seconds": f"{result.elapsed_s:.9f}",
                    "workdir": str(result.workdir),
                    "message": result.message,
                }
            )

        case_failed = (
            leroy_result.status != "OK"
            or ours_result.status != "OK"
        )

        if case_failed:
            overall_failure = True
            if not args.continue_on_error:
                print("\nStopping because a run failed.")
                print("Use --continue-on-error if you want all remaining files processed.")
                break

    with summary_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "index",
            "case_id",
            "input_name",
            "input_path",
            "program",
            "executable",
            "input_mode",
            "status",
            "returncode",
            "elapsed_seconds",
            "workdir",
            "message",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    write_text(manifest_path, "\n".join(manifest_lines) + "\n")

    print()
    print(f"Summary written to: {summary_path}")
    print(f"Manifest written to: {manifest_path}")

    if overall_failure:
        print("One or more runs failed; inspect the corresponding stderr.txt/status.txt files.")
        return 1

    print("Batch processing completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
