#!/usr/bin/env python3
"""
Run XML tests under Tests/xml/* using the built grains binary,
then run any Python plotting scripts in each test directory.

This mirrors the behaviour of run_xml_tests.sh but implemented in Python.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# ANSI color helpers
_C_RESET = "\033[0m"
_C_GREEN = "\033[32m"
_C_RED = "\033[31m"

def _green(s: str) -> str:
    return f"{_C_GREEN}{s}{_C_RESET}"

def _red(s: str) -> str:
    return f"{_C_RED}{s}{_C_RESET}"


def find_binary(default_bin):
    bin_env = os.environ.get("BINARY")
    if bin_env:
        return Path(bin_env)
    return Path(default_bin)


def main():
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    bin_dir = "bin" + os.environ.get("GRAINS_FULL_EXT", "")
    default_bin = repo_root / "Main" / bin_dir / "grains"
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", "-b", help="Path to grains binary (overrides default)")
    parser.add_argument("--timeout", "-t", type=int, default=300, help="Per-simulation timeout (seconds)")
    parser.add_argument("--generator", help="Path to generate_tests.py", default=str(script_dir / "generate_tests.py"))
    parser.add_argument("--xml-root", help="XML tests root directory", default=str(script_dir))
    parser.add_argument("--clean", action="store_true", help="Remove all generated/ directories under the xml root and exit")
    parser.add_argument("--plots-only", action="store_true", help="Skip generation and simulation; only run plotting scripts using existing results")
    parser.add_argument("--case", help="Run tests only for this subdirectory (name under xml-root or path)")
    args = parser.parse_args()

    xml_root = Path(args.xml_root)
    if not xml_root.exists():
        print(f"XML root not found: {_red(str(xml_root))}")
        sys.exit(1)

    # If requested, remove all generated directories and the top-level logs/results/plots
    if args.clean:
        removed_generated = 0
        for sub in sorted(xml_root.iterdir()):
            if not sub.is_dir():
                continue
            if sub.name in ("logs", "plots", "results"):
                continue
            gen_dir = sub / "generated"
            if gen_dir.exists() and gen_dir.is_dir():
                try:
                    shutil.rmtree(gen_dir)
                    print(f"Removed: {_green(str(gen_dir))}")
                    removed_generated += 1
                except Exception as e:
                    print(f"Failed to remove {str(gen_dir)}: {_red(str(e))}", file=sys.stderr)

        # remove top-level directories
        removed_roots = 0
        for root_name in ("logs", "plots", "results"):
            root_dir = xml_root / root_name
            if root_dir.exists() and root_dir.is_dir():
                try:
                    shutil.rmtree(root_dir)
                    print(f"Removed: {_green(str(root_dir))}")
                    removed_roots += 1
                except Exception as e:
                    print(f"Failed to remove {str(root_dir)}: {_red(str(e))}", file=sys.stderr)

        print(f"Done. Removed {removed_generated} generated/ directories and {removed_roots} root dirs.")
        return

    binary = Path(args.binary) if args.binary else find_binary(default_bin)
    if not binary.exists() or not os.access(str(binary), os.X_OK):
        print(f"Error: grains binary not found or not executable at {_red(str(binary))}")
        sys.exit(1)

    # roots
    logs_root = xml_root / "logs"
    plots_root = xml_root / "plots"
    results_root = xml_root / "results"
    for p in (logs_root, plots_root, results_root):
        p.mkdir(parents=True, exist_ok=True)

    tests_run = tests_passed = tests_failed = 0

    # determine which subdirectories to run
    if args.case:
        case_path = Path(args.case)
        if case_path.is_dir():
            candidates = [case_path]
        else:
            candidate = xml_root / args.case
            if candidate.is_dir():
                candidates = [candidate]
            else:
                print(f"Case not found: {args.case}")
                sys.exit(1)
    else:
        candidates = [p for p in sorted(xml_root.iterdir()) if p.is_dir()]

    for sub in candidates:
        if not sub.is_dir():
            continue
        name = sub.name
        if name in ("logs", "plots", "results"):
            continue

        print(f"Running: {_green(name)}")

        # If requested, only run plotting scripts for this subdirectory
        if args.plots_only:
            print(f"  {_green('plots-only')}: running plots for existing results")
            test_log_dir = logs_root / name
            test_plots_dir = plots_root / name
            test_results_dir = results_root / name
            test_log_dir.mkdir(parents=True, exist_ok=True)
            test_plots_dir.mkdir(parents=True, exist_ok=True)
            if not test_results_dir.exists():
                print(f"  {_red('MISSING RESULTS DIR')} - expected: {test_results_dir}")
                continue

            # run plotting scripts once for the directory
            for py in sorted(sub.glob("*.py")):
                pybase = py.name
                plot_log = test_log_dir / f"{pybase}.plot.log"
                try:
                    with plot_log.open("wb") as out:
                        cmd = [sys.executable, py.name, "--result-dir", str(test_results_dir), "--plot-dir", str(test_plots_dir)]
                        subprocess.run(cmd, cwd=str(sub), stdout=out, stderr=subprocess.STDOUT, check=True)
                        print(f"  {pybase}: {_green('PLOT GENERATED')}")
                except subprocess.CalledProcessError:
                    print(f"  {pybase}: {_red('FAILED')} (log: {plot_log})", file=sys.stderr)
            continue

        # Clean any existing generated xmls
        gen_dir = sub / "generated"
        gen_dir.mkdir(exist_ok=True)
        for f in gen_dir.glob("*.xml"):
            try:
                f.unlink()
            except Exception:
                pass

        # choose a YAML config if present (prefer .yaml/.yml). If none, generator runs with defaults.
        yamls = sorted(sub.glob("*.yaml")) + sorted(sub.glob("*.yml"))
        if yamls:
            cfg = yamls[0]
        else:
            # JSON support removed: if JSON exists, warn and fall back to defaults
            jsons = sorted(sub.glob("*.json"))
            if jsons:
                print(f"  {jsons[0].name}: {_red('JSON config found but JSON support removed; using defaults instead')}")
            cfg = None

        test_log_dir = logs_root / name
        test_log_dir.mkdir(parents=True, exist_ok=True)
        gen_log = test_log_dir / "generate.log"

        # run the generator
        gen_cmd = [sys.executable, args.generator]
        if cfg is not None and cfg.suffix.lower() in (".yaml", ".yml"):
            gen_cmd += ["--input", str(cfg), "--out-dir", str(gen_dir)]
        else:
            gen_cmd += ["--out-dir", str(gen_dir)]

        with gen_log.open("wb") as out:
            try:
                subprocess.run(gen_cmd, cwd=str(xml_root), stdout=out, stderr=subprocess.STDOUT, check=False)
            except Exception as e:
                out.write(str(e).encode())

        # Prefer XMLs matching the input stem when a YAML config was used,
        # otherwise take any generated XMLs.
        if cfg is not None and cfg.suffix.lower() in (".yaml", ".yml"):
            xmls = sorted(gen_dir.glob(f"{cfg.stem}*.xml"))
            if not xmls:
                print(f"  generated/*.xml: {_red('XML MISSING')} (generator log: {gen_log})")
                continue
        else:
            xmls = sorted(gen_dir.glob("*.xml"))
            if not xmls:
                print(f"  generated/*.xml: {_red('XML MISSING')} (generator log: {gen_log})")
                continue

        for xml in xmls:
            xmlbase = xml.name
            print(f"  {xmlbase}: {_green('XML CREATED')} (log: {gen_log})")
            tests_run += 1

            test_log_dir = logs_root / name
            test_plots_dir = plots_root / name
            test_results_dir = results_root / name
            test_log_dir.mkdir(parents=True, exist_ok=True)
            test_plots_dir.mkdir(parents=True, exist_ok=True)
            test_results_dir.mkdir(parents=True, exist_ok=True)

            log_file = test_log_dir / f"{xmlbase}.log"

            try:
                with log_file.open("wb") as out:
                    subprocess.run([str(binary), str(xml)], stdout=out, stderr=subprocess.STDOUT, timeout=args.timeout, check=True)
                print(f"  {xmlbase}: {_green('SIMULATION FINISHED')}")
                tests_passed += 1
            except subprocess.CalledProcessError:
                print(f"  {xmlbase}: {_red('FAIL')} (log: {log_file})")
                tests_failed += 1
                continue
            except subprocess.TimeoutExpired:
                print(f"  {xmlbase}: {_red('FAIL (timeout)')} (log: {log_file})")
                tests_failed += 1
                continue

        # After all XMLs for this subdirectory have been simulated, run plotting
        # scripts once for the directory (so plots aggregate results from all cases).
        for py in sorted(sub.glob("*.py")):
            pybase = py.name
            plot_log = test_log_dir / f"{pybase}.plot.log"
            try:
                with plot_log.open("wb") as out:
                    cmd = [sys.executable, py.name, "--result-dir", str(test_results_dir), "--plot-dir", str(test_plots_dir)]
                    subprocess.run(cmd, cwd=str(sub), stdout=out, stderr=subprocess.STDOUT, check=True)
                    print(f"  {pybase}: {_green('PLOT GENERATED')}")
            except subprocess.CalledProcessError:
                print(f"  {pybase}: {_red('FAILED')} (log: {plot_log})", file=sys.stderr)

    print("\n==========================================")
    print("  Test Summary")
    print("==========================================")
    print(f"Run:    {tests_run}")
    print(f"Passed: {_green(str(tests_passed))}")
    print(f"Failed: {_red(str(tests_failed))}")


if __name__ == "__main__":
    main()
