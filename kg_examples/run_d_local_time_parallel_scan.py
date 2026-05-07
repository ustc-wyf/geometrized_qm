#!/usr/bin/env python3
"""Parallel short-window runner for the D-branch local-time solver.

This script does not change the physical model.  It only launches independent
solver runs in separate processes, so it is useful for parameter/debug scans.
One single time evolution is still sequential because step n+1 depends on step n.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def parse_csv_floats(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def parse_csv_ints(text: str) -> list[int]:
    return [int(item.strip()) for item in text.split(",") if item.strip()]


def safe_name(value: float | int) -> str:
    return str(value).replace("+", "").replace("-", "m").replace(".", "p")


def load_summary(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def run_one(job: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(job["run_dir"])
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "run.log"
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log_file:
        proc = subprocess.run(
            job["cmd"],
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    seconds = time.perf_counter() - started
    summary = load_summary(run_dir / "summary.json")
    final = summary.get("final", {}) if isinstance(summary, dict) else {}
    params = summary.get("params", {}) if isinstance(summary, dict) else {}
    return {
        "name": job["name"],
        "returncode": proc.returncode,
        "seconds": seconds,
        "sec_per_step": seconds / max(int(job["steps"]), 1),
        "run_dir": str(run_dir.resolve()),
        "log": str(log_path.resolve()),
        "summary": str((run_dir / "summary.json").resolve()),
        "initial_time_old": job["initial_time_old"],
        "dt_old": job["dt_old"],
        "resolution": job["resolution"],
        "steps": job["steps"],
        "stopped_reason": params.get("stopped_reason"),
        "steps_completed": params.get("steps_completed"),
        "disc_min_trusted": final.get("disc_min_trusted"),
        "rho_tilde_min_trusted": final.get("rho_tilde_min_trusted"),
        "local_time_patch_count": final.get("local_time_patch_count"),
        "local_time_uncovered_trusted_measure_fraction": final.get(
            "local_time_uncovered_trusted_measure_fraction"
        ),
        "local_time_n_tau_rel_delta_sum": final.get("local_time_n_tau_rel_delta_sum"),
        "local_time_measure_rel_delta": final.get("local_time_measure_rel_delta"),
    }


def build_jobs(args: argparse.Namespace) -> list[dict[str, Any]]:
    solver = Path(__file__).with_name("simulate_d_tridomain_full_dynamics.py").resolve()
    jobs: list[dict[str, Any]] = []
    combos = itertools.product(args.initial_times_old, args.dt_olds, args.resolutions)
    for index, (initial_time_old, dt_old, resolution) in enumerate(combos):
        name = (
            f"{args.name_prefix}_i{index:03d}"
            f"_t{safe_name(initial_time_old)}"
            f"_dt{safe_name(dt_old)}"
            f"_n{resolution}"
            f"_s{args.steps}"
        )
        run_dir = args.output / name
        cmd = [
            args.python,
            str(solver),
            "--output",
            str(run_dir),
            "--physical-optical",
            "--ell-over-planck",
            str(args.ell_over_planck),
            "--normalize-kg",
            "--initial-time",
            str(initial_time_old),
            "--initial-time-old-units",
            "--geometry-closure",
            "inert_tridomain",
            "--matter-variable-mode",
            "local_time",
            "--local-time-patchwise",
            "--local-time-stationary-candidates",
            "--local-time-norm-floor",
            str(args.local_time_norm_floor),
            "--resolution",
            str(resolution),
            "--dt",
            str(dt_old),
            "--dt-old-units",
            "--steps",
            str(args.steps),
            "--active-dilation",
            str(args.active_dilation),
            "--trusted-erosion",
            str(args.trusted_erosion),
            "--stop-mask",
            args.stop_mask,
            "--no-render",
            "--fast-profile",
            "--diagnostics-every",
            str(args.diagnostics_every),
            "--interface-every",
            str(args.interface_every),
            "--metric-extension-every",
            "0",
            "--geometry-every",
            str(args.geometry_every),
        ]
        if args.quick_diagnostics:
            cmd.append("--quick-diagnostics")
        if args.skip_fields_npz:
            cmd.append("--skip-fields-npz")
        if args.extra_args:
            cmd.extend(args.extra_args)
        jobs.append(
            {
                "name": name,
                "cmd": cmd,
                "run_dir": run_dir,
                "initial_time_old": initial_time_old,
                "dt_old": dt_old,
                "resolution": resolution,
                "steps": args.steps,
            }
        )
    return jobs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run independent D-branch local-time short windows in parallel for faster debugging."
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--name-prefix", default="d_local_time_scan")
    parser.add_argument("--initial-times-old", type=parse_csv_floats, default=parse_csv_floats("8"))
    parser.add_argument("--dt-olds", type=parse_csv_floats, default=parse_csv_floats("2.5e-5"))
    parser.add_argument("--resolutions", type=parse_csv_ints, default=parse_csv_ints("64"))
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--local-time-norm-floor", type=float, default=1.0e-12)
    parser.add_argument("--active-dilation", type=int, default=2)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument("--stop-mask", choices=["support", "trusted"], default="trusted")
    parser.add_argument("--geometry-every", type=int, default=100000)
    parser.add_argument("--diagnostics-every", type=int, default=0)
    parser.add_argument("--interface-every", type=int, default=0)
    parser.add_argument("--quick-diagnostics", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--skip-fields-npz", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Extra arguments appended to each solver command. Put this option last.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(args)
    if args.dry_run:
        summary = {
            "dry_run": True,
            "max_workers": args.max_workers,
            "jobs": [{"name": job["name"], "cmd": job["cmd"], "run_dir": str(job["run_dir"])} for job in jobs],
        }
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return

    started = time.perf_counter()
    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, int(args.max_workers))) as pool:
        future_to_name = {pool.submit(run_one, job): job["name"] for job in jobs}
        for future in concurrent.futures.as_completed(future_to_name):
            result = future.result()
            results.append(result)
            print(json.dumps(result, ensure_ascii=False), flush=True)
    elapsed = time.perf_counter() - started
    results.sort(key=lambda item: item["name"])
    summary = {
        "elapsed_seconds": elapsed,
        "max_workers": args.max_workers,
        "job_count": len(jobs),
        "successful_count": sum(1 for item in results if item["returncode"] == 0),
        "results": results,
    }
    summary_path = args.output / "parallel_scan_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"summary": str(summary_path.resolve()), **summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
