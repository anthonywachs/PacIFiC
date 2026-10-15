# XML-based Simulation Tests

End-to-end validation tests that run the `grains` binary and check physics against analytical solutions.

## Structure

```
xml/
|-- template.xml          # XML template populated by generate_tests.py
|-- generate_tests.py     # YAML -> XML generator
|-- run_tests.py          # Orchestrates generation, simulation, and plotting
|-- collidingSpheres/     # Two spheres colliding head-on
|-- cylinderWall/         # Cylinder bouncing off a wall at varying angles
|-- movingBox/            # Box in free fall -- time integration convergence study
|-- smoke/                # Minimal smoke tests (fast sanity checks)
|-- results/              # Simulation output (.dat files), created at runtime
|-- plots/                # Generated plots (.eps), created at runtime
\-- logs/                 # Run logs, created at runtime
```

## Adding a New Test

1. Create a subdirectory (e.g. `myTest/`).
2. Add a `myTest.yaml` describing the test matrix (see existing YAMLs for reference).
   - `Constant` -- settings shared across all cases.
   - `Fork` -- list of cases, each with a `name` and `overrides`.
3. Add a `myTest.py` post-processing script that accepts `--result-dir` and `--plot-dir`.

## Running Tests

```bash
# All tests
python3 run_tests.py

# Single test
python3 run_tests.py --case movingBox

# Re-plot without re-running simulations
python3 run_tests.py --case movingBox --plots-only

# Regenerate XMLs only (no simulation)
python3 generate_tests.py -i movingBox/MovingBox.yaml -d movingBox/generated
```

A custom binary path can be set with `--binary` or the `BINARY` environment variable.
