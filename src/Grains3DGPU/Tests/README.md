# Tests

Unit tests using Google Test. Cover geometry primitives (Vector3, Matrix3, Quaternion), collision detection (GJK, OBB, OBC), rotation math consistency, and contact force tables.

## Prerequisites

Build the Grains library first:
```bash
cd Grains && make
```

## Building and running

```bash
cd Tests
mkdir build && cd build
source ../../Env/grainsGPU.env.sh
cmake ..
make -j$(nproc)
./grains_tests
```

Or from the repo root:
```bash
make build-tests && make run-tests
```

## XML-based tests

End-to-end simulation tests live under `xml/`. See `xml/README.md`.
