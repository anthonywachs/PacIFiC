# Welcome to PacIFiC!
[![Copr build status](https://copr.fedorainfracloud.org/coprs/colive/PacIFiC/package/pacific/status_image/last_build.png)](https://copr.fedorainfracloud.org/coprs/colive/PacIFiC/package/pacific/)

PacIFiC is a high-performance MPI parallel c/c++ software to compute particle-laden flows at the particle scale. PacIFiC stands for "PArtiCles In FluId Computations".

PacIFiC is open-sourced under the MIT license, and is developed by the research group of Prof. Anthony Wachs at the University of British Columbia, Vancouver, Canada with the support of IFP Energies nouvelles, France.

```mermaid
graph TD;
  pacific-->Octree;
  pacific-->Cartesian;
  pacific-->Grains3D;
  Octree-.->basilisk;
  Octree-->FDM1[DLMFD];
  Octree-->LBM;
  Octree-->MovingCutCell;
  Octree-->VTKHyperTree;
  Octree-->eulerian_caps;
  Octree-->lagrangian_caps
  Cartesian-->FLUID;
  FLUID-->FDM2[DLMFD];
  FLUID-->DirectionSplitting;
  Cartesian-->MacWorld;
  MacWorld-->MAC;
  MacWorld-.->HYPRE;
  MacWorld-.->PETSc;
```


## Install

For those who wish to *just use* our software, the easiest and most hassle-free method by far is to use our linux packages. For those who wish to contribute to PacIFiC, instructions on locally building the suite yourself follow in [building](#building)

### Fedora

On an RPM-based distro you may install from our [copr repo](https://copr.fedorainfracloud.org/coprs/colive/PacIFiC/):
```bash
sudo dnf copr enable colive/PacIFiC 
sudo dnf install -y pacific-openmpi pacific-openmpi-devel
```
Installation on Enterprise Linux requires the `copr` plugin in addition to CodeReady Builder (CRB) and EPEL release repositories, if they have not already been enabled.

Since RPM-based distros require installation under `/usr/lib64/openmpi/` or similar, you will need to run
```bash
source /etc/profile.d/modules.sh
module load mpi/openmpi-x86_64
```
Then our MPI-linked libraries and binaries become availible under `$PATH`, `$LD_LIBRARY_PATH` as well as the include path from `mpicc --showme`. You should now be able to run ``grains3d --version``. 

### Ubuntu

On Ubuntu, you may install from our [ppa archive](https://launchpad.net/~colive/+archive/ubuntu/pacific).
```bash
sudo apt-add-repository ppa:colive/pacific
sudo apt update
sudo apt install pacific-basilisk pacific-tools libpacific-mpi-dev
```

## Building 

### Requirements

The following depenencies are required:

 * [OpenMPI](https://www.open-mpi.org/)
 * [Xerces-C++](https://xerces.apache.org/xerces-c/)
 * [zlib](https://www.zlib.net/)

In addition, the following toolchain is required or reccomended
 * [GCC](https://gcc.gnu.org/)  
 * [CMake](https://cmake.org/)
 * [Make](https://www.gnu.org/software/make/)

On RPM-based distributions (e.g. RedHat, Fedora) these can be obtained using 
```bash
sudo dnf install -y gcc g++ make cmake git pkgconfig openmpi-devel xerces-c-devel zlib-ng-devel hdf5-openmpi-devel petsc-openmpi-devel
```
On apt-based distributions (e.g. Ubuntu, Debian), use
```bash
sudo apt-get install -y build-essential make cmake git pkg-config libopenmpi-dev libxerces-c-dev zlib1g-dev libhdf5-openmpi-dev petsc-dev
```

### Building with CMake

The PacIFiC project requires the CMake meta-build tool and one actual build tool. 

First, configure and generate the build system by running

```bash
cmake --preset Default
```

To build all targets in the project

```bash
cmake --build --preset Default
```

To build just one target in the project (for example, Grains) run 

```bash
cmake --build --preset Default --target Grains3D_main
```

Once the `--build` step is run, the resulting executable `grains3d` may be found in `build/release/src/Grains3D/Main` which itself links against the shared libraries such as `libGrains3D.so` that also live in the build tree. 

You may run the executable using the absolute or relative path e.g. by typing 

```bash
./build/release/src/Grains3D/Main/grains3d
```

from the project root. But typing `grains3d` alone *does not work* since the build folder is not part of the `$PATH` variable that your system uses to find binaries. If you are someone who is developing Grains3D (and frequently rebuilding) you may find it useful to not type this everytime. In this case, a file is generated that may be sourced as 

```bash
source ./build/release/pacific-env.sh
```
which temporarily adds the build folder to `$PATH`. Once your shell is closed (or you type `exit`) then `grains3d` will no longer run without typing the full path.

If you are not frequently rebuilding and wish to install the entire project to a system directory, where `grains3d` will be found in `$PATH` every type, then you may install your build to a system directory with

```bash
sudo cmake --install build/release 
```
which by default chooses `/usr` as the installation prefix. Since `/usr/bin` is found in path, you may once again use it from anywhere.

### Building with Make

To build the project with make, first source the configuration file with 

```bash
source env/config.env.sh
```

Next, build the project 

```bash
make all
```

Individual targets (e.g. `grains3d` or `fluid`) may be built with 
```bash
make grains3d
```

Once they are built, the resulting binaries may be found in `build/PacIFiC-{$PACIFIC_ARCH}`. So long as you have run `source` in your current shell session, they will be added and findable via your `$PATH` variable and thus you should now be able to simply type e.g.
```bash
grains3d --help
```

At the moment there is no `install` step for the make build system: In order to use the software you must always first run `source`.


## Examples

Some examples may be found in `examples`: These can be relocated, and by default, use `FetchContent` to include copies of PacIFiC into their own build tree automatically. They therefore do not require installation of PacIFiC, although they still assume the project dependencies are installed.

## Building on HPC

### Sockeye

CMake presets are provided to help build on `sockeye.arc.ubc.ca`. It is reccomended to build inside of the compute nodes, which are offline. Since these do not have internet access, we must use git submodules to recursively clone the submodules in `third_party` to ensure they are available at configure and build time. 

Start by cloning the submodule into scratch
```bash
cd ~/scratch/user
git clone https://github.com/anthonywachs/PacIFiC.git --recurse-submodules
```
Next, start an interactive shell 
```bash
srun --pty -A st-wachs-1 -p interactive_cpu -t 00:10:00 -N 1 -c 8 --mem=16G bash -l
```
Inside the compute node, run
```bash
module load gcc cmake
cmake --preset SockeyeRelease
cmake --build --preset SockeyeRelease
cmake --install build/sockeye_release
exit
```
Now that you are off the compute node, the installed files (including XercesC) now live inside `~/scratch/user/PacIFiC/install/sockeye_release`. You may relocate this if you'd like with 
```bash
mv ~/scratch/user/PacIFiC ~/project/user
```
By default, the binaries and libraries are not yet availible on `$PATH` or `$LD_LIBRARY_PATH`. You may either `export` these manually to add them, or alternatively, write a `lmod` file such as 
```lua
whatis("Name : PacIFiC")
whatis("Version : 0.0.1")
whatis("Target : skylake_avx512")
whatis("Short description : Particles In Fluid Computations ")
help([[Name   : PacIFiC]])
help([[Version: 0.0.1]])
help([[Target : skylake_avx512]])
help([[.]])
depends_on("gcc/9.4.0")
depends_on("zlib-ng/2.0.7")
depends_on("openmpi/4.1.1-cuda11-3")
depends_on("hdf5/1.10.7-additional-bindings")
prepend_path{"PATH","/home/user/project/user/PacIFiC/install/sockeye_release/bin",delim=":"}
prepend_path{"CMAKE_PREFIX_PATH","/home/user/project/user/PacIFiC/install/sockeye_release/.",delim=":"}
append_path{"LD_LIBRARY_PATH","/home/user/project/user/PacIFiC/install/sockeye_release/lib64",delim=":"}
setenv("UBC_CLUSTER","sockeye")
```
and saving this to `$HOME/project/user/modules/PacIFiC/0.0.1.lua`. This can then be used in slurm jobs by calling 
```bash
module use "$HOME/project/user/modules"
module load PacIFiC/0.0.1
```

## Documentation

Additional documentation may be found by visiting the [documentation website](https://anthonywachs.github.io/PacIFiC/).
