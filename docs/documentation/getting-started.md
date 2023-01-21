# Getting Started

## Build Environment

MicroFC can be built in multiple ways on various operating systems. Please select your desired configuration from the list bellow:

<details>
  <summary><h2>*nix</h2></summary>

- **On supported clusters:** Load environment modules

```console
$ . ./mfc.sh load
```

- **Via [Aptitude](https://wiki.debian.org/Aptitude):**

```console
$ sudo apt update
$ sudo apt upgrade
$ sudo apt install tar wget make cmake gcc g++ \
                   python3 python3-dev         \
                   "openmpi-*" libopenmpi-dev
```

- **Via [Pacman](https://wiki.archlinux.org/title/pacman):**

```console
$ sudo pacman -Syu
$ sudo pacman -S base-devel coreutils  \
                 git ninja gcc-fortran \
                 cmake openmpi python3 \
                 python-pip openssh    \
                 python-virtualenv vim \
                 wget tree
```

If you wish to build MicroFC using [NVidia's NVHPC SDK](https://developer.nvidia.com/hpc-sdk),
first follow the instructions [here](https://developer.nvidia.com/nvidia-hpc-sdk-downloads).

</details>

<details>
  <summary><h2>Windows</h2></summary>

On Windows, you can either use Intel Compilers with the standard Microsoft toolchain or the
[Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/) for a Linux experience.

 <details>
   <summary><h3>Windows + Intel (Native)</h3></summary>

Install the latest version of:
- [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/)
- [Intel® oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
- [Intel® oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)

Then, in order to initialize your development environment, open a terminal window and run:
```console
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```

To follow this guide, please replace `./mfc.sh` with `mfc.bat` when running any
commands. `./mfc.sh` is intended Unix-like systems. You will also have access to the `.sln`
Microsoft Visual Studio solution files for an IDE (Integrated Development 
Environment).

  </details>

  <details>
     <summary><h3>Windows + WSL</h3></summary>

Install the latest version of the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/)
as well as a distribution such as Ubuntu which can be found [here](https://apps.microsoft.com/store/detail/ubuntu/9PDXGNCFSCZV). Acquiring an   interactive session is as simple as typing `wsl` in your
command prompt, or alternatively, selecting the distribution from the dropdown menu
available in the [Microsoft Terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701).

You can now follow the appropriate instructions for your distribution.

  </details>

</details>

<details>
  <summary><h3>MacOS (x86 and Apple Silicon)</h3></summary>

  - **MacOS v10.15 (Catalina) or newer [ZSH]** (Verify with `echo $SHELL`)

```console
$ touch ~/.zshrc
$ open ~/.zshrc
```

  - **Older than MacOS v10.15 (Catalina) [BASH]** (Verify with `echo $SHELL`)
  
```console
$ touch ~/.bash_profile
$ open ~/.bash_profile
```
  
An editor should open. Please paste the following lines into it before saving the file. If you wish to use a version of GNU's GCC other than 11, modify the first assignment. These lines ensure that LLVM's Clang, and Apple's modified version of GCC, won't be used to compile MFC. Further reading on `open-mpi` incompatibility with `clang`-based `gcc` on macOS: [here](https://stackoverflow.com/questions/27930481/how-to-build-openmpi-with-homebrew-and-gcc-4-9). We do *not* support `clang` due to conflicts with our Silo dependency.

```console
# === MFC MPI Installation ===
export MFC_GCC_VER=11
export OMPI_MPICC=gcc-$MFC_GCC_VER
export OMPI_CXX=g++-$MFC_GCC_VER
export OMPI_FC=gfortran-$MFC_GCC_VER
export CC=gcc-$MFC_GCC_VER
export CXX=g++-$MFC_GCC_VER
export FC=gfortran-$MFC_GCC_VER
# === MFC MPI Installation ===
```

**Close the open editor and terminal window**. Open a **new terminal** window before executing the commands bellow.

```console
$ brew install wget make python make cmake coreutils gcc@$MFC_GCC_VER
$ HOMEBREW_MAKE_JOBS=$(nproc) brew install --cc=gcc-$MFC_GCC_VER --verbose --build-from-source open-mpi
```

They will download the dependencies MicroFC requires to build itself. `open-mpi` will be compiled from source, using the version of GCC we specified above with the environment variables `HOMEBREW_CC` and `HOMEBREW_CXX`. Building this package might take a while.

</details>

## Fetching MicroFC

```console
$ git clone https://github.com/MFlowCode/MicroFC.git
$ cd MicroFC
```

## Building MicroFC

MicroFC can be built with support for various (compile-time) features:

| Feature   | Enable    | Disable      | Default | Description                                                     |
| :-------: | :-------: | :----------: | :-----: | --------------------------------------------------------------- |
| **MPI**   | `--mpi`   | `--no-mpi`   | On      | Lets MFC run on multiple processors (and nodes) simultaneously. |
| **GPU**   | `--gpu`   | `--no-gpu`   | Off     | Enables GPU acceleration via OpenACC.                           |
| **Debug** | `--debug` | `--no-debug` | Off     | Requests the compiler build MFC in debug mode.                  |

_⚠️ The `--gpu` option requires that your compiler supports OpenACC for Fortran for your target GPU architecture._

When these options are given to `mfc.sh`, they will be remembered when you issue future commands. You can enable and disable features at any time by passing any of the arguments above. For example, if you have previously built MFC with MPI support and no longer wish to run using MPI, you can pass `--no-mpi` once, for the change to be permanent.

MFC is composed of three codes, each being a separate _target_. By default, all targets (`pre_process`, `simulation`, and `post_process`) are selected. To only select a subset, use the `-t` (i.e `--targets`) argument. For a detailed list of options, arguments, and features, please refer to `./mfc.sh build --help`.

Most first-time users will want to build MFC using 8 threads (or more!) with MPI support:
```console
$ ./mfc.sh build -j 8
```

Examples:

- Build MFC using 8 threads with MPI and GPU acceleration: `$ ./mfc.sh build --gpu -j 8`.
- Build MFC using a single thread without MPI, GPU, and Debug support: `$ ./mfc.sh build --no-mpi`.
- Build MFC's `simulation` code in Debug mode with MPI and GPU support: `$ ./mfc.sh build --debug --gpu -t simulation`.

## Running the Test Suite

Run MFC's test suite with 8 threads:

```console
$ ./mfc.sh test -j 8
```

Please refer to the [Testing](testing.md) document for more information.

## Running an Example Case

MicroFC has example cases in the `examples` folder. You can run such a case interactively using 2 tasks by typing:

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -n 2
```

Please refer to the [Running](running.md) document for more information on `case.py` files and how to run them.
