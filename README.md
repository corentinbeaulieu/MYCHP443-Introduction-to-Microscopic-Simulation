# Microscopic Simulation

The repository contains my version of the Microscopic Simulation project given by Michel Masella at CHPS master, UVSQ.

The code is written in [zig](https://www.ziglang.org). Please make sure you have the [`0.12.0` release](https://github.com/ziglang/zig/releases/tag/0.12.0) version of the compiler
in order to build this project.

## Compilation

```sh
$ zig build -Doptimize=ReleaseSafe
```

The executable is named `simulation` and located in `zig-out/bin` directory.

## Usage

```
$ ./simulation [option] input_position_file [input_speed_file] [output_file]
options:
    -h, --help
            Display this help and exit

    -d, --delta-t <TIME>
            Choose delta t value (default: 1 fs)

    -t, --total-iterations <TIME>
            Number of iterations to do (default: 100)

    -s, --save-step <N>
            Save positions every N step (default: 10)

    -m, --thermostat-step <N>
            Frequency to update the temperature using Beredsen thermostat (default: 1)
```

The input files must be formatted as follow
```
first line
mass x y z
    .
    .
    .
```
For reference, please look at `particule.xyz` file.

## Results

During execution, statistics about the simulated system are printed (Energies, Temperature, Sum of forces...).
A color scheme is present to differentiate abnormal values from satisfying ones.
At this end, a summary of the density observed during the simulation is printed (exam question).

Finally, the program outputs the system in a PDB format (default file: test.pdb) on a regular basis.
This file can then be used to visualize the simulation with [vmd](https://www.ks.uiuc.edu/Research/vmd/) for instance.

## Credit
All the guidelines and tips to have a working simulation code were given by Michel Masella
