# Multi-Echo GRE Pulseq for MATLAB (gre-pulseq-matlab)

This is an advanced MRI sequence development project based on MATLAB and the [Pulseq](https://github.com/pulseq/pulseq) framework. The repository is specifically designed to generate and export multi-echo 3D Gradient Echo (GRE) pulse sequences with high flexibility and hardware safety.

## Relationship to Shuffle3DGRE & What's New

`gre-pulseq-matlab` was originally inspired by and built upon the excellent [Shuffle3DGRE](https://github.com/rita-schmidt/Shuffle3DGRE) project by Rita Schmidt. While inheriting its advanced K-space shuffling and trajectory reordering concepts, this repository introduces architectural refactoring and various new extensions.

**Key Enhancements in this Repository:**
* **Highly Modularized Architecture:** The sequence generation is fully refactored into a strict modular pipeline (`prep/` and `check/` modules), significantly improving readability, maintainability, and reusability for custom sequence modifications.
* **Multi-Echo Readout Support:** Introduces dedicated readout strategies (Bipolar, Monopolar, and Separate TRs).
* **Accelerated Imaging:** Expanded support for parallel imaging (e.g., 2D CAIPIRINHA acceleration masks).

## Core Sequence Features

The codebase supports various spatial encoding and readout modes to cater to different scanning requirements:

* **Readout Modes:**
  * **Bipolar:** Alternating gradient polarities for minimum ESP.
  * **Monopolar:** Flyback gradients inserted between echoes to enforce identical sampling directions for all echoes.
  * **Separate:** Each echo has an independent RF excitation (only one echo is acquired per TR).
* **3D Spatial Encoding & Undersampling:** Supports flexible sampling orders in the PE (Y) and 3D (Z) dimensions, including standard Cartesian ordering, CAIPIRINHA patterns with acceleration factors, and random Shuffling trajectories.

## Dependencies

Before running the code, please ensure you have the following installed and added to your MATLAB path:
* **MATLAB**
* [**Pulseq**](https://github.com/pulseq/pulseq) (MATLAB core library)
* [**Shuffle3DGRE**](https://github.com/RitaSchmidt/Shuffle3DGRE) (Optional: Required only if using specific shuffling sampling trajectories).

*(Note: Dependencies can usually be initialized via Git submodules. Refer to the `.gitmodules` file in the root directory).*

## Quick Start

1. **Set up paths**: Ensure all folders in the root directory and the `pulseq` library are added to your MATLAB search path.
2. **Configure and run the main script**: Open `GRE_3D.m` (or `GRE_3D_Separate.m`) and modify the parameters in the `Setup` struct according to your experimental needs:
   ```matlab
   Setup.nRO = 40;
   Setup.nPE = 40;
   Setup.n3D = 20;
   Setup.TE = [2.1, 3.6, 5.1]*1e-3;
   Setup.TR = 25e-3;
   Setup.bBipolarROGrads  = false; % Set to true for bipolar, false for monopolar
3. **Export the sequence**: Upon running the script, the code will automatically generate the sequence, run timing/safety checks, and output a .seq file ready to be executed on the scanner.
