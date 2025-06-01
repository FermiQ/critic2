# `critic2.F90`

## Overview

`critic2.F90` is the main program file for the critic2 application. It serves as the entry point and controls the overall execution flow. Its primary responsibilities include initializing various modules, parsing command-line arguments, setting up the user interface (either GUI or text-based), and managing the main execution loop. It also handles the finalization and cleanup tasks before the program exits.

## Key Components

### Program

*   **`critic`**: The main program unit.
    *   **Behavior**:
        1.  Initializes parameters and timing.
        2.  Parses command-line arguments using `stdargs` from `tools_io`.
        3.  Initializes global settings (`global_init`), system module (`systemmod_init`), and Spglib hall mapping.
        4.  Parses global control options (e.g., `testing`, `quiet`, `help`).
        5.  Prints an initial banner and configuration information unless in `quiet` mode.
        6.  If GUI is enabled and requested (`usegui` is true and `HAVE_GUI` is defined), it starts the GUI using `gui_start` from `gui_main`.
        7.  Otherwise, it runs the text-based interface by calling `critic_main` from the `global` module.
        8.  Performs cleanup operations: `grid1_clean_grids`, `history_end`, `global_end`.
        9.  Prints a success message with warning and comment counts and timing information unless in `quiet` mode.
        10. Optionally pauses execution on Windows if not using the GUI to allow users to see the output.

## Important Variables/Constants

*   **`optv`** (`character(len=:), allocatable`): Stores command-line options.
*   **`ghome`** (`character(len=:), allocatable`): Stores the GFORTRAN_HOME path or equivalent.
*   **`fileroot`** (from `global` module, `character(len=mlen)`): The base name for input/output files.
*   **`quiet`** (from `global` module, `logical`): If true, suppresses most screen output. Set by the `-q` option or if `usegui` is true.
*   **`testing`** (from `global` module, `logical`): If true, enables testing mode. Set by the `-t` option.
*   **`usegui`** (from `tools_io` module, `logical`): If true, indicates that the GUI should be used.

## Usage Examples

The program is typically run from the command line:

```bash
critic2 [options] <input_file>
```

Common options:
*   `-h`: Display help message.
*   `-q`: Quiet mode (suppress output).
*   `-t`: Enable testing mode.

The input file (specified by `fileroot`) contains commands and data for the critic2 application.

## Dependencies and Interactions

*   **Uses Modules**:
    *   `gui_main` (conditionally, if `HAVE_GUI` is defined): Provides `gui_start` for launching the graphical user interface.
    *   `spglib`: Provides `spg_build_hall_mapping` for space group symmetry calculations.
    *   `systemmod`: Provides `systemmod_init` for initializing system-related functionalities.
    *   `grid1mod`: Provides `grid1_clean_grids` for cleaning up 1D grid data.
    *   `global`: Provides various global variables (`fileroot`, `quiet`, `testing`), banner/help functions (`initial_banner`, `help_me`, `config_write`), core logic (`critic_main`), and initialization/finalization routines (`global_init`, `global_end`).
    *   `config`: Provides `getstring` and `istring_datadir` for accessing configuration data, likely the data directory path.
    *   `tools_io`: Provides utilities for I/O (`uout`, `ucopy`), command-line argument parsing (`stdargs`), history tracking (`history_init`, `history_end`), timing (`start_clock`, `print_clock`, `tictac`), string conversion (`string`), and GUI flag (`usegui`). Also provides `ncomms`, `nwarns` for counting comments and warnings.
    *   `param`: Provides `param_init` for parameter initialization.
*   **Preprocessor Directives**:
    *   `#ifdef HAVE_GUI ... #endif`: Conditionally compiles GUI-related code.
    *   `#ifdef _WIN32 ... #endif`: Conditionally compiles code for Windows to pause execution.

```
