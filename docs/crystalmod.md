# `crystalmod.f90`

## Overview

The `crystalmod.f90` module is a cornerstone of the critic2 application, providing the fundamental tools and data structures for representing and manipulating crystal and molecular structures. It defines the `crystal` derived type, which encapsulates all structural information, including lattice parameters, atomic coordinates, symmetry operations, and molecular fragments. The module and its associated submodules offer a comprehensive suite of procedures for various crystallographic computations, structure editing, environment analysis, symmetry determination, and input/output operations for a wide range of file formats.

Submodules extend the functionality:
*   `crystalmod@proc.f90`: Handles initialization, creation, and basic manipulation of `crystal` objects, including coordinate transformations.
*   `crystalmod@env.f90`: Manages atomic environments, neighbor lists, and distance calculations.
*   `crystalmod@mols.f90`: Deals with identifying and analyzing molecular fragments within crystal structures.
*   `crystalmod@complex.f90`: Implements more complex operations like Ewald sums, promolecular density calculations, and packing analysis.
*   `crystalmod@edit.f90`: Provides routines for modifying crystal structures, such as cell transformations, atom deletion/movement, and reordering.
*   `crystalmod@symmetry.f90`: Contains procedures for determining and applying symmetry operations, site symmetry, and space group information, often interfacing with Spglib.
*   `crystalmod@powderproc.F90`: Calculates powder diffraction patterns and radial distribution functions.
*   `crystalmod@vibrations.F90`: Handles molecular and crystal vibrations, including reading force constants and calculating phonon frequencies.
*   `crystalmod@write.f90`: Responsible for writing crystal and molecular structures to various file formats and generating reports.

## Key Components

### Derived Types/Classes

*   **`vibrations`**: Represents molecular or crystal vibrations.
    *   **Fields/Attributes**:
        *   `file` (`character(len=mlen)`): Source file of vibration data.
        *   `ivformat` (`integer`): Format of the vibration data.
        *   `hasfc2` (`logical`): True if second-order force constants (FC2) are available.
        *   `fc2` (`real*8, allocatable(:,:,:,:)`): 2nd-order FC matrix (3,3,nat,nat).
        *   `isinit` (`logical`): True if frequencies/eigenvectors are available.
        *   `nqpt` (`integer`): Number of q-points.
        *   `qpt` (`real*8, allocatable(:,:)`): q-point coordinates (3,nqpt) (fractional).
        *   `nfreq` (`integer`): Number of frequencies.
        *   `freq` (`real*8, allocatable(:,:)`): Frequencies (nfreq,nqpt) (cm-1).
        *   `vec` (`complex*16, allocatable(:,:,:,:)`): Phonon eigenvector (3,nat,nfreq,nqpt).
    *   **Procedures**: `end`, `print_summary`, `print_fc2`, `print_freq`, `print_eigenvector`, `read_file`, `calculate_q`, `calculate_thermo`.

*   **`crystal`**: Represents a crystal or molecular structure.
    *   **Key Fields/Attributes**:
        *   `isinit` (`logical`): True if the crystal structure has been initialized.
        *   `havesym` (`integer`): Indicates if symmetry was determined (0: no, 1: full).
        *   `file` (`character(len=mlen)`): Source file name.
        *   `isformat` (`integer`): File format identifier.
        *   `nspc` (`integer`): Number of species.
        *   `spc` (`type(species), allocatable(:)`): Array of species information.
        *   `nneq` (`integer`): Number of non-equivalent atoms.
        *   `at` (`type(neqatom), allocatable(:)`): Array of non-equivalent atoms.
        *   `ncel` (`integer`): Number of atoms in the main cell.
        *   `atcel` (`type(celatom), allocatable(:)`): List of atoms in the main cell.
        *   `aa`, `bb` (`real*8(3)`): Cell lengths (bohr) and angles (degrees).
        *   `omega` (`real*8`): Unit cell volume.
        *   `gtensor`, `grtensor` (`real*8(3,3)`): Metric and reciprocal metric tensors.
        *   `m_x2c`, `m_c2x`, etc. (`real*8(3,3)`): Various crystallographic/Cartesian conversion matrices.
        *   `spgavail` (`logical`): True if Spglib symmetry information is available.
        *   `spg` (`type(SpglibDataset)`): Spglib's symmetry dataset.
        *   `neqv`, `ncv` (`integer`): Number of symmetry operations and centering vectors.
        *   `cen` (`real*8, allocatable(:,:)`): Centering vectors.
        *   `rotm` (`real*8(3,4,48)`): Symmetry operations.
        *   `ismolecule` (`logical`): True if the system is a molecule in a box.
        *   `molx0`, `molborder` (`real*8(3)`): Centering vector and border for molecular systems.
        *   `nstar` (`type(neighstar), allocatable(:)`): Neighbor stars (connectivity).
        *   `nmol` (`integer`): Number of molecules in the unit cell.
        *   `mol` (`type(fragment), allocatable(:)`): Molecular fragments.
        *   `iperiod` (`integer`): Periodicity identifier (e.g., 0D, 1D, 2D, 3D).
        *   `vib` (`type(vibrations)`): Molecular/crystal vibrations data.
    *   **Key Bound Procedures (many more exist, categorized by submodule)**:
        *   `init`, `end`, `struct_new` (from `@proc`): Construction, destruction, initialization.
        *   `x2c`, `c2x`, `distance`, `shortest` (from `@proc`): Basic crystallographic operations.
        *   `build_env`, `list_near_atoms`, `nearest_atom` (from `@env`): Atomic environments.
        *   `fill_molecular_fragments`, `calculate_periodicity` (from `@mols`): Molecular fragment analysis.
        *   `ewald_energy`, `promolecular_array3` (from `@complex`): Complex calculations.
        *   `newcell`, `cell_standard`, `delete_atoms` (from `@edit`): Structure modification.
        *   `sitesymm`, `calcsym`, `clearsym` (from `@symmetry`): Symmetry operations.
        *   `powder`, `rdf` (from `@powderproc`): Powder diffraction and RDF.
        *   `report`, `write_cif`, `write_vasp` (from `@write`): Output and reporting.

*   **`xrpd_peaklist`**: Represents X-ray powder diffraction (XRPD) peak information.
    *   **Key Fields/Attributes**:
        *   `npeak` (`integer`): Number of peaks.
        *   `lambda` (`real*8`): Wavelength of the radiation (Angstrom).
        *   `fpol` (`real*8`): Polarization correction factor.
        *   `th2ini`, `th2end` (`real*8`): Initial and final reflection angles (degrees).
        *   `th2` (`real*8, allocatable(:)`): Reflection angles (2-theta, degrees).
        *   `ip` (`real*8, allocatable(:)`): Peak intensities.
        *   `hvec` (`integer, allocatable(:,:)`): Reflection indices (h,k,l).
    *   **Procedures**: `end`, `from_crystal`, `from_peaks_file`, `from_profile_file`, `write`, `calculate_profile`.

### Module Procedures (Standalone)

*   `search_lattice`: (from `@symmetry`) Searches for lattice points.
*   `pointgroup_info`: (from `@symmetry`) Gets holohedry and Laue class from H-M symbol.
*   `crosscorr_gaussian`: (from `@powderproc`) Calculates Gaussian cross-correlation between XRPD patterns.
*   `vcpwdf_compare`: (from `@powderproc`) Compares crystals allowing for cell deformations.
*   `gaussian_compare`: (from `@powderproc`) Compares a crystal to an XRPD pattern using GPWDF.
*   `david_sivia_calculate_background`: (from `@powderproc`) Calculates background for a diffraction pattern.

## Important Variables/Constants

*   **`iperiod_3d_crystal`, `iperiod_2d`, `iperiod_mol_single`, etc.** (`integer, parameter`): Define system periodicity types (e.g., 3D crystal, 2D slab, single molecule).
*   **`iperiod_vacthr`** (`real*8, parameter`): Threshold for vacuum detection (15.0 Bohr).
*   **`holo_unk`, `holo_tric`, `holo_mono`, etc.** (`integer, parameter`): Holohedry identifiers.
*   **`holo_string`** (`character(len=12), parameter, public(0:7)`): String representations of holohedries.
*   **`laue_unk`, `laue_1`, `laue_2m`, etc.** (`integer, parameter`): Laue class identifiers.
*   **`laue_string`** (`character(len=12), parameter, public(0:11)`): String representations of Laue classes.
*   **`xrpd_lambda_def`, `xrpd_fpol_def`, `xrpd_sigma_def`, etc.** (`real*8, parameter`): Default parameters for powder X-ray diffraction routines.

## Usage Examples

```fortran
use crystalmod
use crystalseedmod, only: crystalseed ! For creating an initial structure
implicit none

type(crystal) :: my_crystal
type(crystalseed) :: my_seed

! ... (Initialize my_seed with atomic positions, cell parameters, etc.) ...

! Create a new crystal structure from the seed
call my_crystal%struct_new(my_seed, crashfail=.true.)

! Check if initialization was successful
if (my_crystal%isinit) then
    ! Perform operations, e.g., calculate cell volume
    print *, "Cell volume (Bohr^3): ", my_crystal%omega

    ! Write the structure to a CIF file
    call my_crystal%write_cif("my_structure.cif", usesym0=.true.)

    ! Calculate powder pattern (list of peaks)
    ! ... (define th2ini, th2end, lambda_val, fpol_val) ...
    ! call my_crystal%powder(mode=1, th2ini0=th2ini, th2end0=th2end, &
    !                        lambda0=lambda_val, fpol=fpol_val, &
    !                        th2p=my_th2_peaks, ip=my_intensities)
end if

! Clean up
call my_crystal%end()
```

## Dependencies and Interactions

*   **Uses Modules**:
    *   `spglib`: For space group symmetry determination (`SpglibDataset`, `spg_standardize_cell`, etc.).
    *   `types`: Provides basic data types like `neqatom`, `celatom`, `neighstar`, `species`, `thread_info`.
    *   `fragmentmod`: Provides the `fragment` type for representing molecular fragments.
    *   `param`: Provides global parameters like `maxzat0` (max atomic number), `mlen` (max string length).
    *   `crystalseedmod`: Used by `struct_new` to initialize a `crystal` from a `crystalseed`.
    *   `grid1mod`: Used for atomic density grids in `struct_new` (via `grid1_register_ae`).
    *   `global`: Provides constants like `symprec`, `atomeps_structnew`.
    *   `tools_math`: For mathematical utilities (matrix operations, geometry).
    *   `tools_io`: For I/O operations and string/number conversions.
    *   `tools`: For various utility functions (e.g., `wscell`, sorting).
    *   `json_module`: For writing crystal information in JSON format.
    *   `graphics`: (in `@write`) For generating 3D models (`grhandle`).
    *   `iso_c_binding`: (in `@symmetry`, `@powderproc`) For C interoperability, likely with Spglib and NLopt.
    *   `hdf5`: (in `@vibrations`) For reading HDF5 files (conditional compilation).
    *   `nlopt.f`: (in `@powderproc`) For NLopt optimization library (conditional compilation).
*   **Interactions**:
    *   The `crystal` type is central and is passed to numerous procedures throughout the critic2 codebase.
    *   Many procedures within `crystalmod` and its submodules call each other to perform complex tasks. For instance, structure writers in `@write` rely on data accessors and transformation routines in `@proc`. Symmetry operations in `@symmetry` are used by editing functions in `@edit`.

```
