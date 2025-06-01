# `fieldmod.f90`

## Overview

The `fieldmod.f90` module is designed to handle scalar fields within the critic2 application. It defines the `field` derived type, which encapsulates the data and methods necessary to evaluate a scalar field (such as electron density, electrostatic potential, or custom fields) and its derivatives at any point in space. This module supports various sources for field data, including promolecular densities, grid-based data from other quantum chemistry codes (like WIEN2k, ELK, VASP, Gaussian Cube files), and custom mathematical expressions (ghost fields).

A key functionality of this module is the analysis of the topology of scalar fields, including the location and characterization of critical points (CPs) and the tracing of gradient paths.

The primary submodule associated with `fieldmod.f90` is:
*   `fieldmod@proc.f90`: Implements the core logic for field evaluation, numerical differentiation, critical point searching algorithms (Newton-Raphson), and gradient path integration.

## Key Components

### Derived Types/Classes

*   **`field`**: Represents a scalar field.
    *   **Key Fields/Attributes**:
        *   `c` (`type(crystal), pointer`): Pointer to the associated crystal structure.
        *   `id` (`integer`): Unique identifier for the field.
        *   `isinit` (`logical`): True if the field is initialized.
        *   `type` (`integer`): Type of the field (e.g., `type_promol`, `type_grid`, `type_wfn`, `type_ghost`). See "Important Variables/Constants" for specific types.
        *   `usecore` (`logical`): If true, augment the field with core densities.
        *   `numerical` (`logical`): If true, use numerical differentiation.
        *   `exact` (`logical`): If true, use exact (analytical if available) calculations, otherwise approximate.
        *   `typnuc` (`integer`): Signature of nuclear critical points (-3 NATT, -1 NNIM, +1 NBIM, +3 NCPS).
        *   `name` (`character(len=mlen)`): Name of the field.
        *   `file` (`character(len=mlen)`): Source file for the field data.
        *   `elk` (`type(elkwfn), allocatable`): Data for ELK densities.
        *   `wien` (`type(wienwfn), allocatable`): Data for WIEN2k densities.
        *   `pi` (`type(piwfn), allocatable`): Data for PI wavefunctions.
        *   `grid` (`type(grid3), allocatable`): Data for generic grid-based fields.
        *   `wfn` (`type(molwfn), allocatable`): Data for GTO/STO atom-centered wavefunctions.
        *   `dftb` (`type(dftbwfn), allocatable`): Data for DFTB wavefunctions.
        *   `fr` (`type(fragment), allocatable`): Fragment for fragment-based promolecular density.
        *   `zpsp` (`integer, allocatable(:)`): Pseudopotential charges for core density augmentation.
        *   `expr` (`character(len=mmlen)`): Mathematical expression for ghost fields.
        *   `sptr` (`type(c_ptr)`): Pointer to the parent system (for ghost fields).
        *   `fcp_deferred` (`logical`): True if calculation of CPs on nuclei was deferred.
        *   `ncp` (`integer`): Number of non-equivalent critical points.
        *   `cp` (`type(cp_type), allocatable(:)`): Array of non-equivalent critical points.
        *   `ncpcel` (`integer`): Number of critical points in the complete list (unit cell).
        *   `cpcel` (`type(cp_type), allocatable(:)`): Array of all critical points in the unit cell.
    *   **Key Bound Procedures (many implemented in `@proc`)**:
        *   `end`: Deallocates data and uninitializes the field.
        *   `set_default_options`: Sets default options for the field.
        *   `set_options`: Sets field options from a command string.
        *   `field_new`: Creates a new field from a `fieldseed`.
        *   `load_promolecular`: Loads a promolecular density field.
        *   `load_as_fftgrid`: Loads a field by transforming a 3D grid (e.g., for derivatives).
        *   `load_ghost`: Loads a field defined by a mathematical expression.
        *   `grd`: Calculates the field value and its derivatives (up to 2nd order) at a given point.
        *   `grd0`: Calculates only the field value at a given point.
        *   `der1i`, `der2ii`, `der2ij`: Perform numerical differentiation.
        *   `typestring`: Returns a string identifying the field type.
        *   `printinfo`: Prints information about the field.
        *   `write_json`: Writes field information in JSON format.
        *   `init_cplist`: Initializes the critical point list, typically with atomic positions.
        *   `init_cplist_deferred`: Calculates field properties at nuclear positions (deferred from `init_cplist`).
        *   `nearest_cp`: Finds the nearest critical point to a given position.
        *   `identify_cp`: Identifies a critical point by its position.
        *   `testrmt`: Tests for discontinuities at muffin-tin radii (for WIEN2k/ELK fields).
        *   `benchmark`: Tests the speed of field evaluation.
        *   `newton`: Performs a Newton-Raphson search for a critical point.
        *   `addcp`: Adds a new critical point to the list after validation.
        *   `sortcps`: Sorts the critical point list.
        *   `gradient`: Traces a gradient path from a given point.

### Module Procedures

*   `realloc_field`: Utility subroutine to reallocate an array of `field` type.

## Important Variables/Constants

*   **`type_uninit`, `type_promol`, `type_grid`, `type_wien`, `type_elk`, `type_pi`, `type_wfn`, `type_dftb`, `type_promol_frag`, `type_ghost`** (`integer, parameter`): Public constants defining the different types of scalar fields supported.
*   **`ndif_jmax`** (`integer, parameter`): Likely related to the maximum order or number of points for numerical differentiation stencils (value: 10).
*   **`flooreps`** (`real*8, parameter, private` in `@proc`): Small epsilon value (1d-4) used for boundary checks around the unit cell.
*   **`derw`, `derw2`, `big`, `safe`** (`real*8, parameter, private` in `@proc`): Parameters for Richardson's extrapolation in numerical differentiation.
*   **`hini`, `errcnv`** (`real*8, parameter, private` in `@proc`): Initial step size and convergence error for numerical differentiation.

## Usage Examples

```fortran
use fieldmod
use crystalmod, only: crystal
use fieldseedmod, only: fieldseed
use types, only: scalar_value
implicit none

type(crystal) :: crys
type(field) :: fld
type(fieldseed) :: fs
type(scalar_value) :: sv
real*8 :: point_cart(3)

! ... (Initialize 'crys' and 'fs' appropriately) ...
! For example, load a crystal structure into 'crys'
! and set 'fs' to load a promolecular field:
! fs%iff = ifformat_promolecular
! fs%fid = "my_promol_field"

! Create and initialize the field
call fld%field_new(fs, crys, id=1, sptr=c_null_ptr, errmsg=err_msg)
if (len_trim(err_msg) > 0) then
    print *, "Error initializing field: ", trim(err_msg)
    stop
end if

! Define a point in Cartesian coordinates (Bohr)
point_cart = (/ 0.1d0, 0.2d0, 0.3d0 /)

! Evaluate the field and its first two derivatives at the point
call fld%grd(point_cart, nder=2, res=sv)

! Print results
print *, "Field value at point: ", sv%f
if (sv%avail_der1) then
    print *, "Gradient at point: ", sv%gf
endif
if (sv%avail_der2) then
    print *, "Laplacian at point: ", sv%del2f
    print *, "Hessian matrix:"
    print '(3(F12.6,X))', sv%hf(1,:)
    print '(3(F12.6,X))', sv%hf(2,:)
    print '(3(F12.6,X))', sv%hf(3,:)
endif

! Clean up
call fld%end()
! ... (call crys%end()) ...
```

## Dependencies and Interactions

*   **Uses Modules**:
    *   `fieldseedmod`: To get initial field parameters via `fieldseed` type.
    *   `crystalmod`: Heavily relies on the `crystal` type for structural information.
    *   `fragmentmod`: For `type_promol_frag` using the `fragment` type.
    *   `elk_private`, `wien_private`, `pi_private`, `grid3mod`, `wfn_private`, `dftb_private`: These modules provide the underlying data structures and routines for specific field types (e.g., `elkwfn`, `grid3`, `molwfn`).
    *   `param`: Provides constants like `mlen`, `mmlen`, `maxzat0`.
    *   `types`: Provides `cp_type` for critical points, `scalar_value` for field evaluation results, `gpathp` for gradient paths, `thread_info`.
    *   `hashmod`: For hashing functionality (not explicitly shown in provided snippets but likely used internally).
    *   `iso_c_binding`: For C interoperability, particularly for `c_ptr` used with ghost fields.
    *   `arithmetic` (in `@proc`): For evaluating mathematical expressions in ghost fields.
    *   `global` (in `@proc`): For global parameters like `CP_hdegen`, `nav_step`.
    *   `tools_math` (in `@proc`): For mathematical utilities like `rsindex`.
    *   `tools_io` (in `@proc`): For I/O and string/number conversions.
*   **Interactions**:
    *   The `field` type is fundamental for any analysis involving scalar fields in critic2.
    *   The `grd` and `grd0` methods are the primary interfaces for evaluating the field.
    *   Critical point search and gradient path tracing rely on repeated calls to `grd`.
    *   The module interacts with various private modules (`elk_private`, `wien_private`, etc.) that handle the specifics of reading and interpreting data from different quantum chemistry codes.

```
