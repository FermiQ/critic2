# `arithmetic.F90`

## Overview

The `arithmetic.F90` module provides a powerful expression evaluator capable of parsing and computing mathematical and physical property expressions. It allows for the definition of variables, the use of various mathematical functions, and the evaluation of scalar fields (and their derivatives) defined elsewhere in the critic2 system. This module is designed to be thread-safe.

The core functionality is implemented in the `arithmetic@proc.F90` submodule, which handles tokenization of expressions, parsing using a shunting-yard-like algorithm, and evaluation of the resulting postfix expression. It supports a wide range of unary and binary operators, mathematical functions (e.g., `sin`, `cos`, `exp`, `log`), and special chemical functions (e.g., Thomas-Fermi kinetic energy `gtf`, ELF `elf`, LOL `lol`). It can also access structural variables related to the system's geometry (e.g., distance to nearest nucleus `dnuc`, atomic coordinates `xc`, `yc`, `zc`).

## Key Components

### Module Procedures

*   **`eval(expr, errmsg, x0, sptr, periodic)`**: `recursive module function real*8`
    *   **Overview**: Evaluates a given arithmetic expression string.
    *   **Arguments**:
        *   `expr` (`character(*), intent(in)`): The arithmetic expression to evaluate.
        *   `errmsg` (`character(len=:), allocatable, intent(inout)`): Error message string; populated if evaluation fails.
        *   `x0` (`real*8, intent(in), optional`): 3D Cartesian coordinates (Bohr) at which to evaluate fields or structural variables referenced in `expr`.
        *   `sptr` (`type(c_ptr), intent(in), optional`): C pointer to the parent `system` object, necessary if `expr` involves fields or structural variables.
        *   `periodic` (`logical, intent(in), optional`): If true (default), evaluates fields considering periodic boundary conditions. If false, treats the system as non-periodic.
    *   **Returns**: (`real*8`) The result of the expression evaluation.

*   **`eval_grid(n, expr, sptr, f, iok)`**: `module subroutine`
    *   **Overview**: Evaluates an arithmetic expression over a 3D grid.
    *   **Arguments**:
        *   `n` (`integer, intent(in)(3)`): Dimensions of the grid.
        *   `expr` (`character(*), intent(in)`): The arithmetic expression to evaluate.
        *   `sptr` (`type(c_ptr), intent(in)`): C pointer to the parent `system` object.
        *   `f` (`real*8, intent(out)(n(1),n(2),n(3))`): Output grid containing the evaluated expression.
        *   `iok` (`logical, intent(out)`): True if evaluation was successful, false otherwise.

*   **`fields_in_eval(expr, errmsg, n, idlist, sptr)`**: `module subroutine`
    *   **Overview**: Parses an expression and returns a list of field identifiers (names or IDs) found within it.
    *   **Arguments**:
        *   `expr` (`character(*), intent(in)`): The expression to parse.
        *   `errmsg` (`character(len=:), allocatable, intent(inout)`): Error message string.
        *   `n` (`integer, intent(out)`): Number of field identifiers found.
        *   `idlist` (`character(len=mlen), allocatable, intent(inout)(:)`): Array populated with the field identifiers.
        *   `sptr` (`type(c_ptr), intent(in)`): C pointer to the parent `system` object.

*   **`setvariable(ikey, ival)`**: `module subroutine`
    *   **Overview**: Defines or updates a user variable that can be used in expressions.
    *   **Arguments**:
        *   `ikey` (`character(*), intent(in)`): Name of the variable.
        *   `ival` (`real*8, intent(in)`): Value to assign to the variable.

*   **`isvariable(ikey, ival)`**: `module function logical`
    *   **Overview**: Checks if a variable with the given key is defined.
    *   **Arguments**:
        *   `ikey` (`character(*), intent(in)`): Name of the variable.
        *   `ival` (`real*8, intent(out)`): If the variable exists, its value is returned here.
    *   **Returns**: (`logical`) True if the variable is defined, false otherwise.

*   **`clearvariable(ikey)`**: `module subroutine`
    *   **Overview**: Deletes a user-defined variable.
    *   **Arguments**:
        *   `ikey` (`character(*), intent(in)`): Name of the variable to delete.

*   **`clearallvariables()`**: `module subroutine`
    *   **Overview**: Deletes all user-defined variables.

*   **`listvariables()`**: `module subroutine`
    *   **Overview**: Prints a list of all currently defined user variables and their values.

*   **`listlibxc(doref, doname, doflags)`**: `module subroutine`
    *   **Overview**: Lists available Libxc functionals if critic2 is compiled with Libxc support.
    *   **Arguments**:
        *   `doref` (`logical, intent(in)`): If true, include references for each functional.
        *   `doname` (`logical, intent(in)`): If true, include the long name of each functional.
        *   `doflags` (`logical, intent(in)`): If true, include flags describing the properties of each functional.

## Important Variables/Constants (Private in `@proc`)

The `@proc` submodule defines numerous private integer parameters representing different functions and operators used during parsing and evaluation (e.g., `fun_plus`, `fun_sin`, `fun_elf`, `svar_dnuc`). These are used internally by the tokenizer and evaluator.

*   **Mathematical Operators**: `fun_uplus`, `fun_uminus`, `fun_power`, `fun_plus`, `fun_minus`, `fun_prod`, `fun_div`, `fun_modulo`.
*   **Logical Operators**: `fun_leq`, `fun_geq`, `fun_equal`, `fun_neq`, `fun_and`, `fun_or`, `fun_great`, `fun_less`.
*   **Mathematical Functions**: `fun_abs`, `fun_exp`, `fun_sqrt`, `fun_floor`, `fun_ceiling`, `fun_round`, `fun_log`, `fun_log10`, `fun_sin`, `fun_asin`, `fun_cos`, `fun_acos`, `fun_tan`, `fun_atan`, `fun_atan2`, `fun_sinh`, `fun_cosh`, `fun_erf`, `fun_erfc`, `fun_min`, `fun_max`.
*   **Chemical/Physical Functions**: `fun_xc` (Libxc functional), `fun_gtf` (Thomas-Fermi KE), `fun_vtf` (TF Virial Potential Energy), `fun_htf` (TF Total Energy), `fun_gtf_kir` (TF KE with Kirzhnits correction), `fun_elf` (Electron Localization Function), `fun_lol` (Localized-Orbital Locator), `fun_mep` (Molecular Electrostatic Potential), `fun_rdg` (Reduced Density Gradient), and others.
*   **Structural Variables**: `svar_dnuc` (distance to nucleus), `svar_xnucx` (x-coord of nearest nucleus, cryst.), `svar_x` (current x-coord, default units), `svar_xc` (current x-coord, Cartesian), etc.
*   **Token Types**: `token_num`, `token_fun`, `token_op`, `token_lpar`, `token_rpar`, `token_comma`, `token_field`, `token_structvar`.

## Usage Examples

```fortran
use arithmetic
use systemmod, only: system ! Assuming 'syl' is an initialized system object
use iso_c_binding, only: c_loc
real*8 :: result, x_coord(3)
character(len=200) :: expr_str, err_msg

! Define a variable
call setvariable("my_var", 10.5d0)

! Simple arithmetic
expr_str = "2 * (my_var + sin(0.5))"
result = eval(expr_str, err_msg)
if (len_trim(err_msg) == 0) then
    print *, "Result of '", trim(expr_str), "': ", result
else
    print *, "Error: ", trim(err_msg)
end if

! Expression involving a field (e.g., electron density '$rho')
! Assume 'syl' is an initialized system object pointer and x_coord is set
x_coord = (/ 0.0d0, 0.0d0, 0.0d0 /) ! Cartesian coordinates in Bohr
expr_str = "$rho + 0.1 * $rho:x" ! Density + 0.1 * x-derivative of density
result = eval(expr_str, err_msg, x0=x_coord, sptr=c_loc(syl))
if (len_trim(err_msg) == 0) then
    print *, "Result of '", trim(expr_str), "' at (0,0,0): ", result
else
    print *, "Error: ", trim(err_msg)
end if

! Expression involving a structural variable
expr_str = "@dnuc * 2.0" ! Distance to nearest nucleus times two
result = eval(expr_str, err_msg, x0=x_coord, sptr=c_loc(syl))
if (len_trim(err_msg) == 0) then
    print *, "Result of '", trim(expr_str), "' at (0,0,0): ", result
else
    print *, "Error: ", trim(err_msg)
end if
```

## Dependencies and Interactions

*   **Uses Modules**:
    *   `iso_c_binding`: For `c_ptr` when interacting with the `system` object.
    *   `systemmod` (via `sptr`): To access field data (`system%f`) and crystal structure information (`system%c`) when evaluating expressions containing field references (`$field_id`) or structural variables (`@var_name`).
    *   `hashmod` (via `param%vh` in `@proc`): Used internally for storing and retrieving user-defined variables.
    *   `param` (in `@proc`): Accesses the variable hash `vh`.
    *   `tools_io` (in `@proc`): For string manipulation and number/type checking.
    *   `fieldmod` (in `@proc`): To call `f%grd()` for field evaluations.
    *   `xc_f90_lib_m` (conditionally in `@proc`): If compiled with Libxc, this module is used to evaluate exchange-correlation functionals.
*   **Interactions**:
    *   The `eval` function is the primary public interface.
    *   It relies heavily on the `system` module (passed via `sptr`) to resolve field identifiers and access their `grd` methods for evaluation at specific points.
    *   User-defined variables set by `setvariable` are stored globally (likely via `param%vh`) and can be accessed in subsequent `eval` calls.

```
