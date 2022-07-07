# NewOASIM library

NewOASIM library is a re-implementation of OASIM software developed by NASA. It comes with two shared objects: `liboasim.so`, a Fortran 2018 implementation and a Python wrapper `liboasim-py.so` which can be used with `ctypes` module in Python.

## Main Features

The library is conceived to be used _inside_ a parallel environment (aka MPI or OpenMP) and therefore it uses a single thread. 

- **oasim_lib**: is an object that read the configuration file (see below) and store the latitude and longitude mesh of the whole system.

- **calc_unit**: is a calculator which depends on the oasim_lib object it refers to, and process the input. Every calc_unit is independent from another. This allows the user to fix an environment using a single oasim_lim object and to have several calc_units performing different computations in parallel. 

## Requirements

- `gcc` and/or `intel` Fortran 2018 capable compiler;

- `cmake` version 3.0 or later;

- `makefile`.

## Compiling NewOASIM

- run one of the `build_<rel>_<comp>` script (where `<rel>` is either `debug` or `release` and `<comp>` is either `gcc` or `intel`) 

```bash
./build_<rel>_<comp>
```

- enter the corresponding folder

```bash
cd builds/<rel>_<comp>
```

- launch the Makefile

```bash
make
```

The results of the compilation are in the folder `builds/<rel>_<comp>/OASIMlib` in which there are

- `liboasim.so`;

- `liboasim-py.so`.

## Minimal Working Examples

Two minimal working examples are provided:

- `builds/<rel>_<comp>/test/f90/test.x` which tests the functionality of `liboasim.so`;

- `builds/<rel>_<comp>/test/python/oasim.py` which tests the functionality of `liboasim-py.so`.

## How to use OASIMlib

### Configuration file

As mentioned before, one of the main ingredient for OASIMlib is the configuration file. It it a `.yaml` file with the following fileds:

```yaml
input:                                 # input section
    atmospheric_data_file:             # string: specify the path of the atmosperich file (see below)
    absorption_data_file:              # string: specify the path of the absorption file (see below)
    cloud_slingo_file:                 # string: specify the path ot the cloud data file (see below)
compute:                               # compute section
    integration_step_secs:             # real: specify the tentative step of integration (see below)
    max_integration_steps:             # integer: specify the maximum number of steps of integrations (see below)
    zenit_avg:                         # boolean: specify whether to use average zenith angle instead of integration
    local_time:                        # boolean: specify whether to use local time
    am:                                # am paramer of the model
    vi:                                # vi parameter of the model
output:                                # output section
    bin_file:                          # string: specify the path of file containing the output binning
```

#### Integration

The of the sine of the solar zenith angle is performed if the `zenith_avg` parameter is `false`. Otherwise the sine of the zenith angle at the beginning of the interval is used. In the case of integration one can specify two parameters: `integration_step_secs` specify the width of the integration interval in seconds; the second parameter `max_integration_steps` is used to limit the maximum number of integration step. Therefore the actual number of integration steps is the minimum between the floor of the one obtained dividing the interval by `integration_step_secs` and `max_integration_steps`.

#### Data file format

In all these files lines starting with `#` are not parsed and considered as comment

- `atmospheric_data_file`: 2 columns for range of wavelength and 6 columns of data;

- `absorption_data_file`: 2 columns for range of wavelength and 2 columns of data;

- `cloud_slingo_file`: 2 columns for range of wavelength and 6 columns of data;

- `bin_file`: 2 columns for range of wavelength.

Notice that for the internal computation the program uses the binning specified in the `bin_file` path; this means that  all the other binning are converted at the beginning to the output specification.

### Fortran

As shown in the f90 test, it is necessary to import the module

```fortran
use oasim
```

The main oasim_lib object can be initialized in the following way

```fortran
type(oasim_lib) :: lib
...

lib = oasim_lib(config_path, lat, lon, error)
```

where

- `config_path` (`character (len=*)`): specify the path of the `.yaml` configuration file;
- `lat` (`real (kind=real_kind), dimension(:)`) is the *1d* mesh specifying the latitude mesh;
- `lon` (`real (kind=real_kind), dimension(:)`) is the *1d* mesh specifying the longitude mesh;
- `error` (`logical`): is the output to signal errors in the initialization process.

The calc_unit object can be initialized in the following way

```fortran
type (calc_unit) :: calc
...

calc = calc_unit(p_size, lib)
```

where

- `p_size` (`integer`): specify the size of the batch of points to be processed;
- `lib` (`type(oasim_lib), target`): specify the pointer to the oasim_lib object it depends on.

The calculation method (subroutine) is called `monrad`:

```fortran
calc%monrad(points, year, day, sec_b, sec_e, slp, wsm, oz, wv, rh, ccov, cdre, taua, asymp, ssalb, edout, esout, error)
```

the input are the following (here `bin_n` is the number of the bin specified in the `bin_file`):

- `points` (`integer, dimension(:)`): specify the indices of the points on the mesh to work on (the number of point must coincide with `size_p`);
- `year` (`integer`): specify the year of the simulation;
- `day` (`integer`): specify the day of the year of the simulation;
- `sec_b` (`real (kind=real_kind)`): specify the beginning second of the day of the simulation;
- `sec_e` (`real (kind=real_kind)`): specify the end second of the day of the simulation;
- `sp` (`real (kind=real_kind), dimension(:)`): surface pressure [Pa], must have `size_p` length;
- `msl` (`real (kind=real_kind), dimension(:)`): air pressure at mean sea level [Pa], must have `size_p` length;
- `ws10` (`real (kind=real_kind), dimension(:)`): wind Speed [m/s], must have `size_p` length;
- `tc03` (`real (kind=real_kind), dimension(:)`): total column ozone [kg m^-2] , must have `size_p` length;
- `t2m` (`real (kind=real_kind), dimension(:)`): 2 metre temperature [K], must have `size_p` length;
- `d2m` (`real (kind=real_kind), dimension(:)`): 2 metre dewpoint temperature [K], must have `size_p` length;
- `tcc` (`real (kind=real_kind), dimension(:)`): total cloud cover [0-100], must have `size_p` length;
- `tclw` (`real (kind=real_kind), dimension(:)`): total column cloud liquid water [kg m^-2], must have `size_p` length;
- `cdrem` (`real (kind=real_kind), dimension(:)`): cloud droplet effective radius [um], must have `size_p` length;
- `taua` (`real (kind=real_kind), dimension(:)`): aerosol optical thickness, must have `size_p` x `bin_n` length;
- `asymp` (`real (kind=real_kind), dimension(:)`): aerosol asymmetry parameter, must have `size_p` x `bin_n` length;
- `ssalb` (`real (kind=real_kind), dimension(:)`): aerosol single scattering albedo [-], must have `size_p` x `bin_n` length;

the output are the following

- `edout` (`real (kind=real_kind), dimension(:)`): contains ..., must have `size_p` length;
- `esout` (`real (kind=real_kind), dimension(:)`): contains ..., must have `size_p` length;
- `error` (`logical`): it signals error in the computation procedure.

### Python

In the Python wrapper the  methods are almost one-to-one copy of the Fortran one, therefore we will explain just the arguments that differ from the Fortran implementation. Here all the array are from `numpy`.

```python
class oasim_lib:
    def __init__(self, lib_path, config_path, lat, lon):
        ...
```

here `lib_path` is the path of the the `lib-oasim-py.so` library.

```python
class calc_unit
    def __init__(self, p_size, lib):
        ...
    def monrad(self, points, year, day, sec_b, sec_e, slp, wsm, oz, wv, rh, ccov, cdre, taua, asymp, ssalb):
        ...
```

Here the input is identical while output of the method corresponds to `edout, esout`.
