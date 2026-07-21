<img width="500" alt="image" src="https://github.com/myroms/roms/assets/23062912/2fa815f4-df51-4671-b9ed-188995b87a7b">

# Auxiliary ROMS Programs:

- **`types.F`**: It checks the precision and range of floating-point
  variables in a particular compiler architecture. To compile and
  link, use:

  ```
  ifort -o types.x types.F
  gfortran -o types.x types.F
  ```

- **`yaml_parser_test.F`**: It tests the **yaml_parser.F** available in
  **ROMS**. It requires the Fortran 2003 standard and a few features of
  the Fortran 2008 release, which are available in modern compilers.
  We have reports from users having issues compiling with older versions of
  **gfortran**. To compile and link, use:

  ```
  ifort -o yaml_parser_test.x yaml_parser_test.F
  gfortran -o yaml_parser_test.x yaml_parser_test.F
  ```
  It is easy to run by specifying the desired **YAML** files available in **ROMS**:

  ```
  > yaml_parser_test.x

  Enter YAML filename:  ../../ESM/coupling_esmf.yaml
   ...
  Enter YAML filename:  ../../ESM/roms_cmeps.yaml
   ...
  Enter YAML filename:  ../External/varinfo.yaml

  Enter YAML filename:
  ```
