# DOCUMENTATION .//analysis.py

---

This python-script is the main workhorse to provide a comprehensive analysis tool for data from Geant4 or other simulation software as well as real data.

---

## Operation
Inputs:

- param_analysis.csv
  - this provides the program all its necessary parameters to function properly
  - parameters
    - *dummy* - a placeholder and should not be removed
    - *file_1...8* - the paths to the data, that should be used. Not needed line should not be removed, but just kept empty
    - *color scheme* - defines the color scheme for the spectra and should only include on of the following strings: c_rainbow, c_complementary, c_violetorange.
    - *plot title 1* - the title for the first spectrum

Outputs:

- procedure_<>.txt
  - this includes all the used parameters, mirroring param_analysis.csv, as well as all other used ones to have a complete record

---

## param_analysis.csv