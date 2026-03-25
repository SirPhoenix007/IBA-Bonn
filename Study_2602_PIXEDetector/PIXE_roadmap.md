# PIXE Roadmap

## Measurement Campaigns

#### Amersham 2084 Source

**Set C**: 600 seconds per element Amersham measurement

**Set D**: 3600 seconds per element Amersham measurement

**Set G**: Al+Si Wafer with Amersham

**Set H**: Coffee Bean with Amersham

**Set I**: Different Energy Gain values

#### Amptek Mini-X Source

**Set E**: Source further away

**Set F**: Using both anodes on *CalSam1*.

---


---

## Evaluations

### Energy Calibration

**Added EG values**: The value set for EG is manually adjustable, therefore more values in between can be realised to use with different expected values.

**Voigt profile**: As the line shape of a spectral line is Lorentzian in nature and the detector line shape is Gaussian, their convolution, namely a Cauchy-Voigt profile shall be tested for fitting the peaks.

**Background (reduction)**: *Workflow*: raw data -> background reduction -> translate into energies.

### Bremsstrahlung

*The Bremsstrahlung needs to be accounted for.*

### Peak Intensity Comparison

*There is data to compare against on whether the intensities are "properly" detected.*

### Multi-Voigt/Gauss summation for recreation of peak landscape?

**K-/L-lines**

Resulting plots:

- all ratios per element
- each ratio for all elements
  - compressed
  - with x-axis = Z for fit of paper

---



## Presenting Results

**Workflow Representation**

raw data -> (background reduction) -> peak detection -> energy calibration -> spectra

**Types of Spectra**

*A*: simple, calibrated

*B*: incl. background reduction

*C*: incl. Bremsstrahlung fit

*D*: incl. background reduction, incl. peaks

*E*: incl. Bremsstrahlung fit, incl. peaks

*F*: incl. missing peaks

*G*: x-ray spectrum of sample with highlighted, measured peaks


**Beam symboling script**

See notes for details on how.
