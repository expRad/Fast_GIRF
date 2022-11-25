# Fast measurement of the gradient system transfer function (GSTF)

This repository contains instructions and source code to reproduce the results presented in:

> Fast measurement of the gradient system transfer function at 7 T
> Scholten H, Lohr D, Wech T, Köstler H.
> Magnetic Resonance in Medicine. 2022. *Early view*. [DOI](https://doi.org/10.1002/mrm.29523)

Please cite this work if you use the content of this repository in your project.

The code is written in MATLAB (R2019b). In order for all scripts to run faultless, all folders in this repository must be added to the MATLAB path.

### Download instructions for measurement data:

*The measurement data of the triangular and trapezoidal gradients is avaliable [here](https://doi.org/10.5281/zenodo.7361610). After you have downloaded and extracted the .zip archives from zenodo, move the .mat files contained in triangle_measurements.zip into the folder* triangle_measurements, *and the .mat files from trapezoid_measurements.zip into the folder* trapezoid_measurements. *(Or replace the respective folders from this repository with the extracted ones.)*

## Determination of the GSTF

Run the script `/GSTF_calculation/main_H_fast.m ` to calculate the transfer function H<sub>fast</sub> as described in our paper. In the second section of the code ( `%% Select data files` ), you can specify the axis for which you wish to calculate the GSTF by commenting in or out the respective lines defining `meas_name` and `ax_name`.

In order to calculate the transfer function H<sub>ref</sub>, use the script `/GSTF_calculation/main_H_ref.m`.

The required gradient data (ie, waveforms of the input triangles and the respective measurement raw data) are to be contained in `/triangle_measurements/`. The results are automatically saved to ` /GSTF_calculation/results/`. With the script ` /GSTF_calculation/results/plotGSTFs.m`, you can re-create Figure 4 and Figure 5 from our paper.

## GSTF-based gradient predictions

For the comparison of the measured time-courses of a trapezoidal test gradient with the GSTF-based predictions, run the script ` /gradient_calculation/GradMeas_vs_GSTFpredict.m`. Here you can also specify the gradient axis in the second section of the code ( `%% Set file paths` ) by setting the variable `ax` to the respective value. The script reproduces Figure 6 and Figure 7 (or Figures S1 and S2, or Figures S3 and S4, depending on the selected axis) of our paper.

## Gradient measurements with GSTF-based pre-emphasis

To calculate the pre-emphasized input waveforms of the trapezoidal test gradient as described in our paper, run the script ` /pre_emphasis/pre_emphasis.m`. Here, you also need to specify the gradient axis in the second section of the code. The script creates Figure 8 of the paper for the specified gradient axis, and Figure 9 (or S5, or S6, respectively). The raw data of the corresponding gradient measurements are to be located in `/trapezoid_measurements/`.
