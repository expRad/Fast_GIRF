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

## Instructions for demo GSTF calculation

Do you want to process your own measurement data with our method? You can do so by following these steps:
* Download the demo dataset *triangle_data_demo.zip* from [zenodo](https://doi.org/10.5281/zenodo.7361610) and place the contents in the folder `/triangle_data_demo/ `. This is actually just a subset of */triangle_measurements/measurement_H_ref_x.mat*, containing only the data measured according to scheme (i) (cf. Figure 2(A) of our paper), ie 10 triangles played out during signal acquisition, similar to the approaches published by [Vannesjo et al.](https://doi.org/10.1002/mrm.24263) and [Campbell-Washburn et al.](https://doi.org/10.1002/mrm.25788)

* The data has the following structure (which you need to replicate with your measurement data):
  * The input data *input_H_demo.mat* contains three entries:
    * The array *grad_input* contains the input waveforms of the gradient pulses, here the 10 triangles, in mT/m on a time grid with 1us resolution. In our case, the input data starts at the start of every TR interval and ends at the end of the readout with the largest delay period within a TR interval. But in principle, the input waveforms during the readout are sufficient if only measurement scheme (i) is used.
    * The parameter *lengthADC* specifies the readout length in us.
    * The  array *shift* specifies the starting point of each readout within a TR interval in us. In the case of our 10 triangles measured with scheme (i), *shift* has 10 identical entries.

  * The output data *measurement_H_demo.mat* contains six entries:
    * The parameter *dwelltime* specifies the measurement dwell time in seconds.
    * The parameter *FOV* specifies the field of view in mm.
    * The parameter *numSlices* specifies the number of slices in which the signal was measured according to the thin-slice method (this has to be at least 2 slices for our code to work).
    * The parameter *orientation* specifies the slice orientation: 'dTra' for transversal slices, 'dSag' for sagittal slices, or 'dCor' for coronal slices. 
    * The array *PosSlices* contains the slice positions in mm in the order that corresponds to the raw data of the measurement.
    * The array *kspace* contains the raw measurement data. It is a 6-D array with the following dimensions: [data points per readout, coil elements, phase-encoding steps in direction 1, phase-encoding steps in direction 2, slices, 2 times number of triangles]. The two phase-encoding directions refer to the approach published by [Rahmer et al.](https://doi.org/10.1002/mrm.27902) If you used the thin-slice method without phase encoding, those array dimensions are just size 1. In the last dimensions, it is important that every second entry corresponds to the same triangle as the preceding entry, but measured with inverted sign.

* After you have saved your data in the described format, you can calculate your gradient system transfer function with the script `/GSTF_calculation/main_H_demo.m `. 
  * If you have named the *.mat* files with your input and output data differently, change lines 13 (`meas_name`) and 22 (`input_name`) accordingly. 
  * In line 46 (`prep.numTriang`), you can specify the number of triangle pulses you have measured (counting only the ones with positive sign, in our case 10).
  * If you want the result to be saved to `/GSTF_calculation/results/ `, set `doSaveGIRFs` in line 33 to 1 and, if needed, adjust `name2save` in line 35.
