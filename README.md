Interpret
=========

INTERferometry Phase RETrieval - perform phase retrieval, unwrapping and density retrieval of experimental interferograms. Requires Matlab 2012b or later.

To set up
---------

1. Add `Interpret.m` to the Matlab root directory, containing the code:

`mroot = fileparts(which('Interpret.m'));`
`addpath([mroot '/Interpret/InterpretSource'])`
`addpath([mroot '/Interpret'])`
`cd([mroot '/Interpret'])`
`IntGui`

2. Create a folder in the Matlab root directory called `Interpret`
3. Move all of the files to the `Interpret` folder
4. Run by typing `Interpret` in the command window

Set up a file naming convention
-------------------------------

If you have a structured series of files you would like to access using a known naming convention, you can edit `ExperimentDatabase.csv`. The first column is the identifier, typically the experiment name. The second column is the file path format. The file path format is structured as follows, where the length of each string is the amount of zero padding applied to the number:

- YY corresponds to the year text box
- MM corresponds to the month text box
- DD corresponds to the day text box
- RR corresponds to the run text box
- SS corresponds to the shot text box
- * corresponds to a wildcard
- ^ corresponds to the shot text box with no zero padding

See `ParseExperimentPath` for more details

Open any file
-------------

Click the `Other` button and select a file

Load a calibration
------------------

The drop-down list is populated from the `CalibrationDatabase.csv` file, and loads the variables defined by their name and contents. These variables are displayed in the table, and by editing the `.csv` file Interpret will have access to different variables at runtime.

Select a calibration and click `Apply`. Any changes made are applied and saved automatically, unless the `Auto-apply updates` checkbox is unticked.

Create a new calibration
------------------------

Click the `New` button and type a name for the calibration. This calibration is selected and populated with default values, then saved automatically.

Set the reference field to `None` if no reference exists.


Retrieve a phase profile
------------------------

FFT:

- Click and drag to choose a Fourier region to analyse.
- Save a Fourier region by populating the xfft, yfft, wfft and hfft fields in the calibration table

Unwrap a phase profile
----------------------

1D unwrap:

- Click the `Top`/`Bottom`/`Left`/`Right` buttons to perform simple unwraps from the edge of the image

Volkov unwrap

- Click the `Volkov unwrap` button to perform a fast 2D unwrap. Sometimes trying multiple unwraps in a row can improve the unwrapping

Goldstein unwrap

- Click the `Goldstein unwrap` button to perform a slow Goldsten unwrap
- Click a pixel in the image where it is known that no phase ambiguities occur
- Wait...

Phase analysis buttons
----------------------

`Subtract gradient`

- Subtract a smooth gradient from the phase profile. Useful if no reference exists, but not very accurate

`Invert phase`

- Negates the phase profile

`Region of interest`

- Click and drag to crop the phase image. Useful if phase errors only occur at the edges

`Smooth phase`

- Smooth the phase profile with a 2D cubic spline smoothing filter

`Adjust contrast`

- Bring up the `imcontrast` tool for the phase plot

Density retrieval
-----------------

Choose a laser wavelength and a refractive medium, and a density retrieval method. The micperpix calibration variable needs to be set for accurate density retrieval.

Abel inversion

- Upon clicking, the user is prompted to manually select a symmetry axis. This can also be set using the `ymid` calibration field.
- The Abel inversion is performed line-by-line, smoothed according to the `asmooth` calibration field. `asmooth = 1` corresponds to no smoothing, and `asmooth = 0` to maximum smoothing. This variable should be used in a logarithmic sense, i.e. one should try values of 0.1, 0.01, 0.001 etc.
