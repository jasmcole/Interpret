Interpret
=========

**INTER**ferometry **P**hase **RET** rieval - perform phase retrieval, unwrapping and density retrieval of experimental interferograms. Requires Matlab 2012b or later. Tested on OSX.

![MainWindow](/Images/MainWindow.png)

## Contents
1. [Installation](#Installation)
2. [Set up a file naming convention](#set-up-a-file-naming-convention)
3. [Open an interferogram](#open-an-interferogram)
4. [Apply a calibration](#apply-a-calibration)
5. [Retrieve a phase profile](#retrieve-a-phase-profile)
6. [Unwrap a phase profile](#unwrap-a-phase-profile)
7. [Phase processing](#phase-processing)
8. [Density retrieval](#density-retrieval)
9. [Density analysis](#density-analysis)
10. [Batch process](#batch-process)

Installation
------------

1. Download the zipped Interpret repository and extract
2. Move `Interpret_placemeinroot.m` to the Matlab root directory, and rename it to 'Interpret.m'
3. Rename the Interpret-master folder to Interpret, and place it in the Matlab root directory.
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

See `ParseExperimentPath.m` for more details

Open an interferogram
---------------------

<img src="/Images/OpenFile.png" alt="LoadFile" style="width: 200px;"/>

If loading a file from a structured data directory, enter the correct numbers in the text boxes and click `Load File`.

If loading any other file, click `Other`.

Interpret uses the Matlab `imread` command for most images. For `.raw` images it will try a standard set of image sizes. For more details see `ReadRAW16bit.m`.

`Adjust contrast`

Opens an imcontrast widget to adjust the colourmap of the interferogram axes.

Apply a calibration
-------------------

<img src="/Images/Calibration.png" alt="Calibration" style="width: 200px;"/>

The drop-down calibration list is populated from the `CalibrationDatabase.csv` file, and loads the variables defined by their name and contents. These variables are displayed in the table, and by editing the `.csv` file Interpret will have access to different variables at runtime.

Select a calibration and click `Apply`. Any changes made are applied and saved automatically, unless the `Auto-apply updates` checkbox is unticked.

#### Create a new calibration

<img src="/Images/NewCalibration.png" alt="Calibration" style="width: 200px;"/>

Click the `New` button and type a name for the calibration. You can copy values from another calibration, or use default values.

Default calibration entries are
- `Name` - the name of the calibration
- `micperpix` - the number of microns per pixel in the interferograms
- `rotation` - the rotation applied to the image
- `x`/`y`/`w`/`h` - the ROI of the image. Clicking any of these entries will bring up an interactive rectangle over the uncropped image. Resize as necessary, then double-click inside the rectangle to set the ROI
- `xfft`/`yfft`/`wfft`/`hfft` - the ROI applied in the Fourier plane. See details below
- `reference` - the path to the reference interferogram. Selecting this entry will open a file selector dialog.
- `hsmooth1`/`hsmooth2` - smoothing parameters for the Hilbert transform
- `asmooth` - smoothing parameter for the inverse Abel transform
- `ymid` - the pixel coordinate of the symmetry axis of the interferogram
- `nfringes` - the number of fringes in the image, used in the wavelet transforms
- `CWT_fringelo`/`CWT_fringehi` - the range over which to fit fringe in the CWT methods. Set to 0.5/2 means scan over the range of 50% to 200% of the mean fringe wavelength.
- `CWT_ntheta` - the number of fringe angles to try in the CWT2D methods
- `CWT_nlambda` - the number of fringe wavelengths to try in the CWT methods
- `Moire_p_um` - the grating pitch if doing Moire deflectometry
- `Moire_d_mm` - the grating separation if doing Moire deflectometry


Retrieve a phase profile
------------------------

There are several methods to retrieve the phase shift.

#### FFT

- Fast, noise tolerant, can suffer from ringing.
- Works in the Fourier plane - a 2D FFT of the interferogram is displayed in the phase diagnostic axes.
- Click and drag a rectangle to choose a Fourier region to analyse. Typically one of the lobes either side of the origin.

<img src="/Images/FFT.png" alt="Calibration" style="width: 400px;"/>

- (optional) Save a Fourier region by populating the xfft, yfft, wfft and hfft fields in the calibration table.

#### Hilbert
- Fast, noise tolerant, sometimes doesn't completely remove fringes.
- Works row-by-row through the image, so the fringes need to be mostly **vertical**. If they are not, use the `rotation` attribute to rotate the image by 90 degrees.
- Needs to work on a sinusoidal signal oscillating about 0, therefore the mean of the image must be subtracted. The mean is a heavily smoothed copy of the image using smoothing parameter `hsmooth1`. This should be close to zero.
- After subtraction the data can be smoothed using `hsmooth2`. A value of 1 indicates no smoothing, values closer to 0 indicate more smoothing.
- The filtered interferogram is displayed in the phase diagnostic axes. If the fringes are broken, try changing `hsmooth1`.

#### CWT
- Slow, can be susceptible to noise, gives most accurate results.
- Again works row-by-row so the fringes should be close to vertical.
- Need to set the `CWT_nfringes` parameter to the number of fringes visible in the image.
- The phase diagnostic axes show the best CWT fit to the local fringe wavelength (blue dots). If the blue dots are all at the lower or upper limits of the plot, change `CWT_nfringes` and 'CWT_fringelo'/`CWT_fringehi` accordingly.
- If there is wide variation in fringe wavelength try increasing `CWT_nlambda`. Otherwise keep this at 2 to decrease errors in the fit.

#### CWT2D
- Fast, susceptible to noise, necessary when fringes are at a large angle.
- Works in any orientation.
- Try increasing `CWT_nlambda` and `CWT_ntheta` if the fit is poor.

Unwrap a phase profile
----------------------
After retrieval the phase will be wrapped within ± π. For large phase shifts it will need to be unwrapped. There are various techniques available, in practice some will work better than others on different images.

#### 1D unwrap
- Fast, works line-by-line
- Click the `Top`/`Bottom`/`Left`/`Right` buttons to perform simple unwraps from the edge of the image.

#### Volkov unwrap
- Fast, based on 2D FFT.
- Click the `Volkov unwrap` button. Sometimes trying multiple unwraps in a row can improve the unwrapping.

#### Goldstein unwrap
- Slow, robust.
- Click the `Goldstein unwrap` button to perform a slow Goldsten unwrap
- Click a pixel in the image where it is known that the phase isn't wrapped.

#### Costantini unwrap
- Slowish
- Uses a network flow method.

Phase processing
----------------

There are several functions available for modifying the retrieved phase map.

`Subtract gradient`

Subtract a smooth gradient from the phase profile. Useful if no reference exists, but not very accurate.

`Smooth phase`

Smooth the phase profile with a 2D cubic spline smoothing filter.

`Invert phase`

Negates the phase.

`Remove hot pixels`

Sometimes, especially after a Volkov unwrap, there will be small residual phase errors. Use this function to remove such errors.

`Region of interest`

Click and drag to crop the phase image. Useful if phase errors only occur at the edges

`Adjust contrast`

Bring up the `imcontrast` tool for the phase plot

Density retrieval
-----------------

Choose a laser wavelength and a refractive medium, and a density retrieval method. The micperpix calibration variable needs to be set for accurate density retrieval.

It is currently assumed that the symmetry axis is **vertical**. Therefore the phase map must be rotated if not, using the `Left`/`Right` buttons. If the symmetry axis is not perfectly vertical, the `rotation` property should be adjusted. Currently this requires a full phase retrieval again - should be updated in future.

If the `ymid` property is set to 0 the user is prompted to select a symmetry axis, otherwise the axis is set to the pixel coordinate in `ymid`.

#### Abel inversion

- The Abel inversion is performed line-by-line, smoothed according to the `asmooth` calibration field. `asmooth = 1` corresponds to no smoothing, and `asmooth = 0` to maximum smoothing. This variable should be used in a logarithmic sense, i.e. one should try values of 0.1, 0.01, 0.001 etc.

#### Trapezoidal fit

- Attempts to fit a trapezoidal density profile to the phase map. Currently untested, should be avoided.

Density analysis
----------------

Once a density map has been calculated it is displayed in the density axes, with a lineout down the symmetry axis plotted in the density diagnostic axes.

`Density lineout`

Select a point in the density image to plot orthogonal lineouts through that point.

`Get average`

Click, from left-to-right, to points in the lineout to define a range over which to analyse. The mean density and standard deviation are indicated on the plot and in the status box.

<img src="/Images/MeanDensity.png" alt="MeanDensity" style="width: 400px;"/>

`Save to workspace/file`

Save the phase and density profiles to the Matab workspace in the form of a struct called `InterpretData`, or to a `.mat` file.

Batch process
-------------

If reading data from structured directories, it is possible to batch process a number of shots. The two text boxes mark the region of shots to process (e.g., from shot 2 to 15).

By default, for every image loaded, the currently applied calibration is used, the currently selected phase retrieval method is used, and the currently selected density retrieval method is used.

If additional commands are required, e.g. unwrapping, click the `Edit Batch Commands` button. This opens the file `EditableBatchCommands.m` which returns a cell array of strings. Each string corresponds to a button click, and clicks are executed in order.

After completion, a struct `InterpretDataBatch` is saved to the Matlab workspace.
