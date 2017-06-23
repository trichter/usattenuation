## USAttenuation: Intrinsic and scattering attenuation for the contiguous U.S.

This repository contains the processing, results and visualization for the following publication:

Tom Eulenfeld and Ulrich Wegler (2017), Crustal intrinsic and scattering attenuation of high-frequency shear waves in the contiguous United States, *Journal of Geophysical Research: Solid Earth*, 122, doi:[10.1002/2017JB014038](https://dx.doi.org/10.1002/2017JB014038). [[pdf](http://www.geophysik.uni-jena.de/igwphymedia/_users/eule/Eulenfeld_Wegler_2017_US_intrinsic_and_scattering_attenuation_Draft.pdf)]

View results comfortably here: [*Qviewer*](https://trichter.github.io/usattenuation).
To run the *Qviewer* application offline, clone or download the repository and open the file `visualization/Qviewer/index.html` in a browser. It will be a much smoother experience.

Content of folders||
--- | ---
`processing`     | Python scripts to reproduce figures of the publication and everything in the results folder
`results`        | Results in JSON format, figures of the paper and electronic supplement, georeferenced TIF files
`visualization`  | Javascript app to explore results, example Jupyter notebook on how to use JSON and georeferenced TIF files
`data`           | Downloaded data (waveforms, event and station metadata)
`tmp`            | Temporary files
**Processing**   |
`us_prepare.py`  | Download event and station metadata
`us_qopen.py`    | Download waveforms and invert each event for intrinsic and scattering Q
`us_imaging.py`  | Interpolate between station locations and plot results
`us_plot_fits.py`| Reprocess some events and plot fits
`us_compare_mags.py`| Compare estimated moment magnitudes with magnitudes from reference catalogs
`config_us.json` | Configuration file for `qopen` module
**Results**      |
`results.json`   | Determined intrinsic and scattering attenuation and site amplification at station locations
`results_uncorrected.json` | Same as `results.json`, but site amplifications are not corrected
`eventresults.json`  | Results for the analyzed earthquakes (source displacement spectrum, estimated seismic moment and moment magnitude)
`eventresults_uncorrected.json`  | Same as `eventresults.json`, before correction of site amplifications
`figs_paper/*.pdf`   | Figures used in the article
`figs_suppl/*.pdf`   | Figures used in the electronic supplement of the article
`tif/*.tif`          | Georeferenced TIF files
**Visualization**    |
`Qviewer/index.html` | Javascript application to explore all results comfortably. It allows to view results at different stations as a function of frequency.
`how_to.ipynb`       | Jupyter notebook to demonstrate how to use the JSON and georeferenced TIF files


## Notes

The scripts use the *Qopen* method. Please browse the corresponding [Qopen Github repository](https://github.com/trichter/qopen) and the reference mentioned therein for more information.

To run the scripts download or clone the repository and install Python3 and `cartopy gdal matplotlib numpy obspy scipy seaborn shapely shove qopen`. This can be easily achieved with anaconda:

```
conda config --apend channels conda-forge
conda create -n usattenuation python=3 cartopy obspy gdal seaborn shapely joblib statsmodels
source activate usattenuation
pip install shove qopen
```

After that, create the tmp folder. If you want to run the `us_prepare.py`, `us_qopen.py` or `us_plot_fits.py` scripts it is also necessary to create the data folder and subfolders:

```
mkdir tmp
mkdir data data/events data/stations data/waveforms
```

Before running `us_qopen.py` or `us_plot_fits.py`, you need to download the metadata with `us_prepare.py` script. The other scripts can be run independently and recreate the files in `results` folder and the images for the *Qviewer* application in `visualization` folder.
