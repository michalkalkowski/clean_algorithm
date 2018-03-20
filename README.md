CLEAN
==============================

An implementation and demonstration of the CLEAN algorithm for estimating times of arrival in multi component signals.

This is an implementation of the CLEAN algorithm which identifies individual waveforms in a multi-component signal by simple spectral summations.

The details of the algorithm can be found in (and are not recalled here):
[1] Gough, P.T., 1994. **A fast spectral estimation algorithm based on the FFT**. *IEEE Transactions on Signal Processing* 42, 1317–1322. https://doi.org/10.1109/78.286949
[2] Holmes, C., Drinkwater, B.W., Wilcox, P.D., **Post-processing of the full matrix of ultrasonic transmit–receive array data for non-destructive evaluation**, 2005. *NDT & E International* 38, 701–711. https://doi.org/10.1016/j.ndteint.2005.04.002
[3] Hunter, A.J., Drinkwater, B.W., Zhang, J., Wilcox, P.D., 2011. **A STUDY INTO THE EFFECTS OF AN AUSTENITIC WELD ON ULTRASONIC ARRAY IMAGING PERFORMANCE**. *Rev. QNDE* pp. 1063–1070. https://doi.org/10.1063/1.3592054

Two functions are provided. `extract_CLEAN` performs the iterative search and
identifies individual components (wave packets) present in the measured signal.
`plot_components` can be used to visualise the outcome of `extract_CLEAN`.

Project Organization
------------

    ├── README.md          <- The top-level README for developers using this project.
    ├── data               <- Original, raw, external data 
    │
    ├── clean/  <- Python module with the source
    code
    │
    ├── notebooks          <- Jupyter notebooks for expploration,
    demonstration, visualisation
    │
    └── setup.py           <- allows for installing the module
--------
