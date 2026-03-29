# NRMM-FRET-hist-fit
The MATLAB script `script_hist_fit.m` is designed to fit single molecule FRET histograms with a **realistic-noise mixture model (NRMM)**.
Adjust user parameters in the section "USER PARAMETERS" prior running the script.

## Input
Histogram data must first be formatted to a suitable .mat file. Run script `import_from_txt.m` to convert a text file to a suitable .mat file prior runing this script. 

Column **1** of the text file must contains the **bin centers** and the column **2** the **bin counts**.

## Histogram peak model
The model used to fit one FRET histogram peak includes:
- the broadening inherent to the **flexibility of the molecular structure**,
- the broadening induced by the **experimental noise** including photon emission and camera detection.

Due to the complexity - but high fidelity - of the model, we use Monte Carlo (MC) simulations to pre-calculate a look-up table (LUT) of probability masses for fixed grids of parameters.
Model parameters are the **FRET mean** and the **structure flexibility** of the peak, both varying between 0 and 1. 

Model constants, however, include many experimental mensurations such as donor and acceptor **background photon counts**, **mean donor emission** in the absence of acceptor, **detection parameters** and **correction factors** *[cross-talk coefficients will be added in the future]*.

### Probability look-up table
Practically, the LUT is calculated once for the parameters grids and saved to the analysis `.mat` file.
Calculation progress shown in the prompt as:
```
grid index xxxx / xxxx done
```
If the model constants are changed in `` or ``, the table will be recalculated at the next run.

Below is shown a heatmap of the slice $Flx=0.5$ of an example LUT.
Probability masses are calculated for apparent FRET values ranging from 0 to 1.

<img width="370" height="241" alt="image" src="https://github.com/user-attachments/assets/b1fe6733-38e0-4d4b-bf23-e4f2e36b1c36" />

### Structure flexibility
The flexibility parameter maps the standard deviation from the equilibrium distance between the FRET pair between 0 (rigid) and 1 (loose).
This is achieved by using a sigmoïdal relationship:

$$
Flx = \frac{\sigma_r^k}{S^k + \sigma_r^k},
$$

where $S$ and $k$ are pre-defined constants.

For $S=2$ and $k=2$, the flexibility reaches 0.5 for a standard deviation of 2Å and asymptotically converges towards 1.

<img width="359" height="229" alt="image" src="https://github.com/user-attachments/assets/17c386e5-bd94-4c0e-baa1-44948166b355" />

## Inference methods
Three inferrence methods are available:
- **Bayesian non-parametric inferrence** of NPMM using Dirchlet Process mixture models (DPMM),
- **Expectation-maximization (EM)** of NPMM.
- EM of **Gaussian mixture models (GMM)**.

When chosen, several runs of DPMM are executed in parallel and the best outcome is selected .
A figure is created for each run and exported to `.png` and `.fig` files once the inference is completed.

## Output
Inference results are saved to a `_results_[method].mat` file and a graphical summary is exported to `.png` and `.fig` image files, all in a subfolder named after the source file.

<img width="1169" height="947" alt="histogram_eD135e_L43E_100mM_pre_DPMM(4)" src="https://github.com/user-attachments/assets/0f534fc3-2580-4e52-8b50-2097e5a2e967" />

The graphical summary shows from left to right and top to bottom:
1. a slice (at $Flx=0.5$) of the LUT
1. the flexibility in function of the standard FRET distance deviation of the structure, 
1. the experimental FRET histogram with the fitted mixture, 
1. a superposition of 2D representations of the pre-calculated PMF of an histogram peak for some values of the FRET mean and at 0.5 flexibility,
1. the prior distribution of the model parameters (used in DPMM)
1. the experimental FRET histogram with the individual components of fitted mixture, 
1. iteration traces of the inferred number of peaks (K) and the posterior log-probability of the inferred model
1. distribution of the means of the inferred FRET peaks accross the iteration where K=Kopt
1. distribution of K accross the iterations
1. distribution of posterior log-prob accross the iteration
1. distribution of the standard distance deviation of the inferred peaks
1. distribution of the relative weight of the peaks in the inferred mixture.

## Parallel computing
For speed, LUT and DPMM calculations uses the parallel computing tool. It is highly recommended to install and use it.

