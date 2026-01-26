# Centered algorithm and (re)calibration

Examples and helper functions for centered algorithm and (re)calibration.

## Main example files

### `1.1 linear_regression_simulated_data.qmd`

Detailed examples and some visualizations for properly centering variables in a **linear regression** model using simulated data, variable types including:

1.  Single dichotomous predictor
2.  Single categorical predictor (after transforming it into dummy variables)
3.  Single continuous predictor as it is
4.  Single continuous predictor using restricted cubic spline (rcs), including introduction of rcs
5.  Interaction of two categorical variables
6.  Interaction of one categorical variable and one continuous variabl (no rcs)
7.  Interaction of two continuous variables (no rcs)
8.  Interaction of one categorical variable and one continuous variabl (with rcs)
9.  Interaction of two continuous variables (one with rcs, one without)
10. Interaction of two continuous variables (with rcs)

### `1.2 linear_regression_simulated_data_recalibration.qmd`

An example using simulated data where a linear regression model (original and the centered version) is developed using a `development dataset`, and then applied and recalibrated on a `test dataset`. **This is the main code example for doing proper centering and recalibration**, including using the helper functions in `R/util.R` and `R/calibration.R`.

The following helper functions are used in this example:

-   **`R/util.R`**: `step_dummy`,`step_rcs`,`step_interaction`,`get_mean`,`step_center`
-   **`R/calibration.R`**: `calibration`

For details about these functions, see below.

### `1.3 linear_regression_real_data.qmd`

An example using real data on centering and recalibration of a linear regression model.

-   Dataset: [Medical Cost Personal Datasets](https://www.kaggle.com/datasets/mirichoi0218/insurance) from Kaggle. The dataset has been downloaded and stored in the `/data` folder. Please refer to the `.qmd` file for details.
-   Helper functions used: same as in the previous example.

### `2.1 logistic_regression_real_data.qmd`

An example using real data on centering and recalibration of a logistic regression model.

-   Dataset: [Heart Disease Dataset](https://www.kaggle.com/datasets/johnsmith88/heart-disease-dataset) from Kaggle. The dataset has been downloaded and stored in the `/data` folder. Please refer to the `.qmd` file for details.
-   Helper functions used: same as in the previous example.
-   Note to **Doug**: I have used a different way to split the data to exaggerate the difference between development and external dataset, see the last plot.

### `3.1 Cox_real_data.qmd`

An example using real data on centering and recalibration of a logistic regression model.

-   Dataset: Breast Cancer Survival Data from Rotterdam and Germany, see `?CalibrationCurves::trainDataSurvival`. Please refer to the `.qmd` file for details.
-   Helper functions used: same as in the previous example.

## The `/data` folder

-   `data/insurance.csv`: The [Medical Cost Personal Datasets](https://www.kaggle.com/datasets/mirichoi0218/insurance) from Kaggle, see `1.3 linear_regression_real_data.qmd`.

-   `data/heart.csv`: The [Heart Disease Dataset](https://www.kaggle.com/datasets/johnsmith88/heart-disease-dataset) from Kaggle, see `2.1 logistic_regression_real_data.qmd`.

## The `/R` folder

Some helper functions.

### `R/util.R`

-   `get_rcs`: Implement the formula of rcs components
-   `step_dummy`: Get dummy variables
-   `step_rcs`: Unpack the rcs components
-   `step_interaction`: Get interaction terms
-   `step_center`: Center all variables
-   `get_mean`: Get mean values
-   `get_logLik`: Get log likelihood
-   `root.search`: Root searching, not used in current examples

### `R/calibration.R`

-   `plot.calibration`: Recreate the calibration plot and summary stats from `rms::val.prob` for comparing calibration with different shifts.
-   `calibration`: A wrapper function to call different functions for calibration plots
-   `calibrationCon`: Generate calibration plot for continuous outcome (linear regression)

## The `/doc` folder

Discussion note and mathematical details.

## The `/old` folder

Backup files for older versions. Can be ignored.