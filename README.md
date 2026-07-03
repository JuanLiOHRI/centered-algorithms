# Centering, recalibrating, and calibration curves (2026-07-03)

## Introduction

This repo produces a pipeline for 

1. Centering the following models (including __non-linear terms__: polynomial, restricted cubic spline, and interaction): 
    - linear regression
    - logistic regression
    - survival models: Cox proportional hazards model, and Fine-Gray subdistribution hazard model

    In brief, the workflow is to fit the original model as usual, and then pass the fitted model object and the derivation dataset into `center_fit` to transform data and perform centering that follows the _"transform then center"_ principle.

2. Make sure the model output
    - is most light-weighted without individual outcome or predictor data (_i.e._ `x = FALSE, y = FALSE` when fitting the models), this, to satisfy the strict privacy requirements when working in a secured environment (e.g., ICES).
    - but is also completed that the model can be applied to and evaluated using external dataset.
    
    Several functions for survival models have been modified to meet these two requirements.

## Linear regression: `1_Example_center_recalibrate_linear.qmd`

Detailed example using the [Medical Cost Personal Datasets](https://www.kaggle.com/datasets/mirichoi0218/insurance) dataset. 

✅ __All tasks have been completed here.__

- [x] Model fitting and centering (`center_fit`, which calls `center_prepare_0` and several utility functions)
- [x] Model transfering and predicting on the target data (`center_predict`, which calls `center_prepare`)
- [x] Model recalibrating (both math and code are clear): updating both predicter means for centerings and outcome means as the intercept (`center_recalibrate`, which calls `center_prepare`)
- [x] Calibration curves and metrics: may not be necessary for linear regression, but nonetheless done. (the wrapper function `calibration`, which calls `calibration.con`; and `calibration.OvsP` for the barplot of subgroups)
- [x] Calibration results confirm that recalibration improves model performance on the target dataset.

## Logistic resression: `2_Example_center_recalibrate_logistic.qmd`

Detailed example using the [Heart Disease Dataset](https://www.kaggle.com/datasets/johnsmith88/heart-disease-dataset) dataset. 

✅ __All tasks have been completed here.__

- [x] Model fitting and centering (`center_fit`, which calls `center_prepare_0` and several utility functions)
- [x] Model transfering and predicting on the target data (`center_predict`, which calls `center_prepare`)
- [x] Model recalibrating (both math and code are clear): only update the intercept by a shift `delta_recal = logit(outcome_mean_target) - logit(outcome_mean_original)` (`center_recalibrate`, which calls `center_prepare`)
- [x] Calibration curves and metrics: (the wrapper function `calibration`, which calls `CalibrationCurves::val.prob.ci.2` or `rms::val.prob`; and `calibration.OvsP` for the barplot of subgroups)
- [x] Calibration results confirm that recalibration improves model performance on the target dataset.

## Cox proportional hazards model

### Notes for `survival::coxph`

1. In R’s `survival::coxph` function, predictors are internally mean-centered for numerical stability. However, dummy variables are not automatically centered. We can ensure that dummy variables are also centered by setting `nocenter = NULL`. Therefore, the centering pipeline (`center_fit` + `center_predict`) is not strictly required. However, it is still recommended to use these centering functions to make sure all essential information are outputted to work with the modified function below.

2. For recalibration, unlike linear and logistic regression models, the `coxph` model object (and the baseline hazard) can not be directly modified. Thus the function `center_predict` is used to implement the model shift `delta_recal`.

3. The commonly-used function `CalibrationCurves::valProbSurvival` requires `x=TRUE, y=TRUE` and uses the fitted `coxph` model object. Therefore, it needs to be modified to work with the recalibrated model.

4. Further, the function `riskRegression::Score` called within `CalibrationCurves::valProbSurvival` for brier scores also uses the fitted `coxph` model object. For simplicity, I decides to instead modify the `survival::brier` function (again, calls the fitted `coxph` model object) to achieve almost the same results (`survival::brier` doesn't produce confidence interval).

5. Calibration and recalibration of the survival model requires a __specific time point__. When working inside the secured environment, one may want to prespecify a `times` vector when running `center_fit` to collect the "outcome mean in the derivation dataset" at various time horizon, which is the event probability based on the Kaplan-Meier curve at specific time.

### Functions modified

- [x] `survival::brier` --> `brier.cox`: calls the appropriate prediction and outcome instead of the fitted `coxph` model object 
- [x] `CalibrationCurves::valProbSurvival` --> `valProbSurvival.2`: calls the appropriate prediction and outcome instead of the fitted `coxph` model object and __doesn't__ require the model fitted with `x=TRUE, y=TRUE`.
- [x] `3.1_Test_new_functions_Cox.qmd` __validates the modified functions__ by comparing their outputs with the ones generated by their counterparts.

### Example: `3_Example_center_recalibrate_cox.qmd`

Detailed examples using the Breast Cancer Survival Data from Rotterdam and Germany (see `?CalibrationCurves::trainDataSurvival`). 

✅ __Almost all tasks have been completed here.__

- [x] Model fitting and centering (`center_fit`, mainly to update with `nocenter = NULL` and collect some essential information about the derivation dataset. Requires `times` for outcome mean in the derication dataset)
- [x] Model transfering and predicting on the target data (`center_predict`, which calls `center_prepare`)
- [x] Model recalibrating (both math and code are clear): only update the intercept by a shift `delta_recal = clog_log(outcome_mean_target) - clog_log(outcome_mean_original)` (`center_recalibrate`, which calls `center_prepare`)
- [x] Calibration curves and metrics: (the wrapper function `calibration`, which calls `valProbSurvival.2`; and `calibration.OvsP` for the barplot of subgroups)
- [ ] From the calibration results, recalibration makes performance worse than the original model (it might be a bit overshoot) at almost all time points that tested.

## Fine-Gray subdistribution hazard model

### Notes for `survival::finegray` + `survival::coxph`

1. Function `survival::finegray` creates a weighted, expanded dataset `data_fg`. Then, although `survival::coxph` with `nocenter = NULL` still does the internal centering, the mean values used are of the expanded `data_fg`, instead of the original `data` as we intended. Therefore, the centering pipeline (`center_fit` + `center_predict`) __must__ be used.

2. Function `CalibrationCurves::valProbSurvival` (and the modified `valProbSurvival.2`) doesn't work with the fine-gray model. Further develpment/modification is needed.

3. The alternative `riskRegression::FGR` (wrapper of `cmprsk::crr`), `riskRegression::Score`, and `riskRegression::plotCalibration` don't work with the fine-gray model fitted using the `survival` package either. In fact, there is a huge methodology difference between `survival::finegray` + `survival::coxph` vs `riskRegression::FGR`. My previous effort to make them work together led to nowhere. 

    - Also, `riskRegression::FGR` has some syntax issues with `poly`, `rcs`, and interaction.

    - Outputs of `survival::finegray` + `survival::coxph` and `riskRegression::FGR` are close but not identical (because the difference in their underlying methods).

4. 🚧 The current plan is to further modify `valProbSurvival.2` --> `valProbSurvival.fg` (the brier function may need to be modified as well). I'll fit a simple model using the `riskRegression` functions, and use their output as the benchmark to validate the `valProbSurvival.fg` function (in `4.1_Test_new_functions_FG.qmd`).

### Example: `4_Example_center_recalibrate_cox_fine_gray.qmd`

Detailed examples using the Monoclonal gammopathy data (see `?survival::mgus2`).

🚧 __Still Work In Progress__

- [x] Model fitting and centering (`center_fit`, which calls `center_prepare_0` and several utility functions. Requires `times` for outcome mean in the derication dataset)
- [x] Model transfering and predicting on the target data (`center_predict`, which calls `center_prepare`)
- [?] Model recalibrating: math and code should be similar to cox but may require some modifications
- [ ] Calibration curves and metrics

## The `/doc` folder

Discussion note and mathematical details.

## The `/old` folder

Backup files for older versions. Can be ignored.