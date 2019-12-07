% calculate confidce intervals of low folate experiment
[result_ratio_best_fit_LF result_ratio_confidece_intervals_LF] = compute_SHMT1_vs_SHMT2_confidence_intervals('LF');
% calculate confidce intervals of high folate experiment
[result_ratio_best_fit_HF result_ratio_confidece_intervals_HF] = compute_SHMT1_vs_SHMT2_confidence_intervals('HF');
% print results
fprintf('best fit ratio(LF)=%f  confidence intervals(LF)=[%f %f]\n', result_ratio_best_fit_LF, result_ratio_confidece_intervals_LF(1), result_ratio_confidece_intervals_LF(2));
fprintf('best fit ratio(HF)=%f  confidence intervals(HF)=[%f %f]\n', result_ratio_best_fit_HF, result_ratio_confidece_intervals_HF(1), result_ratio_confidece_intervals_HF(2));