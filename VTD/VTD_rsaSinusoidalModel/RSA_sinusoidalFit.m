%% Fitting a sine wave to a single-beat Respiratory Sinus Arrhythmia data
function [A, B, phi, y, SSres, SStot, Rsq, phases_clean, dIBI_clean, thetas_fit, amplitudes_fit] = RSA_sinusoidalFit(phases,deltaIBI)

    %% INPUTS
    % phases = ...     % Nx1 vector of breathing phase (radians)
    % deltaIBI = ...     % Nx1 vector ΔIBI...
    %                     ... ΔIBI = (the R-R interval that starts with your target beat) - (the R-R interval that ends with your target beat)

    %% OUTPUTS
    % A = amplitude of the fitted sine wave (distance between the smallest and largest point on the y-axis)
    % B = vertical offset of the fitted sine wave
    % phi = phase offset of the fitted sine wave (horizontal offset)
    % y = value series fitted wave (predicted values)
    % 

    % phases_clean = inputted phase values after cleaning for outliers (threshold 3*SD)
    % amplitudes_clean = inputted ΔIBI values after cleaning for outliers (threshold 3*SD)
    % thetas_fit = fitted phase values
    % amplitudes_fit = fitted ΔIBI values

    %% sinusoidal model fitting
    % Outlier exclusion
    z_amp = zscore(deltaIBI);
    keep_idx = abs(z_amp) < 3;  % threshold: 3 std deviations
    phases_clean = phases(keep_idx);
    dIBI_clean = deltaIBI(keep_idx);

    % Define model function
    sinusoid = @(b, theta) b(1) * sin(theta + b(2)) + b(3);  % b = [A, phi, B]
    
    % Define loss function for fitting
    lossfun = @(b) sum((dIBI_clean - sinusoid(b, phases_clean)).^2);
    
    % Initial guess: [A, phase shift, offset]
    b0 = [max(dIBI_clean) - min(dIBI_clean), 0, mean(dIBI_clean)];
    
    % Optimize
    b_fit = fminsearch(lossfun, b0);
    
    % Predicted values
    thetas_fit = linspace(0, 2*pi, 1000);
    amplitudes_fit = sinusoid(b_fit, thetas_fit);

    A = b_fit(1);
    phi = b_fit(2);  % in radians
    B = b_fit(3);

    y = sinusoid(b_fit, phases_clean);
    SSres = sum((dIBI_clean - y).^2);
    SStot = sum((dIBI_clean - mean(dIBI_clean)).^2);
    Rsq = 1 - SSres / SStot;

end
