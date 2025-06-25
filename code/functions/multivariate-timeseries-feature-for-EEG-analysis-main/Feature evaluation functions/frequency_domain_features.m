function features = frequency_domain_features(signal)
    % Compute the power spectral density of the signal using Welch's method
    [psd,f] = pwelch(signal, [], [], [], 256);
    
    % Compute the total power and the average power in different frequency bands
    total_power = sum(psd);
    delta_power = sum(psd(f >= 0 & f <= 4));
    theta_power = sum(psd(f >= 4 & f <= 8));
    alpha_power = sum(psd(f >= 8 & f <= 13));
    beta_power = sum(psd(f >= 13 & f <= 30));
    gamma_power = sum(psd(f >= 30 & f <= 100));
    average_delta = delta_power / total_power;
    average_theta = theta_power / total_power;
    average_alpha = alpha_power / total_power;
    average_beta = beta_power / total_power;
    average_gamma = gamma_power / total_power;



    % Compute the spectral entropy
    spectral_entropy_val = entropy(psd);

    % Compute the spectral flatness
    spectral_flatness_val = geomean(psd) / mean(psd);
    
    % Concatenate the computed features into a single feature vector
    features = [average_delta average_theta average_alpha average_beta average_gamma spectral_entropy_val spectral_flatness_val];
end