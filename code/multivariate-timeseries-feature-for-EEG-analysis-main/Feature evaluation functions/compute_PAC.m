function [ft, ft_dscrp] = compute_PAC(y1, y2,filters)

% This function computes the PAC features for a pair of signals (y1, y2)
% using the provided frequency bands and sampling frequency.
% It returns the PAC features (ft) and their descriptions (ft_dscrp).
%y1 = y1(2:end) - y1(1:end-1);
%y2 = y2(2:end) - y2(1:end-1);

y_filt_delta = filtfilt(filters{1}, y1);

y_filt_theta = filtfilt(filters{2}, y1);

y_filt_alpha = filtfilt(filters{3}, y1);


y_delta_analytic = hilbert(y_filt_delta);
y_theta_analytic = hilbert(y_filt_theta);
y_alpha_analytic = hilbert(y_filt_alpha);

y2_filt_theta = filtfilt(filters{2}, y2);
y2_filt_alpha = filtfilt(filters{3}, y2);
y2_filt_beta = filtfilt(filters{4}, y2);

y2_theta_analytic = hilbert(y2_filt_theta);
y2_alpha_analytic = hilbert(y2_filt_alpha);
y2_beta_analytic = hilbert(y2_filt_beta);

amp_theta = real(y2_theta_analytic);
amp_alpha = real(y2_alpha_analytic);
amp_beta = real(y2_beta_analytic);

phase_delta = angle(y_delta_analytic);
phase_theta = angle(y_theta_analytic);
phase_alpha = angle(y_alpha_analytic);

idx_cycle_delta = find_cycles(phase_delta);
idx_cycle_theta = find_cycles(phase_theta);
idx_cycle_alpha = find_cycles(phase_alpha);

zt_DT = amp_theta(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
zt_DA = amp_alpha(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
zt_DB = amp_beta(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
zt_TA = amp_alpha(1:idx_cycle_theta).*exp(j*phase_theta(1:idx_cycle_theta));
zt_TB = amp_beta(1:idx_cycle_theta).*exp(j*phase_theta(1:idx_cycle_theta));
zt_AB = amp_beta(1:idx_cycle_alpha).*exp(j*phase_alpha(1:idx_cycle_alpha));

DTcouple = abs(sum(zt_DT)/sum(amp_theta(1:idx_cycle_delta)));
DAcouple = abs(sum(zt_DA)/sum(amp_alpha(1:idx_cycle_delta)));
DBcouple = abs(sum(zt_DB)/sum(amp_beta(1:idx_cycle_delta)));
TAcouple = abs(sum(zt_TA)/sum(amp_alpha(1:idx_cycle_theta)));
TBcouple = abs(sum(zt_TB)/sum(amp_beta(1:idx_cycle_theta)));
ABcouple = abs(sum(zt_AB)/sum(amp_beta(1:idx_cycle_alpha)));
ft = [DTcouple DAcouple DBcouple TAcouple TBcouple ABcouple];
ft_dscrp = {'DTcouple', 'DAcouple', 'DBcouple', 'TAcouple', 'TBcouple', 'ABcouple'};

end

function idx = find_cycles(phase)
stp = phase(1);
for i = length(phase):-1:2
    if phase(2)>stp
        if phase(i) >= stp && phase(i-1)<= stp
            idx = i;
            return
        end
    end
    if phase(2)<stp
        if phase(i) <= stp && phase(i-1)>= stp
            idx = i;
            return
        end
    end
end
end