%% Frequency Evolution – Experimental

close all
clearvars

%% Signal Processing Inputs

fs = 2048*4;                        % Sampling frequency (Hz)
dt = 1/fs;                          % Time step (s)
nfft = 8192;

%% High-pass Filter: Drift Removal

fc = 2;                             % Cut-off frequency (Hz)
[b,a] = butter(4,fc/(fs/2),'high');

%% Import and Process Experimental Data

for n = 1:3

    if n == 1
        output(n).disc = "No";
        data = tdmsread("Submission\Experimental Data\Harmonic\Sine_39Hz_A0_1.8_100s.tdms");
        %data = tdmsread("Experimental Data\Swept\Swept_Sine_10-200Hz_Sr_2_A0_0.7.tdms");

    elseif n == 2
        output(n).disc = "Borderline";
        data = tdmsread("Submission\Experimental Data\Harmonic\Sine_39Hz_A0_1.91_100s.tdms");
        %data = tdmsread("Experimental Data\Swept\Swept_Sine_10-200Hz_Sr_2_A0_1.5.tdms");

    elseif n == 3
        output(n).disc = "Persistent";
        data = tdmsread("Submission\Experimental Data\Harmonic\Sine_39Hz_A0_1.98_100s.tdms");
        %data = tdmsread("Experimental Data\Swept\Swept_Sine_10-200Hz_Sr_2_A0_1.8.tdms");
    end

    [vel,base_accel,drive_sig] = extract_struct(data);

    % Time vector.
    t_series = (0:length(vel)-1)' * dt;

    % High-pass filter velocity.
    vel_hp = filtfilt(b,a,vel);

    % Integrate velocity to displacement.
    disp = cumtrapz(t_series,vel_hp);

    % Second high-pass filter to remove residual drift.
    disp = filtfilt(b,a,disp);

    % Remove residual offset.
    disp = disp - mean(disp);

    output(n).t_series = t_series;
    output(n).vel = vel_hp;
    output(n).disp = disp;
    output(n).acc_mag = max(abs(base_accel));

    % Velocity STFT.
    [s_vel,f,t_welch] = stft(vel,fs, ...
        Window=hann(nfft), ...
        OverlapLength=nfft/2, ...
        FFTLength=nfft, ...
        FrequencyRange="onesided");

    output(n).s_vel_abs = 20*log10(abs(s_vel));
    output(n).f = f;
    output(n).t = t_welch;

    % Base acceleration STFT.
    [s_acc,f,t_welch] = stft(base_accel,fs, ...
        Window=hann(nfft), ...
        OverlapLength=round(0.5*nfft), ...
        FFTLength=nfft, ...
        FrequencyRange="onesided");

    output(n).s_accel_abs = 20*log10(abs(s_acc));

end

%% Excitation Spectrograms

for n = 1:3

    figure
    imagesc(output(n).t,output(n).f,output(n).s_accel_abs)

    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    xlim([0 100])
    ylim([0 200])
    clim([-50 110])

    set(gcf,'Units','centimeters','Position',[0,0,65/4,45/4])
    % colorbar
    % saveas(gcf,'Experimental_FINAL_' + output(n).disc + '_Harmonic_Excitation.png')

end

%% Snap-through Spectrograms

for n = 1:3

    figure('Name',output(n).disc + " Snap-through")
    imagesc(output(n).t,output(n).f,output(n).s_vel_abs)

    axis xy
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    xlim([0 100])
    ylim([0 200])
    clim([-150 60])

    set(gcf,'Units','centimeters','Position',[0,0,65/4,45/4])
    % colorbar
    % saveas(gcf,'Experimental_FINAL_' + output(n).disc + '_Harmonic.png')

end

%% Response Time Histories

for n = 1:3

    figure
    T = tiledlayout(2,1);

    nexttile
    plot(output(n).t_series,output(n).vel)
    ylabel("Velocity (m/s)")
    xlabel("Time (s)")
    xlim([0 100])

    nexttile
    plot(output(n).t_series,output(n).disp*1000)
    ylabel("Position (mm)")
    yline(3.84)
    yline(-3.84)
    xlabel("Time (s)")
    xlim([0 100])

    set(gcf,'Units','centimeters','Position',[0,0,65/4,45/4])
    % saveas(gcf,'Experimental_Time_FINAL' + output(n).disc + '_Harmonic.png')

end

%% ------------------------------------------------------------------------
%% Functions
%% ------------------------------------------------------------------------

function [vel,base_accel,drive_sig] = extract_struct(data)

    %% Accelerometer Correction

    correction_ratio = 10.34/10.22;
    base_accel = data{1,1}.Accelerometer * correction_ratio;

    %% Vibrometer Conversion

    % DAQ recorded in mm/s; convert to SI base units (m/s).
    vel = data{1,1}.Vibrometer * 1e-3;

    %% Extract Valid Indices

    % Remove NaNs based on the drive signal.
    drive_sig = data{1,1}.("Drive signal");

    idx = ~isnan(drive_sig);
    cutoff_idx = 4*find(isnan(drive_sig),1,'first');

    vel = vel(1:cutoff_idx);
    base_accel = base_accel(1:cutoff_idx);

    drive_sig = drive_sig(idx);

    %% Remove Gross Mean/DC Offset

    vel = detrend(vel);
    vel = vel - mean(vel);

    base_accel = detrend(base_accel);
    base_accel = base_accel - mean(base_accel);

    drive_sig = detrend(drive_sig);
    drive_sig = drive_sig - mean(drive_sig);

end
