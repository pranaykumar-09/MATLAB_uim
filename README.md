# MATLAB_uim


---

## **About: UIM Extraction from a Real-World Signal**

This repository implements a signal processing pipeline designed to extract **unintentional modulation (UIM)** components from a real-valued signal captured from a radio transmitter. These UIM components are subtle, device-specific imperfections embedded in the phase of a transmitted waveform, often caused by hardware noise, clock jitter, or local oscillator drift. The extracted UIM phase (UIMIP) can be used for device fingerprinting or signal authentication.



---

## **1. Signal Preprocessing**

The input is a `.dat` file containing real-valued samples. These values typically come from ADC recordings of RF signals. Because many RF capture systems introduce a strong DC offset, we first remove the mean value from the signal. This ensures that the Hilbert transform (used later to construct a complex signal) is not corrupted by a large constant bias. Centering the signal around zero is critical for clean phase estimation.

```matlab
x = x - mean(x);
```

---

## **2. Analytic Signal Construction (Hilbert Transform)**

The real-valued signal is then passed through the Hilbert transform to generate its **analytic signal**. This is a complex signal of the form:

$$
x_a[n] = x[n] + j \cdot \hat{x}[n]
$$

where $\hat{x}[n]$ is the 90° phase-shifted version of $x[n]$. This complex signal enables direct calculation of the **instantaneous phase** and frequency of the original waveform.

```matlab
analytic_signal = hilbert(x);
```

---

## **3. Sliding Window Segmentation**

The analytic signal is divided into overlapping time segments using a sliding window of fixed length $N_x$. Each segment is independently processed to estimate its local phase and frequency characteristics. This segmentation is crucial for tracking how phase and modulation evolve over time.

```matlab
for idx = 1:step:N - Nx + 1
    seg = analytic_signal(idx : idx + Nx - 1);
```

---

## **4. Unwrapped Phase Estimation**

For each segment, the **instantaneous phase** is computed using the angle of the complex samples. The `unwrap()` function removes artificial $\pm\pi$ discontinuities, yielding a smooth and continuous phase trajectory within the segment.

```matlab
phase_unwrapped = unwrap(angle(seg));
```

---

## **5. Modified Kay Frequency Estimation**

The **modified Kay algorithm** is used to estimate the local frequency from the unwrapped phase. This involves computing a weighted sum of the first-order phase differences using an optimal weight vector derived from the inverse of a tridiagonal covariance matrix. This approach improves the robustness of frequency estimation in the presence of phase noise.

```matlab
D = diff(phase_unwrapped).';
H = ones(Nx - 1, 1);
C = 2*eye(Nx - 1) - diag(ones(Nx - 2, 1), 1) - diag(ones(Nx - 2, 1), -1);
w_vector = inv(C) * H;
omega_est = sum(w_vector .* D) / sum(w_vector);
f_est = omega_est / (2*pi) * fs;
```

---

## **6. Maximum Likelihood Phase Estimation**

With the frequency estimate in hand, the **instantaneous phase offset** is estimated using a maximum likelihood approach. The idea is to fit a cosine waveform at the estimated frequency to the real part of the segment and use its projection to compute the phase angle. This gives a single phase estimate for each segment.

```matlab
n = 0:Nx-1;
re_part = real(seg);
num_phase = sum(re_part .* sin(2*pi*f_est/fs * n));
den_phase = sum(re_part .* cos(2*pi*f_est/fs * n));
IP_est(idx) = atan2(num_phase, den_phase);
```

---

## **7. Phase Curve Fitting (Intentional Modulation Removal)**

Once all segment-wise phase estimates are collected, the full phase vector is unwrapped again. A low-order polynomial (typically linear or quadratic) is fit to this phase trend using `polyfit`. This polynomial captures the **intentional phase modulation** — such as drift or carrier offset — that the transmitter meant to introduce.

```matlab
t_IP = (0:NT-1) / fs;
IP_unwrapped = unwrap(IP_est);
coeffs = polyfit(t_IP, IP_unwrapped, p_order);
IP_fit = polyval(coeffs, t_IP);
```

---

## **8. UIMIP Extraction (Residual Phase Calculation)**

The **residual phase**, referred to as the **UIMIP**, is computed by subtracting the fitted intentional phase from the total estimated phase. This residual contains the unintentional, hardware-specific imperfections — the true UIM signal.

```matlab
UIMIP = IP_unwrapped - IP_fit;
```

---



---

## **Summary**

This pipeline effectively isolates subtle device-level behaviors embedded in the transmitted signal phase. These residuals (UIMIP) can be used for applications such as:

* Device fingerprinting
* RF transmitter authentication
* Modulation integrity analysis

By working in the phase domain and leveraging carefully designed statistical estimators, this method delivers robust UIM extraction even in noisy conditions.

---

