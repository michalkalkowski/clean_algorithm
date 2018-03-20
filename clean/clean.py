#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==============================================================================
Copyright (C) 2018 Michal Kalkowski (MIT License)
kalkowski.m@gmail.com

This is an implementation of the CLEAN algorithm which identifies individual
waveforms in a multi-component signal by simple spectral summations.

The details of the algorithm can be found in (and are not recalled here):
[1] Gough, P.T., 1994. A fast spectral estimation algorithm based on the FFT.
    IEEE Transactions on Signal Processing 42, 1317–1322.
    https://doi.org/10.1109/78.286949
[2] Holmes, C., Drinkwater, B.W., Wilcox, P.D., Post-processing of the full
    matrix of ultrasonic transmit–receive array data for non-destructive
    evaluation, 2005. NDT & E International 38, 701–711.
    https://doi.org/10.1016/j.ndteint.2005.04.002
[3] Hunter, A.J., Drinkwater, B.W., Zhang, J., Wilcox, P.D., 2011.
    A STUDY INTO THE EFFECTS OF AN AUSTENITIC WELD ON ULTRASONIC ARRAY
    IMAGING PERFORMANCE. pp. 1063–1070. https://doi.org/10.1063/1.3592054

Two functions are provided. `extract_CLEAN` performs the iterative search and
identifies individual components (wave packets) present in the measured signal.
`plot_components` can be used to visualise the outcome of `extract_CLEAN`.
==============================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert


def extract_CLEAN(measured_signal, original_signal, delta_t, threshold=0.4):
    """
    Applies the CLEAN algorithm to extract individual components from a multi-
    component signal. The algorithm is based in the assumption that the
    measured signal is a sum of scaled, delayed and shifted copies of
    the original signal.

    Individual components are extracted by iterating over dominant components
    of the spectrum. The algorithms starts from taking the spectrum of the
    input signal as a starting residue. It determines the amplitude, phase
    and time delay of the dominant component and subtracts the reconstructed
    spectrum of this componend from the residue. A new residue is formed,
    and the search for the dominant component restarts.

    The loop is terminated after the amplitude associated with one individual
    component drops below a given threshold.

    It is assumed that both the original and the measured signals are sampled
    at the same rate.

    Parameters:
    ---
    measured_signal: ndarray, measured signal to decompose
    original_signal: ndarray, original (transmitted) signal,
                              e.g. the excitation applied to the structure
                              under test
    delta_t: float, time increment
    threshold: float, amplitude threshold as a fraction of the maximum
                      amplitude;
                    defaults to 0.4

    Returns:
    ---
    amplitudes: ndarray, amplitudes of the individual components
    delays: ndarray, tiem delays of the individual components
    phases: ndarray, phase shifts of the individual components
    components: ndarray, reconstructed time traces related to each individual
                        component
    """
    original_spectrum = np.fft.fft(original_signal)
    measured_spectrum = np.fft.fft(measured_signal)
    omega = np.fft.fftfreq(len(original_signal), delta_t)
    original_hilbert = hilbert(original_signal)
    original_argpeak = np.argmax(abs(original_hilbert))

    # Initialise lists
    amplitudes = []
    delays = []
    phases = []
    components = []

    # First residue is the spectrum of the measured signal
    residue = np.copy(measured_spectrum)

    # Dummy values initialising the loop
    amplitude = 1
    amp_max = 2

    while amplitude > threshold*amp_max:
        r_t = np.fft.ifft(residue)
        # Remove numerical noise
        r_t = r_t.real
        # Find the dominant wave packet
        r_hilbert = hilbert(r_t)
        r_argpeak = np.argmax(abs(r_hilbert))
        # Extract its amplitude
        amplitude = abs(r_hilbert[r_argpeak])
        # Extract time delay with reference to the envelope of the original
        # signal
        delay = (r_argpeak - original_argpeak)*delta_t
        # Extract phase shift with reference to the original signal
        phase = (np.angle(r_hilbert[r_argpeak]) -
                 np.angle(original_hilbert)[original_argpeak] + 2*np.pi)

        amplitudes.append(amplitude)
        delays.append(delay)
        phases.append(phase)
        # Apply the delay and shift only to the positive half of the spectrum
        if len(original_spectrum) % 2 == 0:
            nyquist_index = int(len(original_spectrum)/2) + 1
        else:
            nyquist_index = int(np.ceil(len(original_spectrum)/2))
        half_spectrum = original_spectrum[:nyquist_index]
        component_spectrum = (amplitude*half_spectrum
                              * np.exp(1j*(-2*np.pi*omega[:nyquist_index]*delay
                                           + phase)))
        # Reconstruct the double-sided spectrum
        if len(original_spectrum) % 2 == 0:
            reconstructed = np.concatenate((
                component_spectrum, component_spectrum.conj()[1:-1][::-1]))
        else:
            reconstructed = np.concatenate((
                component_spectrum, component_spectrum.conj()[1:][::-1]))
        residue = residue - reconstructed
        components.append(np.fft.ifft(reconstructed))
        # Assign maximum amplitude if this is the first iteration
        if len(amplitudes) == 1:
            amp_max = amplitude
    return (np.array(amplitudes[:-1]), np.array(delays[:-1]),
            np.array(phases[:-1]), np.array(components[:-1]))


def plot_components(time, measured_signal, components):
    """
    Plots the measured signal and individual components extracted using
    the CLEAN algorithm.

    ###
    Example
    ###
    >> amps, dls, phss, comps = extract_CLEAN(
                    measured_signal=signal, original_signal=transmitted,
                    delta_t=delta_t, threshold=0.4)
    >> plot_components(time, signal, comps)

    Parameters:
    ---
    time: ndarray, time vector
    measured_signal: ndarray, multi-component measured signal
    components: ndarray, the output of the extract_CLEAN function
    """
    plt.figure()
    plt.plot(time*1e6, measured_signal, c='lightgray', lw=4, label='measured')
    for i, comp in enumerate(components):
        plt.plot(time*1e6, comp.real, label='component {}'.format(i))
    plt.xlabel('time in us')
    plt.ylabel('normalised amplitude')
    plt.legend()
