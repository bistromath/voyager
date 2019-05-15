#!/usr/bin/env python

import scipy.signal
import scipy.fftpack
import numpy as np
import struct
import matplotlib.pyplot as plt
import math
from progress.bar import ChargingBar
import multiprocessing

samp_rate = 2929687.5
decim_rate = 10
start_time = 0
end_time = 60
nsecs = end_time-start_time #total amount of time to read from the file -- integration time
fftlen = 65536*4
num_integrations = int(samp_rate*nsecs/decim_rate/fftlen)
print "Number of integrations: %i" % num_integrations
veclen = fftlen - (fftlen % 4)
input_veclen = int(veclen * decim_rate)
print "Time per integration: %.2fs" % (veclen*decim_rate/samp_rate)
print "FFT length: %i" % veclen
bin_size = samp_rate/decim_rate/veclen
print "Bin size: %.2f" % bin_size

sig_start = 8420546732.5122
sig_freqslope = 0.37172 #signal drift rate in Hz/second, found experimentally
rec_ctr = 8419921875 #center frequency of the recording
offset = sig_start-rec_ctr #offset frequency at the start of the recording
#TODO: sigMF this
samplesize=1

filename = 'blc07_guppi_57650_67573_Voyager1_pol1.sigmf-data'

candidate_drifts = np.linspace(-1.0, 1.0, 20)

def try_acquisition(drift):
    #we'll store integrations here
    f = open(filename, 'rb')
    sumvector = np.zeros(veclen, dtype=np.float64)
    f.seek(int(start_time*samp_rate*samplesize*2))
    for i in range(num_integrations):
        vec = np.fromstring(f.read(input_veclen*2*samplesize), dtype=np.int8).astype(np.float64)/(2**(samplesize*8))
        cvec = vec.view(np.complex128)
        freqvec = np.linspace(offset+drift*i*input_veclen/samp_rate,
                            offset+drift*(i+1)*input_veclen/samp_rate,
                            input_veclen)
        freqshiftvec = np.exp(-2j*math.pi*freqvec*np.arange(input_veclen)/samp_rate)
        rvec = cvec * freqshiftvec
        dvec = scipy.signal.resample_poly(rvec, 1, decim_rate)
        W = scipy.signal.blackman(len(dvec))
        H = scipy.fftpack.fft(dvec*W, len(dvec))
        Ha = np.abs(scipy.fftpack.fftshift(H))
        sumvector += Ha

    Hl = 20*np.log10(sumvector) - 10*np.log10(num_integrations)
    peak = np.max(Hl)
    print "Drift %.2f: %.2fdB" % (drift, peak)
    return (drift, peak, Hl)

p = multiprocessing.Pool(4)
correlations = p.map(try_acquisition, candidate_drifts)
best = max(correlations,key=lambda x: x[1])
best_drift = best[0]
best_corr  = best[1]
Hl = best[2]
print "Best correlation: %.2fdB" % best_corr
print "Best drift rate: %.2f" % best_drift

freqs = np.arange(-samp_rate/decim_rate/2, samp_rate/decim_rate/2, samp_rate/decim_rate/fftlen)

plt.plot(freqs, Hl, color='blue', label="%.2f" % best_drift)
for p in correlations:
    plt.plot(freqs, p[2], alpha=0.2, label="%.2f" % p[0])
plt.ylabel('relative magnitude (dB)')
plt.xlabel('frequency (Hz)')
plt.legend()

#maxfreq = freqs[np.argmax(Hl)]
#print "Maximum frequency: %.4f" % (maxfreq + offset + rec_ctr)
#print "Maximum amplitude: %.4f" % np.max(Hl)
plt.show()
