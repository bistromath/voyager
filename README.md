<h1>Voyager 1 detection from recorded data</h1>
<h3>Nick Foster, 2019</h3>

The master notebook is voyager.ipynb.

This notebook details how to use time integration to locate and detect Voyager 1's signal from a recording. The example recording used here was recorded at the Allen Telescope Array in Hat Creek, CA on July 14, 2010. This recording (in SigMF format) can be found <a href="http://setiquest.info/sigmf/">here</a>. The notebook is general enough to be used with most any recording as long as you provide a valid SigMF header.

The Voyager 1 signal at 13,663,774,381 miles from Earth is far too weak to be plainly detectable in a spectrogram, even when using a radiotelescope array like the ATA, although world-class dishes such as the Goldstone observatory will be able to detect it strongly enough to see without further processing. Instead, to bring the signal out of the noise, we integrate over time. If we integrate for long enough, the power in the frequency bin corresponding to Voyager 1's carrier will increase enough to be detectable.

One problem with this approach is that we don't immediately know exactly where (in frequency) to look for the signal, due to <a href=https://en.wikipedia.org/wiki/Doppler_shift>Doppler shift</a>. There are three major sources of Doppler shift in this instance:

* Voyager 1's velocity relative to the Solar System (currently 38,026.77 mph)
* Voyager 1's velocity relative to Earth's orbit around the sun
* Voyager 1's velocity relative to the observatory as it rises and sets due to Earth's rotation

The first case, Voyager's velocity relative to the sun, is pretty static at this point -- no more thrust is being added to Voyager's trajectory, so it's not changing much (although interestingly there are some very small discrepancies from the expected velocity profile)

As Earth rotates around the sun once per year, its velocity relative to Voyager will change. Half the year it will be moving somewhat towards Voyager, and half the year away. Voyager's trajectory isn't collinear with Earth's orbital plane, so the magnitude of this effect is not exactly Earth's speed around the sun.

Viewed from an Earth observer's perspective, the Sun (or any other object in the sky) is moving toward you as it rises in the sky, and away from you as it sets. This effect causes a diurnal (once per day) change in Doppler shift, and also causes a further complication: the Doppler shift is likely to change over the course of your recording! If not compensated for, this will result in the energy coming from Voyager, which we want to be a nice spike in a single bin (from the carrier tone Voyager continually broadcasts), being smeared across multiple bins and thus much harder to detect.

We need to be able to calculate the magnitudes of the three velocities given above. Luckily, NASA JPL maintains a public-access tool called <a href="https://ssd.jpl.nasa.gov/?horizons">HORIZONS</a> to do exactly this. HORIZONS has accurate ephemerides (a set of numbers giving the position and velocity of an astronomical object) for a whole slew of man-made and astronomical objects. It also calculates relative velocities for any two objects, given the time of interest. It does a whole heck of a lot more than that, actually, but that's all we care about here. Even more luckily, there's a Python component for retrieving and parsing HORIZONS data.

So, we go out to HORIZONS, see where Earth and Voyager were relative to each other on July 14, 2010, and then we calculate an expected Doppler shift as well as the *drift* in that shift over the course of the recording. For the purposes of this notebook, we rely on a linear approximation of Doppler shift over the course of the recording. This is fine for recordings of an hour or so, but if you're really doing all-night observations of an astronomical target you probably want to account for this effect more precisely.

Once we have a Doppler shift profile, we go back to the original data and de-Doppler it so that hopefully the energy lies in a single bin. But, especially for long recordings and weak signals, this isn't quite accurate enough for good results! To maximize the height of the spike in our plot, which is always the goal of any communications engineer, we use the expected Doppler drift rate as the center point in a search pattern: we just try a bunch of candidate drift rates centered around the expected rate, and hope that one of them will be close enough to integrate cleanly without smearing the energy too much. If we've done things correctly, we should see the spike highest at or near the expected drift rate.

The following Python packages are necessary to run this notebook:

* Scipy
* Matplotlib/pyplot
* AstroPy (used for Julian date math)
* AstroQuery (used for downloading and processing ephemerides)
