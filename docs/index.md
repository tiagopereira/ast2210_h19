# Introduction

Theses pages contain the notes and description of the Solar Lab assignment for AST2210. Please read carefully to make sure you understand the background materials, what is expected from this assignment, and the necessary software to run it.

## Required software

This assignment is primarily a computational task, and should be written using the Python programming language. While it is possible to complete most tasks using different programming languages, these notes are centred around Python libraries that simplify your work. To encourage good scientific practices, your report should be in the form of a Jupyter notebook (see [template](template.ipynb)), and you need to use [Jupyterlab](https://jupyterlab.readthedocs.io/en/latest/) or the [Jupyter notebook](https://jupyter.org/).

You will need the following software:

* [Python 3.x](https://docs.python.org/3/)
* [Jupyterlab](https://jupyterlab.readthedocs.io/en/latest/) / [Jupyter notebook](https://jupyter.org/)
* [Astropy](https://www.astropy.org/), [Sunpy](https://www.sunpy.org/), and [IRISpy](http://docs.sunpy.org/projects/irispy/en/latest/) packages

More details on obtaining and installing these packages can be found in the [Preparation](prep.md) section.

## How to read these notes

After the [Preparation](prep.md) section, these notes are written in a notebook fashion to encourage experimentation. Explanations of different procedures are included in between Python code blocks such as the following:

``` python
import numpy as np
print(np.arange(5))
```

You can copy paste the code into notebook cells and run the same commands to ensure you get the same results and get a feel for how it works. To make this task even easier, the sections are also available to download as a Jupyter notebook, which you can open in your computer and run directly.

We recommend you run through all the code examples before getting started on the assignment. The final section details the [assignment](assignment.md) questions and format.


## Data and Observatories

In this project you will work with observations of the Sun taken by two NASA space telescopes: the [Atmospheric Imaging Assembly](http://aia.lmsal.com) (AIA) and the [Interface Region Imaging Spectrograph](https://www.nasa.gov/mission_pages/iris/index.html) (IRIS). During your progress through these notes and for your assignment you will need to download hundreds of megabytes of data, so make sure you plan accordingly to ensure you don't run out of time while waiting for the data to download! The largest data volume will be downloaded before you start and is described in [Preparation](prep.md), but you will also download other that on the go, and the specific amount will depend on your curiosity.

Below is a brief introduction to the two telescopes we will use.

### Atmospheric Imaging Assembly

[AIA](http://aia.lmsal.com) is a part of the much larger [Solar Dynamics Observatory](https://sdo.gsfc.nasa.gov/) (SDO), which was launched in 2010 and is one of the largest and most expensive solar telescopes every launched. AIA is an array 4 telescopes with 20 cm apertures, which give an angular resolution of about 1″ (arcsec). It uses a variety of filters in parallel to gather 8 images of the Sun every 12 seconds. Each image covers the whole Sun and has a resolution of 4096x4096 pixels (16 Megapixels). As of 2018, AIA has collected more than 178 million images of the Sun.

!!! info "Angular vs. spatial resolution"
    In Astronomy, we use angular resolution as the measure of resolving power of a telescope. This is typically measured in arc seconds of a degree (″) and tells us the smallest aperture in the sky that the telescope can distinguish. To convert from angular to spatial resolution we need to know the distance to the far-away object. For the Sun, 1″ translates to about 725 km on the solar surface. This is the equivalent of seeing a human hair held 10 m away.

The different filters on AIA give us a view of different layers of the solar atmosphere. There is a lot of activity in the Sun, and by using narrow-band filters we can study different layers, or temperatures.


AIA channel  |	Source ion        |	Region of solar atmosphere              | Characteristic temperature
-----------: | ------------------ | --------------------------------------- | ----------------------------
    450.0 nm | continuum          | Photosphere 	                        | $5777$ K
    170.0 nm | continuum          | Temperature minimum, photosphere        | $5000$ K
     30.4 nm | He II 	          | Chromosphere & transition region        | $50\,000$ K
    160.0 nm | C IV + continuum   | Transition region & upper photosphere   | $10^5$ & $5000$ K
     17.1 nm | Fe IX 	          | Quiet corona, upper transition region   | $6.3×10^5$ K
     19.3 nm | Fe XII, XXIV       | Corona & hot flare plasma 	            | $1.2×10^6$ & $2×10^7$ K
     21.1 nm | Fe XIV             | Active region corona                    | $2×10^6$ K
     33.5 nm | Fe XVI 	          | Active region corona 	                | $2.5×10^6$ K
      9.4 nm | Fe XVIII 	      | Flaring regions 	                    | $6.3×10^6$ K
     13.1 nm | Fe VIII, XX, XXIII | Flaring regions 	                    | $4×10^5$, $10^7$ & $1.6×10^7$ K

The names of the AIA channels are often the wavelengths in Å (e.g. 171, 304, 1600). You can see live AIA images at the [Sun Today](http://suntoday.lmsal.com/) (and also look at archives).

### Interface Region Imaging Spectrograph

[IRIS](https://www.nasa.gov/mission_pages/iris/index.html) is a NASA small mission explorer that was launched in 2013. While much smaller than AIA, IRIS is able to study in detail the fine structure of the solar atmosphere, from the photosphere, chromosphere, and transition region. Like AIA, IRIS also has a 20 cm aperture telescope, but its angular resolution is 0″.33. This is achieved by having a much smaller field of view (maximum 175″x175″). Besides images, IRIS has a spectrograph that observes in two bands: far UV (133.1 - 140.7 nm) and near UV (278.2 - 283.5 nm). The IRIS imager is called Slit Jaw Imager (SJI) and has the following filters:

IRIS channel |	Source ion      | 	Region of solar atmosphere              | Characteristic temperature
-----------: | ------------------ | --------------------------------------- | ----------------------------
    283.2 nm | continuum          | Photosphere 	                        | $5777$ K
    279.6 nm | Mg II k            | Chromosphere                            | $10\,000$ K
    140.0 nm | Si IV 	          | Transition region                       | $80\,000$ K
    133.0 nm | C II + Fe XXI      | Transition region & flaring corona      | $15\,000$ & $10^7$ K

Although the most detailed information from IRIS comes from its spectra, for this project we will only work with SJI images.

The Institute of Theoretical Astrophysics has been a partner of the IRIS mission since its inception. Members of the solar group are involved with the interpretation of IRIS data, and the Institute has a local mirror of IRIS data, at the [Hinode Science Data Centre Europe](http://sdc.uio.no/sdc/). You can read an article [about our involvement with IRIS](https://www.mn.uio.no/astro/forskning/aktuelt/aktuelle-saker/astronytt/arkiv/2013/iris-science.html)  (in Norwegian). More information about IRIS on the [mission website](http://iris.lmsal.com/).
