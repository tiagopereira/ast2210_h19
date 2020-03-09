# Preparation


## Software

To do this assignment you will need Python 3.x with [Astropy](https://www.astropy.org/), [Sunpy](https://www.sunpy.org/), and [IRISpy](http://docs.sunpy.org/projects/irispy/en/latest/) packages, plus Jupyter and many other dependencies of these packages.

**All the required software is already installed** in the Institute's Linux machines, to which you can connect to. See the  [guide](https://www.mn.uio.no/astro/english/services/it/help/programming/using-python.html) on how to use Python at ITA.

### Installing in your own laptop

Installing in your laptop is optional. You will need to ensure you have all the required software installed. The easiest way to do this is using the [miniconda](https://conda.io/miniconda.html) or [Anaconda](https://www.anaconda.com/download/) Python distributions (*python 3.x versions*). Miniconda is recommended because it is a smaller download, but Anaconda works just as well if you already have it or are more familiar with it.

Once you have conda installed (either a new install or an older version), the recommended way to install the packages is to create a new enviroment (we'll call it `ast2210`) to ensure you have the most recent versions. You can do this by:

``` bash
conda create -n ast2210 -c conda-forge --override-channels --yes \
  python=3.6 jupyterlab sunpy ndcube ipympl widgetsnbextension pip
```

This will download about 160 MB and install all the needed dependencies. Next, you need to activate this environment:

``` bash
source activate ast2210
```

!!! warning
    Every time you want to use the newly installed python packages, you must ensure you are running from the `ast2210` environment. Once active, your prompt will start with `(ast2210)`. If you open a new terminal, you will need to activate the environment again.

IRISPy has to be installed manually from Github. Run the following from the terminal:

``` bash
pip install git+git://github.com/sunpy/irispy.git@master
```

If using Jupyterlab, you will also need to install the extensions:

``` bash
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib
```

For running the examples and writing the assignment you should use Jupyter notebooks. A notebook is a document where you can combine text, images, and source code. If you are unfamiliar with Jupyter, there are [several](https://jupyter-notebook.readthedocs.io/en/stable/) [guides](https://www.youtube.com/watch?v=HW29067qVWk) and [tutorials](https://www.datacamp.com/community/tutorials/tutorial-jupyter-notebook). We recommend you use [Jupyterlab](https://jupyterlab.readthedocs.io/en/latest/), the next-generation version of Jupyter. But you can also use the classical notebook interface.

### Testing installation

To make sure you have all necessary software ready, start Jupyterlab from the terminal (after activating the `ast2210` environment):

```
jupyter lab
```

This will then open up a browser with the Jupyterlab launcher. Choose "Python 3" notebook, and it will start a new notebook. In the first cell enter the following and run:

``` python
import sunpy.map
import irispy.sji
```

!!! success
    If you got no error messages above, your installation is good and you are ready to start!

## Downloading data

To run the examples and complete this assignment you will have to download several hundreds of MB of data. Depending on your internet connection this may take a long time, so you should start early! You will need to download data from two space telescopes: AIA and IRIS.

!!! info "Astronomical file formats"
    All the data files we will use are in the [FITS format](https://en.wikipedia.org/wiki/FITS). FITS is an ancient file format developed for Astronomy, but is so standard that it is still in use by nearly all astronomical observatories. FITS files consist of a header with metadata about the observations, followed by a binary part with data. Depending on the type of observations, the sizes of FITS files can vary from a few KB to several TB. In the future, FITS will probably be replaced by the [ASDF](https://en.wikipedia.org/wiki/Advanced_Scientific_Data_Format) format.

### IRIS data

You will need to download data from two IRIS datasets:

* [2014.09.19 Emerging sunspots](http://sdc.uio.no/search/file/iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits) (18 MB), 2832 SJI, see also [IRIS event page](http://www.lmsal.com/hek/hcr?cmd=view-event&event-id=ivo%3A%2F%2Fsot.lmsal.com%2FVOEvent%23VOEvent_IRIS_20140919_051712_3860608353_2014-09-19T05%3A17%3A122014-09-19T05%3A17%3A12.xml)
* [2015.03.19 Mysterious event at limb](http://sdc.uio.no/search/file/iris_l2_20150319_090911_3860359580_SJI_1330_t000.fits) (572 MB), 1330 SJI, see also [IRIS event page](http://www.lmsal.com/hek/hcr?cmd=view-event&event-id=ivo%3A%2F%2Fsot.lmsal.com%2FVOEvent%23VOEvent_IRIS_20150319_090911_3860359580_2015-03-19T09%3A09%3A112015-03-19T09%3A09%3A11.xml)


### AIA data

There are many ways to download AIA data. As part of the exercises, you will download AIA data from Python, so no data links are provided beforehand. Moreover, it will be part of the exercises to decide how much AIA data to download. At some point you will be free to choose which channels to download, how regularly to sample, etc.

Each AIA file is about 10 MB. This does not seem a lot, but if you download one hour at maximum cadence for one channel, it adds to about 2.9 GB! So you will have to be careful with this.
