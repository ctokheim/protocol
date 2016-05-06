# Protocol

## About

A major goal of the huge public investment in large-scale cancer sequencing has been to find the majority of driver genes.  Robust computational prediction of drivers from small somatic variants is critical to this mission, and it is essential that the best methods be identified.  While many such methods have been proposed, it has been difficult to evaluate them because there is no gold standard to use as a benchmark.  Here we developed an evaluation framework for driver gene prediction methods that does not require a gold standard.  The framework includes a large set of small somatic mutations from a wide range of cancer types and five evaluation metrics.  We propose it can be used to systematically evaluate new prediction methods and compare them to existing methods.  

## Installation

### Python Package Installation

Using the python package installation, all the required python packages will automatically be installed for you.

To install the package into python you can use `pip`. If you are installing to a system wide python then you may need to use `sudo` before the pip command.

```bash
$ pip install https://github.com/KarchinLab/protocol/archive/v1.0.0.tar.gz 
```

The scripts can then be found in `Your_Python_Root_Dir/bin`. You can
check the installation with the following:

```bash
$ which driver_protocol
$ driver_protocol --help
```

### Local installation

Local installation is a good option if you do not have privilege to install a python package and already have the required packages.  The source files can also be manually downloaded from github at https://github.com/KarchinLab/protocol/releases.

**Required packages:**

* numpy
* scipy
* pandas>=0.17.0
* PyYAML
* matplotlib (optional, for plotting)
* seaborn>=0.7.0 (optional, for plotting)

If you don't have the above required packages, you will need to install them. For the following commands to work you will need [pip](http://pip.readthedocs.org/en/latest/installing.html). If you are using a system wide python, you will need to use `sudo` before the pip command.

```bash
$ cd protocol
$ pip install -r requirements.txt
```

Once finished installing requirements, you can then use the script `protocol/protocol.py`. You can check the build worked by the following:

```bash
$ python protocol/protocol.py --help
```
