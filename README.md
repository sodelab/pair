# Pair-PES

This jupyter notebook project uses modern electronic structure methods to build a potential energy surface (PES) of molecular dimers. PES scans can be performed along many characteristic coordinates and combinations of coordinates using a variety of quantum mechanical theory, basis sets, and programs. 

## Getting Started

These instructions will get your copy of the project up and running on your local machine. 

### Prerequisites

In order to run the program you'll obviously need jupyter notebook, but also a few other modules in order to initialize the correct notebook environment. These include:

```
numpy
matplotlib
parsl
ipywidgets
cirpy
ipykernel
electronic structure package, like NWChem or Molpro
```

If all of the following are already installed on your machine, there is no need to proceed through the `Installing` section of the README.

### Installing

In order to install jupyter notebook and the requisite modules, it is advised to use the conda source and package management system. This makes for simple installation and ensures the correct dependencies. 

1. Make sure to download the latest version of the Pair-PES repository, by clicking the Clone or Download link on the main page of this [repository](https://github.com/sodelab/pair). 

2. If you don't already have anaconda installed, download the latest version from the [website](https://www.anaconda.com/download/).  

3. On the terminal, create a new conda environment with at least python 3, since certain modules require more than python 2.7. Also when creating this new environment you may specify the packages to be installed all in one step. Refer to the conda documentation for more details on [managing environments](https://conda.io/docs/user-guide/tasks/manage-environments.html). Here, a new conda environment, named `myenv` was created with python3.6 and all the necessary modules to run the Pair-PES notebook. 

    `conda create --name myenv python=3.6 numpy matplotlib ipywidgets ipykernel`
    
4. In order to install the cirpy and nwchem packages, type the following command into the command line:

    `conda install --name myenv -c mcs07 cirpy`
    
    `conda install --name myenv -c insilichem nwchem`
   
5. Next, activate the new environment and install the remaining necessary packages (parsl). Enter the following commands:

    `source activate myenv`
    
    `pip install parsl`

6. To ensure that jupyter runs the correct version of python and in the appropriate environment, execute the following on the command line. This line installs an ipython kernel under the name `Python (3.6)`. 

    `python -m ipykernel install --name myenv --display-name "Python (3.6)"`
    
7. You should now be able to launch the jupyter notebook by typing `jupyter notebook` in the command line. Once you are at the jupyter web browser home tab, navigate to the directory containing the `pes-2b-ipynb` file and select it. This will take the main notebook. Make sure that the kernel at the far right, is the one you just created (`Python (3.6)`). If it is not, you can change it in the Kernel tab at the top.

8. There are instructions and directions throughout the Pair-PES notebook about its functionality; however, if you'd like to simply see the data output, you can run the notebook all the way through. To do this, at the top of the screen, select the Kernel tab, and then select 'Restart & Run All'. The notebook will take a few minutes to run, but you should see the following figure once it is complete. 

Enjoy!

## Authors

* **Olaseni Sode** - *Initial work* - [oosode](https://github.com/oosode)

## Acknowledgments

* Thanks Yadu Babuji
