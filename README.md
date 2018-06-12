# Pair-PES

This jupyter notebook project uses modern electronic structure methods to build a potential energy surface (PES) of molecular dimers. PES scans can be performed along many characteristic coordinates and combinations of coordinates using a variety of quantum mechanical theory, basis sets, and programs. 

## Getting Started

These instructions will get your copy of the project up and running on your local machine. 

### Prerequisites

In order to run the program you'll obviously need jupyter notebook, but also a few other modules in order to initialize the correct notebook environment. These include

```
numpy
matplotlib
parsl
ipywidgets
cirpy
ipykernel
electronic structure package, like NWChem or Molpro
```

### Installing

In order to install jupyter notebook and the requisite modules, it is advised to use the conda source and package management system. This makes for simple installation and ensures the correct dependencies. 

1. If you don't already have anaconda installed, download the latest version from the [website](https://www.anaconda.com/download/).  

2. On the terminal, create a new conda environment with at least python 3, since certain modules require more than python 2.7. Also when creating this new environment you can specify the packages to be installed all in one step. Refer to the conda documentation for more details on [managing environments](https://conda.io/docs/user-guide/tasks/manage-environments.html). 

    `conda create --name myenv python=3.6 numpy matplotlib parsl ipywidgets ipykernel cirpy nwchem`
    
Here, a new conda environment, named `myenv` was created with python3.6 and all the necessary modules to run the Pair-PES notebook. 

3. Next, you should activate the new environment and ensure that jupyter will run the correct version of python. To do this, execute the following on the command line.

    `source activate myenv`
    `python -m ipykernel install --name myenv --display-name "Python (myenv)"`
    
The first line activates the created environment. And the second line installs the kernel under the name `Python (myenv)`.

4. You should now be able to launch the jupyter notebook by typing `jupyter notebook` in the command line. Once you are at the jupyter web browser home tab, select the file `pes-2b.ipynb`, which will take the main notebook. Make sure that the kernel at the far right, is the one you just created (`Python (myenv)`). If it is not, you can change it in the Kernel tab at the top.

5. There are instructions and directions throughout the Pair-PES notebook about its functionality; however, if you'd like to simply see the data output, you can simply run the notebook all the way through. To do this, at the top of the screen, select the Kernel tab, and then select 'Restart & Run All'. The notebook will take a few minutes to run, but you should see the following figure once the notebook is complete. 

Enjoy!

## Authors

* **Olaseni Sode** - *Initial work* - [SodeLab](https://github.com/sodelab)

## Acknowledgments

* Thanks Yadu Babuji
