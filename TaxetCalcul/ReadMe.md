<!-- PROJECT LOGO -->
![alt text](BioDesignLogo.png)

# ReadMe for TaxetCalcul folder

This folder contains the code employed within the research paper: 

***Size does not matter: allometry reveals consistent pressure and sliding velocity in quadrupedal mammalian elbow***

*Kalenia Marquez-Florez, Loïc Tadrist, Santiago Arroyave-Tobon, Jean-Marc Linares*

Aix Marseille Univ, CNRS, ISM, Marseille, France

[![LinkedIn][linkedin-shield]][linkedin-url]


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li><a href="#prerequisites">Prerequisites</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project
The codes here are part of the results from the research paper: *Size does not matter: allometry reveals consistent pressure and sliding velocity in quadrupedal mammalian elbow*. This repository contains the Python and R scripts used within the project.

<!--
*** Here's a blank template to get started: To avoid retyping too much info. Do a search and replace with your text editor for the following: `github_username`, `repo_name`, `twitter_handle`, `linkedin_username`, `email_client`, `email`, `project_title`, `project_description`
-->

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With

* [![Python][Python.png]][Python-url]
* [![R][R.png]][R-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

Assure that the following packages from <a href="#prerequisites">Prerequisites</a> of Python and R are installed before running the code. Also, be sure to have track of the location of the [".vtk" files](DATA/For_publishing/Distal) of the distal surface of the humeri, as well as the files [New_tax.csv](DATA/For_publishing/New_tax.csv) and [3DfilesDb.csv](DATA/For_publishing/3DfilesDb.csv), as they are the inputs for the Python script.

### Prerequisites
#### *For Python*
We developed the Python script for the version 3.6.9.

* json
  ```sh
  pip3 install json
  ```
* NumPy
  ```sh
  pip3 install numpy
  ```
* cmath
  ```sh
  pip3 install cmath
  ```
* csv
  ```sh
  pip3 install csv
  ```
* Pandas
  ```sh
  pip3 install pandas
  ```
  Be sure to previously install NumPy.

* mlxtend
  ```sh
  pip3 install mlxtend
  ```
  Be sure to previously install the external dependencies for this package, including NumPy, Pandas, Matplotlib and Scikit-learn. 

* cylinder_fitting
  ```sh
  pip3 install cylinder_fitting
  ```
  Be sure to previously install the external dependencies for this package: Numpy and SciPy.
* vtk
  ```sh
  pip3 install vtk
  ```
  Be sure to previously install the external dependencies for this package, including OpenGL, CMake, and Python development headers. For more information on installing VTK, please refer to the [VTK installation instructions](https://gitlab.kitware.com/vtk/vtk/-/blob/master/Documentation/dev/build.md).

* Matplotlib
  ```sh
  pip3 install matplotlib
  ```

  Be sure to previously install the external dependencies for this package, including NumPy, libpng, freetype, and Python development header. For more information on installing Matplotlib, please refer to the [Matplotlib installation instructions](https://matplotlib.org/stable/users/installing/index.html).

* SciPy
  ```sh
  pip3 install scipy
  ```
  Be sure to previously install the external dependencies for this package, including NumPy, BLAS, LAPACK, and CMake. For more information on installing SciPy, please refer to the [SciPy installation instructions](https://scipy.org/install/).

* ete3
  ```sh
  pip3 install ete3
  ```
  Be sure to previously install the external dependencies for this package, including NumPy, SciPy, and Matplotlib. For more information on installing ETE Toolkit, please refer to the [ETE Toolkit installation instructions](http://etetoolkit.org/download/).

#### *For R*
We developed the R script for the version 4.2.2 Patched (2022-11-10 r83330).

* ape
* caper
* dplyr
* mosaic
* phytools
* Polychrome
* reticulate
* rlang
* smatr
* viridis

To install these packages and their dependencies, you can use the following code in R:
  ```sh
  install.packages(c("ape", "caper", "dplyr", "mosaic", "phytools", "Polychrome", "reticulate", "rlang", "smatr", "viridis"))
  ```

<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Kalenia Marquez-Florez - kalenia-maria.marquez-florez@univ-amu.fr

Santigo Arroyave-Tobon - santiago.ARROYAVE-TOBON@univ-amu.fr

Loïc Tadrist - loic.TADRIST@univ-amu.fr

Jean-Marc Linares - jean-marc.linares@univ-amu.fr

Project Link: [https://github.com/sarroyavet/BioDesign_joint_morphogenesis](https://github.com/sarroyavet/BioDesign_joint_morphogenesis)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

This research was supported by French Research National Agency (Agence Nationale de la Recherche, ANR) Grant No. ANR-20-CE10-0008, through the project BioDesign.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/kalenia-márquez-flórez-5b686064/
[product-screenshot]: images/screenshot.png

[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 



[Python.png]: https://img.shields.io/badge/Python-0769AD?labelColor=blue?style=plastic&logo=python&logoColor=white
[Python-url]: https://www.python.org/
[R.png]: https://img.shields.io/badge/R-0769AD?style=plastic&logo=R&logoColor=white
[R-url]: https://www.r-project.org/
