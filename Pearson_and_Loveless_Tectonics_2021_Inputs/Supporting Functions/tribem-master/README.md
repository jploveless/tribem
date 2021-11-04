# tribem
### Boundary element code using triangular dislocation elements in Matlab

Based on algorithms published in:

Meade, B. J. (2007), Algorithms for the calculation of exact displacements, strains, and stresses for triangular dislocation elements in a uniform elastic half space, *Computers and Geosciences*, **33**, 1064–1075, [doi:10.1016/j.cageo.2006.12.003](http://dx.doi.org/10.1016/j.cageo.2006.12.003).

__tribem__ uses [__tridisl__](https://github.com/jploveless/tridisl) as a submodule. To clone __tribem__ using `git` on the command line, run

    $ git clone --recursive https://github.com/jploveless/tribem.git

If not using `git` on the command line, use the following workflow to access the necessary subfunctions: 
- Download __tribem__ as a `.zip` file and unzip
- Download [__tridisl__](https://github.com/jploveless/tridisl) as a `.zip` file and unzip
- Move the contents of the `tridisl` folder to empty `tribem/tridisl` folder
- Add the `tribem` and `tribem/tridisl` directories to your MATLAB path
