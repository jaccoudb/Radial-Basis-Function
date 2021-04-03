Radial Basis Functions (RBF) are frequently used when it is needed to construct an approximation of a certain multivariable function by interpolation between the existing data. This section will give a very brief description of the method to make the reader familiar with main principles and to show how it can be effectively combined with previously described POD theory in the structural context here of interest. 

For more detailed description on RBF:

`Buhmann, M.D.: Radial Basis Functions. Cambridge University Press, Cambridge (2003)`

Let us assume that we wish to approximate some function $f(\mathbf{x})$, where $\mathbf{x}$ is an $M$-dimensional vector, for which the only information that we have is a set of, say N values xi, called “nodes”, for which the values of the function are known. In classical, local methods of interpolation, the problem is solved locally in the neighborhood of the point x for which the value of the function is required. In this case it doesn’t exist just one continuous function defined over the whole domain where the data are situated. Instead, for any arbitrary value of $\mathbf{x}$ the interpolation is performed involving just a few nearby data.

The approach of RBF is different, and it seeks for one continuous function defined over the whole domain and that depends on the entire data set defined by the pairs of given $N$ nodes and their values of the function. Therefore, the approximation of the function is written as a linear combination of some functions $g_i$, that in general case can be non-linear functions, namely

\[f(\mathbf{x})\approx\sum^{N}_{i=1}\alpha_{i}\cdot{g_i(\mathbf{x})}\]

where $\alpha_{i}$ are coefficients of this combination.