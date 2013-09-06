+------------+
|    TCA     |
+------------+


Version 1.0 - September 25th, 2003
----------------------------------






Description
-----------

The tca package is a  Matlab program that implements the tree-dependent
component analysis (TCA) algorithms that extends the  independent 
component analysis (ICA), where instead of looking for a linear transform
that makes the data components independent, we are looking for components
that can be best fitted in a tree structured graphical model. The TCA model
can be applied in any situation where the data can be assumed to have been
transformed by an unknown linear transformation.

In addition, the TCA algorithm can be specialized to provide a principled
way of finding clusters in ICA, where components in the same cluster are
dependent, but independent from components in other clusters.

The TCA algorithm can be applied to non Gaussian temporally independent 
sources or Gaussian stationary sources.

For more information, please read the following paper:

Francis R. Bach, Michael I. Jordan. Beyond independent components: 
trees and clusters, to appear in Journal of Machine Learning Research, 2003.



The tca package is Copyright (c) 2003 by Francis Bach. If you
have any questions or comments regarding this package, or if you want to
report any bugs, please send me an e-mail to fbach@cs.berkeley.edu. The
current version 1.0 has been released on September, 25th 2003. It has been
tested on matlab 6.  Check regularly the following for
newer versions: http://www.cs.berkeley.edu/~fbach

The package contains ica code (jader.m) by Jean-Francois Cardoso
(cardoso@tsi.enst.fr) and some graphical model utilities from 
Kevin Murphy (murphyk@ai.mit.edu).




Installation
------------

1. Unzip all the .m files in the same directory



2. (Optional) if you want a faster implementation which uses pieces of C
code: at the matlab prompt, in the directory where the package is
installed, type:

 >> mex chol_gauss.c

and

 >> mex chol_hermite.c

It should create compiled files whose extensions depends on the platform
you are using:
      Windows: chol_gauss.dll     and  chol_hermite.dll 
      Solaris: chol_gauss.mexsol  and  chol_hermite.dll
      Linux  : chol_gauss.mexglx  and  chol_hermite.dll

To check if the file was correcly compiled, type

 >> which chol_gauss
 >> which chol_hermite

and the name of the compiled versions should appear. If you have any
problems with the C file of if you are using a platform i did not
mention, please e-mail me.





How to use the tca package
---------------------------------

The functions that you should use to run the TCA algorithm are
'tca', 'cluster_tca', 'tca_statio' and 'cluster_tca_statio' (functions with
default setting of parameters) and 'tca_options' (where 
various options can be tried).

A detailed description of its options are described inside
the file and can be reached by simply typing 'help tca_options' at the
matlab prompt. 3 simple demonstration scripts are provided :
'demo_tca1', 'demo_tca2', 'demo_tca3'.

NB-1: all the data should be given in columns, that is, if you have m
components and N samples, the matrix should be m x N.


NB-2: you need to add the three directories 'optimization', 'contrast'
and 'utilities' to the MATLAB path, without relative paths
(which is done automatically in the demo scripts).









