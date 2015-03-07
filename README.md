# PhD-ANMInternalCoordinates

This is an implementation of an anisotropic network model (ANM) for proteins performed in internal coordinates (IC) for PELE++*. The code is based on the formulae shown in the [PhD thesis](http://tdx.cat/handle/10803/81963) of José Ramón López Blanco (JRLB). It has been tuned using the code of iNMA v1.2, developed by JRLB himself,  and  the [imods](http://imods.chaconlab.org/) server from Pablo Chacon group.

The goal of this project was to implement the ANM IC analysis in PELE++, but also to do it in an educational manner, as it was planned as a Bachelor Thesis. Because of this, functions have been implemented the more similar we could to their formulae counterparts. Obviously, this makes our code to be really inneficient compared with JRLB's implementation. I will try to increse its performance in a close future without sacrificing its readability.

Our code was started in late 2013 by me and Alba Rincón. This repository contains a slightly improved version of the one presented  as her Computer Engineering Bachelor Thesis ([published here](http://hdl.handle.net/2099.1/24417)) at the Barcelona Faculty of Computer Sciences (FIB, UPC). 

This collection of classes were developed as part of the PELE++ project in a separate repository. As the last is currently not open source, it has been necessary to create a new repository with the pieces of code that:
- Are part of the ANM IC code.
- Do not contain sensible information about PELE++, but can help to understand the ANM IC code.
- Implement changes/improvements from the base PELE++ code.  

The ANM IC code rights go to their authors, the Barcelona Supercomputing Center (BSC-CNS) and Universitat Politécnica de Barcelona (UPC). Any piece of code containing functions that were not developed in the context of the ANM IC and are part of the PELE++ project are under the PELE++ copyright license (included as 'COPYRIGHT' in the root folder) and may not be used without express consent of their owner/s.  

# Authors
Main authors: Alba Ricón Muñoz and [Víctor A. Gil Sepúlveda](http://victor-gil-sepulveda.github.io/)
Also contributed: [Manuel Rivero](http://garajeando.blogspot.com.es/), Jorge Estrada and Pedro Riera

# Special Thanks
We would like to thank Jose Ramón López Blanco for his help and endless patience (and of course, for the code he led us!). It would have been impossible to do it without him.

\* The C++ rewrite of PELE (J. Comp. Chem., 31:1224-35 (2010))
