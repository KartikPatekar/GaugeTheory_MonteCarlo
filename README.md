# GaugeTheory_MonteCarlo

**Context**
The codes in this repository were used to simulate the dynamics of Gauge fields as described in the [1]. THe simulations were performed for 4 fields (flavors of Higg's field) in 3 dimensions.

[1]. H. D. Scammell, K. Patekar, M. S. Scheurer, and S. Sachdev, “Phases of 𝑆𝑈(2) gauge theory with multiple adjoint Higgs fields in 2+1 dimensions.”, Phys. Rev. B 101, 205124, May 2020



**Code description**
HiggsCompWiseUpdate.cpp: Updates one component of one field at a time.
HiggsFieldWiseUpdate.cpp: Updates entire field in one go.
HiggsParallel.cpp: Parallelize the computations for faster simulations.
