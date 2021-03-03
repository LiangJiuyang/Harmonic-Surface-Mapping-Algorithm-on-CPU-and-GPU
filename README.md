# Harmonic-Surface-Mapping-Algorithm-on-CPU-and-GPU

Note : Recently, we developed an optimized software package of HSMA, including the original HSMA (3D periodic) and our further work of 2D periodic system with dielectric interface (and planar interface). This package has been written in Lammps, or used as an independent package. We plan to publish this package as open source before long. If you want to use this package immediately, please send email to liangjiuyang@sjtu.edu.cn. Belows are old description of original HSMA.

Harmonic Surface Mapping Algorithm, described in paper [Harmonic surface mapping algorithm for fast electrostatic sums](https://arxiv.org/abs/1806.04801) is an efficient implementation for electrostatic pairwise sums of an infinite number of images accelerated by Fast Multiple method(FMM) and graphics processing units(GPU). Numerical calculations of the Madelung constant, electrostatic energy of ions in a metallic cavity, and the time performance for large-scale systems show that the HSMA is accurate and fast, and thus is attractive for many applications.

Main file `GPUmain.cu` is in CUDA programming language and is available on NVIDIA Tesla P100 and other suitable graphics cards.

Another main File `FMMmain.cpp` is in C++ programming language and accelerated by OpenMp. For FMM code, user can download FMM3D package from the website [https://cs.nyu.edu/~harper/kifmm3d/documentation/fmm3d/html/](https://cs.nyu.edu/~harper/kifmm3d/documentation/fmm3d/html/).

Other header files such like `InitialSet.h` and `BulidTestPoint.h` are designed for both two main files. 

The authors thank Mr.Joey Wang from the NVIDIA for the guide on the GPU implementation. 

We hope this code may be helpful.
```
                           Jiuyang Liang
                           Ph.D candidate
                           School of Mathematical Science and Institute of Natural Science
                           Shang Hai Jiao Tong University
                           [Homepage in Github](https://github.com/LiangJiuyang/)
```
[HomePage in Github](https://github.com/LiangJiuyang)
