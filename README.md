# Harmonic-Surface-Mapping-Algorithm-on-CPU-and-GPU
Harmonic Surface Mapping Algorithm, described in paper [Harmonic surface mapping algorithm for fast electrostatic sums](https://arxiv.org/abs/1806.04801) is an efficient implementation for electrostatic pairwise sums of an infinite number of images accelerated by Fast Multiple method(FMM) and graphics processing units(GPU). Numerical calculations of the Madelung constant, electrostatic energy of ions in a metallic cavity, and the time performance for large-scale systems show that the HSMA is accurate and fast, and thus is attractive for many applications.

Main file `GPUmain.cu` is in CUDA programming language and is available on NVIDIA Tesla P100 and other suitable graphics cards.

Another main File `FMMmain.cpp` is in C++ programming language and accelerated by OpenMp. Since FMM code is not an open source, we only present the kernel code of this method in this CPU version. We leave the FMM program interface here and one can contact Pro. Greengard to obtain the FMM code. 

Other header files such like `InitialSet.h` and `BulidTestPoint.h` are designed for both two main files. 

We hope this code may be helpful.
```
                           Jiuyang Liang
                           Ph.D candidate
                           School of Mathematical Science and Institute of Natural Science
                           Shang Hai Jiao Tong University
                           [Homepage in Github](https://github.com/LiangJiuyang/)
```
[HomePage in Github](https://github.com/LiangJiuyang)
