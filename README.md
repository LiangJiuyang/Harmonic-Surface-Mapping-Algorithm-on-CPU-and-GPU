# Harmonic-Surface-Mapping-Algorithm-on-GPU
Harmonic Surface Mapping Algorithm, described in paper [Harmonic surface mapping algorithm for fast electrostatic sums](http://xueshu.baidu.com/swd=paperuri%3A%28b6033f4f4dc26ba2ec88d6860eb5d6d1%29&filter=sc_long_sign&tn=SE_xueshusource_2kduw22v&sc_vurl=http%3A%2F%2Farxiv.org%2Fpdf%2F1806.04801.pdf&ie=utf-8&sc_us=12003589625995823630) is an efficient implementation for electrostatic pairwise sums of an infinite number of images by using Fast Multiple method(FMM) or graphics processing units (GPU). 

Main file 'GPUmain.cu' is in CUDA programming language and is available on NVIDIA Tesla P100 and other suitable graphics cards.

Another main File 'FMMmain.cpp'      we only present the kernel code of this method in this CPU version. 

We hope this code may be helpful.
