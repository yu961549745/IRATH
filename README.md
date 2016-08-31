# IRATH
Maple implement of improved tanh method to find exact wave solutions of nonlinear evolution equations automatically.

This project is a improved platform of the `RAHT`[1] method. You can see its description in `paper.pdf`(In Chinese).

## Abstract
We use the improved tanh function method to find exact wave solutions of nonlinear evolution equations automatically, and implement a Maple package for this method.We improved the classic tanh function method in three aspects: (1) We completely considered the of order balance condition. (2) We proved the correctness of the order transform method which deal with the non-positive integer orders, and proposed the service condition of this method. (3) We proposed a more effective and fast algorithm to solve the polynomial equations. We solved a lot of nonlinear evolution equations utilizing the package IRATH. It successfully recovered all previously known solutions. And more, for some equations, we have found new solutions and more general form of solutions.

## 摘要
本文针对非线性演化方程的精确波解自动化求解问题，采用了改进的双曲正切函数方法进行求解，在Maple上实现了自动求解的软件包IRATH。本文从三个方面对传统双曲正切函数方法进行了改进：(1)对阶数平衡的条件进行了完整的分析；(2)在处理非正整数阶数时，对前人的阶数变换方法进行了严格论证，并讨论了可以使用阶数变换的条件；(3)对于多项式方程组的求解问题，本文提出了求解效率更高、效果更好的实现方法。利用IRATH本文已成功求解大量的非线性演化方程，不仅求出了已知的解，对于某些方程还得到了新的解以及形式更一般的解。

## Files & Usage
+ IRATH.mpl, is the source code of this method.
+ IRATH.mpl, is the library file of this method. You can put is to the system library path (`libname`), and use `with(IRATH)` to call the `findTanhSolutions` method.
+ cmp.mw, is the compare demo of between `RATH` and `IRATH`.
+ make.mw, can make library file from IRATH.mpl.
+ paper.pdf, is the description of this method. 

## Refered Fils
The `maple.ind` and `maple.lib` is the package files of `charsets`, refered to [Dongming Wang](http://www-polsys.lip6.fr/~wang/epsilon/index.html).

The `rath.mla` is the package file of `RATH` method, refered to the [this paper](http://www.sciencedirect.com/science/article/pii/S0010465502005593).

## Reference
[1] Li Z, Liu Y. RATH: A Maple package for finding travelling solitary wave solutions to nonlinear evolution equations[J]. Computer Physics Communications, 2002, 148(2): 256-266.