# PLSS
**Projected Linear Systems Solver**

Matlab and Python implementations from

"PLSS: Projected Linear Systems Solver with Randomization", J.J. Brust and M.A. Saunders (2022)

Content:
  * MATLAB/
    * ALGS/
    * AUXILIARY/
    * EXAMPLES/
    * EXPERIMENTS/
    * external/
  * PYTHON/  

Notes: The Matlab codes include all external libraries and methods in order reproduce the article experiments.
That is, functions from `[1]`,`[2]` and `[3]` are used either to generate problems or to compare. The Python
implementations use numpy https://numpy.org/
    
## Example(s)

### Matlab
You can run a Matlab example from within MATLAB/EXAMPLES/

```
>> example_1

********* Linear Systems Solver ******************** 
*         Alg: PLSS                                  
*         Size: m = 400, n = 400                       
*         tol = 0.00000100                                
*         Maxit = 400                                 
**************************************************** 

Iter 	 norm(Res) 	 norm(pk)   
0 	 9.3808e+00 	 3.9836e+00  
1 	 8.5457e+00 	 2.7639e+00 
2 	 1.0515e+01 	 2.6170e+00 
3 	 9.1554e+00 	 2.2788e+00 
4 	 9.0323e+00 	 2.3500e+00 
5 	 9.1176e+00 	 2.2150e+00 
6 	 8.4172e+00 	 2.1014e+00 
7 	 8.4112e+00 	 2.1293e+00 
8 	 8.3266e+00 	 2.0363e+00 
9 	 7.8896e+00 	 1.9875e+00 
10 	 7.8228e+00 	 2.0224e+00 
20 	 6.4755e+00 	 2.3890e+00 
40 	 6.7713e+00 	 1.5057e+00 
60 	 1.1869e+01 	 3.2345e+00 
73 	 1.8782e-07 	 2.7772e-08 

********* Summary ********************************** 
*         Conv: 1                                   
*         Time (s): 0.02                            
*         norm(Res) = 0.0000                          
*         Iter = 74                                  
*         Num. mult. = 148                            
**************************************************** 
```

### Python 3.9
Likewise, you can run a different example in Python from within PYTHON/
and at the console

```
In[85]: runfile('example_3.py')

********* Linear Systems Solver ********************  
*         Alg: PLSS (Python 3.9)                      
*         Size: m = 10, n = 10                        
*         tol = 0.00000100                                 
*         Maxit = 500                                  
****************************************************  

Iter 	 norm(Res) 	 norm(pk)    
0 	 1.4468e+02 	 1.2560e+01   
1 	 9.0747e+01 	 8.7137e+00  
2 	 5.9964e+01 	 6.8153e+00  
3 	 4.1940e+01 	 5.5617e+00  
4 	 2.9450e+01 	 4.7588e+00  
5 	 2.0575e+01 	 4.0649e+00  
6 	 1.3862e+01 	 3.5855e+00  
7 	 8.8485e+00 	 3.1817e+00  
8 	 4.7557e+00 	 2.7828e+00  
9 	 1.5912e+00 	 2.0086e+00  
10 	 8.3808e-05 	 5.2023e-06  

********* Summary **********************************  
*         Conv: 1                                    
*         Time (s): 0.00                             
*         norm(Res) = 0.0000                           
*         Iter = 12                                   
*         Num. mult. = 25                             
**************************************************** 
```

## Cite
You can cite this work as (bibtex)

```
@TechReport{plss22,
  author      = {Johannes J. Brust and Michael A. Saunders},
  title       = {PLSS: Projected Linear Systems Solver with Randomization},
  institution = {Mathematics Department, University of California, San Diego, CA},
  type        = {Technical Report},
  year        = {2022},
  url         = {TBD}
}
```

## Reference codes
[1] T. A. Davis, Y. Hu, and S. Kolodziej, Suitesparse matrix collection. https://sparse.tamu.edu/, 2015–present.

[2] Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011. Software available at http://www.csie.ntu.edu.tw/~cjlin/libsvm

[3] R. M. Gower and P. Richtarik, Randomized iterative methods for linear systems, SIAM
J. Matrix Anal. Appl., 36 (2015), pp. 1660–1690
