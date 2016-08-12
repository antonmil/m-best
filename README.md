# Joint Probabilistic Matching Using m-Best Solutions

This code accompanies the paper

Joint Probabilistic Matching Using m-Best Solutions
S. H. Rezatofighi, A. Milan, Z. Zhang,  A. Dick, Q. Shi, I. Reid 

Please cite it if you find the code useful
```
    @inproceedings{Rezatofighi:2016:CVPR,
        Author = {Rezatofighi, S. H. and Milan, A. and Zhang, Z. and Shi, Q. and Dick, A. and Reid, I.},
        Booktitle = {CVPR},
        Title = {Joint Probabilistic Matching Using m-Best Solutions},
        Year = {2016}
    }
```

The package contains all code and data to reproduce the results from the paper.

# Requirements and dependencies
    * Matlab
    * Gurobi solver
    

# Usage

There are three sub-folders, one for each example application. Each one contains a file `AddPath.m` where you need to adjust the path to your Gurobi installation.


## Re-ID
Navigate to `./Re-ID` and run `Demo_ReID_mbst.m`.  See the source code for further instruction. The scripts `Main_CVPR_Results_*` should produce numbers and plots as in the paper.

## Sequential Re-ID
Navigate to `./Sequential Re-ID` and run `SeqReID_Demo.m`.

## Feature Matching
You will first need to compile the mex files for BP Matching. To that end, go to 
`./Feature-Matching/Functions_Codes/SourceCodes` and run compileMex.m. Then
Navigate to `./Feature-Matching` and run `demo*.m` for car or motorbike dataset, respectively.


See the respective source code for further information and instructions. 

# Known Issues

    * The SMCM method is disabled for feature matching. Should you wish to compute the results, you will need to compile the package first.
    
    * The results may not correspond 100% to those reported in the paper due to random number generation, in particular for feature point matching.

# License

BSD License