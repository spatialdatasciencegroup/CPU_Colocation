# CPU Colocation

This is the CPU sequential implementation of Colocation Mining Algorithms in [1], [2], and [3].
### How to Compile

1. Install boost. (https://www.boost.org/users/download/) 
2. Change your current directory to be the directory containing main.cpp  and colocationFinder.cpp
3. In the console, compile with: **g++ -std=c++11 -I (boost installation path) main.cpp and colocationFinder.cpp -o CPU-Colocation

### How to Run

1. We have a data folder for you to test out.
2. To run the compiled CPU-Colocation program: ./CPU-Colocation /home/TestCase/config.txt

### Input Data Description
1. Please see the Technical Document pdf file. 

### References

[1] Arpan Man Sainju, Danial Aghajarian, Zhe Jiang, & Sushil K Prasad, (2018). Parallel grid-based colocation mining algorithms on GPUs for big spatial event data. IEEE Transactions on Big Data.

[2] Arpan Ma Sainju, and Zhe Jiang. "Grid-based colocation mining algorithms on gpu for big spatial event data: A summary of results." International Symposium on Spatial and Temporal Databases. Springer, Cham, 2017.

[3] Huang Y, Shekhar S, Xiong H. Discovering colocation patterns from spatial data sets: a general approach. IEEE Transactions on Knowledge and data engineering. 2004 Nov 1;16(12):1472-85.
