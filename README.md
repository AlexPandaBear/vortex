# PIR VORTEX METHOD
Alexandre DUTKA - ISAE-SUPAERO - 05/2020

## Requirements
CMake  
C++ compiler  
Python 3+  
PyBind11  
NumPy  
Matplotlib  

## How to install
```console
user@linux:~/path/to/code$ mkdir build  
user@linux:~/path/to/code$ cd build  
user@linux:~/path/to/code/build$ cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD/install    
user@linux:~/path/to/code/build$ make install  
```

## How to use
```console
user@linux:~/path/to/code$ mkdir data
user@linux:~/path/to/code$ cd python
user@linux:~/path/to/code/python$ export PYTHONPATH=$PWD/../build/install/lib/python
user@linux:~/path/to/code/python$ python -i <script>
```
