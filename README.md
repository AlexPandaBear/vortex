# PIR VORTEX METHOD
Alexandre DUTKA - ISAE-SUPAERO - 05/2020

## Requirements
G++  
Python3+  
CMake  
PyBind11  
NumPy  
Matplotlib  

## How to install
```console
user@linux:~$ mkdir build  
user@linux:~$ cd build  
user@linux:~$ cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD/install  
user@linux:~$ make  
user@linux:~$ make install  
```

## How to use
```console
user@linux:~$ cd python
user@linux:~$ export PYTHONPATH=$PWD/../build/install/lib/python
user@linux:~$ python -i <script>
```
