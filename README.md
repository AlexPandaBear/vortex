# VORTEX METHOD CODE
Alexandre DUTKA - ISAE-SUPAERO - 05/2020

![alt text](logo.png "A quick preview...")

## Requirements
To compile the library :
- CMake  
- A C++ compiler  
- The PyBind11 library  

To use it with Python :
- Python 3+  
- The Numpy and Matplotlib libraries  

To generate the reference documentation :
- Doxygen  

## How to install
Go to the code directory and execute the following commands to create a build directory to compile the library and a data directory to store simulation data.
```console
user@linux:~/path/to/code$ mkdir build  
user@linux:~/path/to/code$ mkdir data  
```

Execute the following commands to compile the C++ library.
```console
user@linux:~/path/to/code$ cd build  
user@linux:~/path/to/code/build$ cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD/install    
user@linux:~/path/to/code/build$ make install  
```

Execute the following command to create the documentation of the library, which will then be available by launching index.html in the build/doc directory.
```console
user@linux:~/path/to/code/build$ make reference_doc
```

## How to use
Execute these commands to get to the python directory and add the library to the PYTHONPATH.
```console
user@linux:~/path/to/code$ cd python
user@linux:~/path/to/code/python$ export PYTHONPATH=$PWD/../build/install/lib/python
```

Use this command to execute the python scripts.
```console
user@linux:~/path/to/code/python$ python <script>
```
