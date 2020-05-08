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
mkdir build  
cd build  
cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD/install  
make  
make install  

## How to use
cd python  
export PYTHONPATH=$PWD/../build/install/lib/python  
python -i <script>  

```console
foo@bar:~$ whoami
foo
```
