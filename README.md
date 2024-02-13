# TAPERSS-Private
Physics-based prediction of 3D RNA structure by fragment assembly.

For usage please see ./docs/TAPERSS.txt.

GSL must be installed for usage of TAPERSS.
On debian based systems you may install with:
```
sudo apt-get install libgsl-dev
```

Ensure you have C++ compiler which allows at least C++17 features. Examples: Clang version 8, or GCC version 8.3
The compiler used maybe configured by editing the "CC =" option in the Makefile.

To compile run "make".
