# Reluctance motor - Finite element project

This project was done in the course `LEPL1110 Introduction to finite element methods`. It simulates in real-time the rotation of a switched reluctance motor by solving a PDE with FEM. The equations and the method are described in the file `Ressources/motor_consignes.pdf`. The linear system is solved using a preconditioned conjugate gradient method (PCG).

## Installation and execution
```
mkdir build
cd build
cmake ..
make
./myFem
```

> Once the code is executed, the list of commands is available by pressing `H` in the GUI.
