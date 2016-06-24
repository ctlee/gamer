Generally, there are two types of file for our program: *.tet and *.tet3.

*.tet file contains the basic information for a tetrahdral mesh: the coordinate and the indices of nodes in each tetrahedron.
*.tet3 file contains more information: the face information in the mesh (a face is a triangle which has one or two adjacent tetrahdra), the boundary information in a tet and the relationship between the faces and the tetrahedra.

The index is 0-based.

Description of tet file:
TET
#nodes #tets
x0 y0 z0
x1 y1 z1
...
x_{n-1} y_{n-1} z_{n-1}
4 i00 i01 i02 i03 b00 b01 b02 b03
4 i10 i11 i12 i13 b10 b11 b12 b13
...
4 i_{m-1}0 i_{m-1}1 i_{m-1}2 i_{m-1}3 b_{m-1}0 b_{m-1}1 b_{m-1}2 b_{m-1}3

Here, the (xk,yk,zk) is the coordinate of the k-th node, (ik0,ik1,ik2,ik3) are the indices of the four nodes in the k-th tetrahedron, and (bk0,bk1,bk2,bk3) can be any intgers.

The routine: tet2tet3 can convert a tet file to a tet3 file, which is the file format of the TetSmooth routine.