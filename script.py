# MIT License
# 
# Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
import os
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from sgrid import Sgrid

points_dims = np.array([10, 10, 2], dtype=np.int32)
points_origin = np.array([0., 0., 0.], dtype=np.float)
spacing = np.array([1., 1., 1.], dtype=np.float)

sgrid = Sgrid(points_dims, points_origin, spacing)
print("cells_dims", sgrid.cells_dims)
print()

print("points_dims", sgrid.points_dims)
print("spacing", sgrid.spacing)
print("points_origin", sgrid.points_origin)
print("points_N", sgrid.points_N)
print("cells_dims", sgrid.cells_dims)
print("cells_N", sgrid.cells_N)
print("faces_dims", sgrid.faces_dims)
print("faces_N", sgrid.faces_N)
print("cell_V", sgrid.cell_V)
print("face_S", sgrid.face_S)
print("neighbors_faces", sgrid.neighbors_faces)
print("neighbors_cells", sgrid.neighbors_cells)
print("normals_neighbors_cells", sgrid.normals_neighbors_cells)
print("normals_neighbors_faces", sgrid.normals_neighbors_faces)

points_array1 = np.random.rand(sgrid.points_N)
points_array2 = np.random.rand(sgrid.points_N)
points_arrays = {"points_array1": points_array1, "points_array2": points_array2}
sgrid.points_arrays = points_arrays

cells_array1 = np.random.rand(sgrid.cells_N)
cells_array2 = np.random.rand(sgrid.cells_N)
cells_arrays = {"cells_array1 ": cells_array1, "cells_array2": cells_array2}
sgrid.cells_arrays = cells_arrays

print("points_arrays", sgrid.points_arrays)
print("cells_arrays", sgrid.cells_arrays)

sgrid.save('file_name.vtu')

# sgridRead = Sgrid('file_name.vtu')
#
# print()
# print("points_dims", sgridRead.points_dims)
# print("spacing", sgridRead.spacing)
# print("points_origin", sgridRead.points_origin)
# print("points_N", sgridRead.points_N)
# print("cells_dims", sgridRead.cells_dims)
# print("cells_N", sgridRead.cells_N)
# print("cell_V", sgridRead.cell_V)
# print("face_S", sgridRead.face_S)
# print("neighbors_faces", sgridRead.neighbors_faces)
# print("neighbors_cells", sgridRead.neighbors_cells)
# print("normals_neighbors_cells", sgridRead.normals_neighbors_cells)
# print("normals_neighbors_faces", sgridRead.normals_neighbors_faces)
# print("points_arrays", sgridRead.points_arrays)
# print("cells_arrays", sgridRead.cells_arrays)
#
# sgridRead.save('read.vtu')


array = np.array([0, 1, 2], dtype=np.int32)
sgrid.set_cells_type('test', array)
sgrid.process_faces_type_by_cells_type('test')
print('sgrid.types_cells[\'test\']', sgrid.types_cells['test'])
if "test_nonbound" in sgrid.types_cells:
    print('sgrid.types_cells[\'test_nonbound\']',  sgrid.types_cells['test_nonbound'])
print('sgrid.types_cells[\'test_bound\']', sgrid.types_cells['test_bound'])
if "test_nonbound" in sgrid.types_faces:
    print('sgrid.types_faces[\'test_nonbound\']',  sgrid.types_faces['test_nonbound'])
print('sgrid.types_faces[\'test_bound\']', sgrid.types_faces['test_bound'])



