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

points_dims = [8, 2, 2]
points_origin = [0., 0., 0.]
spacing = [1., 1., 1.]

sgrid = Sgrid(points_dims, points_origin, spacing)

cells_arrays = {'sells': np.arange(sgrid.cells_N, dtype=np.float64)}
faces_arrays = {'faces': np.arange(sgrid.faces_N, dtype=np.float64)}
points_arrays = {'points': np.arange(sgrid.points_N, dtype=np.float64)}

sgrid.cells_arrays = cells_arrays
sgrid.faces_arrays = faces_arrays
sgrid.points_arrays = points_arrays

print('points_dims', sgrid.points_dims)
print('points_N', sgrid.points_N)
print('points_origin', sgrid.points_origin)
print('spacing', sgrid.spacing)
print()
print('cells_dims', sgrid.cells_dims)
print('cells_N', sgrid.cells_N)
print('cell_V', sgrid.cell_V)
print('types_cells', sgrid.types_cells)
print()
print('faces_dimss', sgrid.faces_dimss)
print('faces_Ns', sgrid.faces_Ns)
print('faces_N', sgrid.faces_N)
print('faces_Ss', sgrid.faces_Ss)
print('faces_axes', sgrid.faces_axes)
print('types_faces', sgrid.types_faces)
print()
print('neighbors_faces', sgrid.neighbors_faces)
print('neighbors_cells', sgrid.neighbors_cells)
print()
print('normals_neighbors_cells', sgrid.normals_neighbors_cells)
print('normals_neighbors_faces', sgrid.normals_neighbors_faces)
print()
print('points_arrays', sgrid.points_arrays)
print('cells_arrays', sgrid.cells_arrays)
print('faces_arrays', sgrid.faces_arrays)
print()
active_cells = np.arange(sgrid.cells_N, dtype=np.uint64)
sgrid.set_cells_type('active', active_cells)
sgrid.process_type_by_cells_type('active')
print('cells:')
print('nonbound', sgrid.types_cells['nonbound'])
print('active', sgrid.types_cells['active'])
print('active_bound', sgrid.types_cells['active_bound'])
print('active_nonbound', sgrid.types_cells['active_nonbound'])
print('faces:')
print('active', sgrid.types_faces['active'])
print('active_bound', sgrid.types_faces['active_bound'])
print('active_nonbound', sgrid.types_faces['active_nonbound'])
print()

sgrid.save_cells('cells.vtu')
sgrid.save_faces('faces.vtu')
