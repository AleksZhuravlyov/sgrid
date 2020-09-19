/* MIT License
 *
 * Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "Sgrid.h"

namespace py = pybind11;
using namespace pybind11::literals;


PYBIND11_MODULE(sgrid_bind, m) {
    py::class_<Sgrid, std::shared_ptr<Sgrid>>(m, "Sgrid")
            .def(py::init<const std::vector<uint16_t>,
                         const std::vector<double>,
                         const std::vector<double>>(),
                 "points_dims"_a, "points_origin"_a, "spacing"_a)

            .def(py::init<std::string>(), "file_name"_a)

            .def("save", &Sgrid::save, "file_name"_a, "isStructured"_a=false)

            .def("calculate_i_cell", &Sgrid::calculateICell, "i_x"_a, "i_y"_a, "i_z"_a)
            .def("calculate_i_x_cell", &Sgrid::calculateIXCell, "i_cell"_a)
            .def("calculate_i_y_cell", &Sgrid::calculateIYCell, "i_cell"_a)
            .def("calculate_i_z_cell", &Sgrid::calculateIZCell, "i_cell"_a)

            .def("calculate_i_x_face", &Sgrid::calculateAxisFace, "i_face"_a)
            .def("calculate_i_x_face", &Sgrid::calculateOffsetFace, "axis"_a)
            .def("calculate_i_face", &Sgrid::calculateIFace,
                 "axis"_a, "i_x"_a, "i_y"_a, "i_z"_a)
            .def("calculate_i_x_face", &Sgrid::calculateIXFace, "i_face"_a)
            .def("calculate_i_y_face", &Sgrid::calculateIYFace, "i_face"_a)
            .def("calculate_i_z_face", &Sgrid::calculateIZFace, "i_face"_a)

            .def("set_cells_type", &Sgrid::setCellsType, "name"_a, "cells"_a)
            .def("process_type_by_cells_type",
                 &Sgrid::processTypesByCellsType, "name"_a)

            .def_readwrite("points_dims", &Sgrid::_pointsDims)
            .def_readwrite("points_N", &Sgrid::_pointsN)
            .def_readwrite("points_origin", &Sgrid::_pointsOrigin)
            .def_readwrite("spacing", &Sgrid::_spacing)

            .def_readwrite("cells_dims", &Sgrid::_cellsDims)
            .def_readwrite("cells_N", &Sgrid::_cellsN)
            .def_readwrite("cell_V", &Sgrid::_cellV)
            .def_property("types_cells",
                          &Sgrid::getTypesCells, &Sgrid::setTypesCells)

            .def_readwrite("faces_dimss", &Sgrid::_facesDimss)
            .def_readwrite("faces_Ns", &Sgrid::_facesNs)
            .def_readwrite("faces_N", &Sgrid::_facesN)
            .def_readwrite("faces_Ss", &Sgrid::_facesSs)
            .def_readwrite("faces_axes", &Sgrid::_facesAxes)
            .def_property("types_faces",
                          &Sgrid::getTypesFaces, &Sgrid::setTypesFaces)

            .def_readwrite("neighbors_faces", &Sgrid::_neighborsFaces)
            .def_readwrite("neighbors_cells", &Sgrid::_neighborsCells)

            .def_readwrite("normals_neighbors_cells", &Sgrid::_normalsNeighborsCells)
            .def_readwrite("normals_neighbors_faces", &Sgrid::_normalsNeighborsFaces)

            .def_property("points_arrays", &Sgrid::getPointsArrays,
                          &Sgrid::setPointsArrays)
            .def_property("cells_arrays", &Sgrid::getCellsArrays,
                          &Sgrid::setCellsArrays)
            .def_property("faces_arrays", &Sgrid::getFacesArrays,
                          &Sgrid::setFacesArrays);

    m.def("save_files_collection_to_file", &saveFilesCollectionToFile,
          "file_name"_a, "files_names"_a, "files_descriptions"_a);
}

