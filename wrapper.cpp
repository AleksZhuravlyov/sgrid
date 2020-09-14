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
            .def(py::init<Eigen::Ref<Eigen::Vector3i>,
                         Eigen::Ref<Eigen::Vector3d>,
                         Eigen::Ref<Eigen::Vector3d>>(),
                 "points_dims"_a,
                 "points_origin"_a,
                 "spacing"_a)
            .def(py::init<std::string>(), "file_name"_a)

            .def("save", &Sgrid::save, "file_name"_a)

            .def("calculate_i_cell", &Sgrid::calculateICell, "i_x"_a, "i_y"_a, "i_z"_a)
            .def("calculate_i_x_cell", &Sgrid::calculateIXCell, "i_cell"_a)
            .def("calculate_i_y_cell", &Sgrid::calculateIYCell, "i_cell"_a)
            .def("calculate_i_z_cell", &Sgrid::calculateIZCell, "i_cell"_a)


            .def_property("points_dims",
                          &Sgrid::getPointsDims, &Sgrid::setPointsDims)
            .def_property("spacing",
                          &Sgrid::getSpacing, &Sgrid::setSpacing)
            .def_property("points_origin",
                          &Sgrid::getPointsOrigin, &Sgrid::setPointsOrigin)

            .def_readwrite("points_N", &Sgrid::_pointsN)

            .def_property("cells_dims",
                          &Sgrid::getCellsDims, &Sgrid::setCellsDims)
            .def_readwrite("cells_N", &Sgrid::_cellsN)
            .def_property("types_cells",
                          &Sgrid::getTypesCells, &Sgrid::setTypesCells)

            .def_property("faces_dims",
                          &Sgrid::getFacesDims, &Sgrid::setFacesDims)
            .def_readwrite("faces_N", &Sgrid::_facesN)
            .def_readwrite("faces_axis", &Sgrid::_facesAxis)
            .def_property("types_faces",
                          &Sgrid::getTypesFaces, &Sgrid::setTypesFaces)

            .def_readwrite("cell_V", &Sgrid::_cellV)
            .def_property("face_S", &Sgrid::getFaceS, &Sgrid::setFaceS)

            .def_property("neighbors_faces",
                          &Sgrid::getNeighborsFaces, &Sgrid::setNeighborsFaces)
            .def_property("neighbors_cells",
                          &Sgrid::getNeighborsCells, &Sgrid::setNeighborsCells)

            .def_property("normals_neighbors_cells",
                          &Sgrid::getNormalsNeighborsCells,
                          &Sgrid::setNormalsNeighborsCells)
            .def_property("normals_neighbors_faces",
                          &Sgrid::getNormalsNeighborsFaces,
                          &Sgrid::setNormalsNeighborsFaces)

            .def_property("points_arrays",
                          &Sgrid::getPointsArrays, &Sgrid::setPointsArrays)
            .def_property("cells_arrays",
                          &Sgrid::getCellsArrays, &Sgrid::setCellsArrays)
            .def_property("faces_arrays",
                          &Sgrid::getFacesArrays, &Sgrid::setFacesArrays);

    m.def("save_files_collection_to_file", &saveFilesCollectionToFile,
          "file_name"_a, "files_names"_a, "files_descriptions"_a);
}

