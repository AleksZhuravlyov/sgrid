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


#include "Sgrid.h"

#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>


#include <pugixml.hpp>

void Sgrid::save(const std::string &fileName, const std::string &type) {

    auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    auto xyz = vtkSmartPointer<vtkDoubleArray>::New();
    xyz->SetNumberOfComponents(3);
    xyz->SetNumberOfTuples(_pointsN);
    xyz->SetName("xyz");
    xyz->SetComponentName(0, "x");
    xyz->SetComponentName(1, "y");
    xyz->SetComponentName(2, "z");
    xyz->SetArray(_nodesCoordinates.data(), _nodesCoordinates.size(), 1);

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    points->SetData(xyz);
    unstructuredGrid->SetPoints(points);

    if (type == "cells") {

        std::vector<std::vector<vtkIdType>> vtkCellNodes;
        for (auto const&[cell, nodes] : _cellsNodes) {
            vtkCellNodes.emplace_back();
            for (auto const &node : nodes)
                vtkCellNodes.back().push_back(node);
        }

        unstructuredGrid->Allocate(_cellsN);
        for (auto const &cellNodes : vtkCellNodes)
            unstructuredGrid->InsertNextCell(VTK_HEXAHEDRON, 8, cellNodes.data());

        for (auto &ent : _cellsArrays) {
            auto array = vtkSmartPointer<vtkDoubleArray>::New();
            array->SetNumberOfComponents(1);
            array->SetNumberOfTuples(_cellsN);
            array->SetName(ent.first.c_str());
            array->SetArray(_cellsArrays.at(ent.first).data(), _cellsN, 1);
            unstructuredGrid->GetCellData()->AddArray(array);
        }

    } else if (type == "faces") {

        std::vector<std::vector<vtkIdType>> vtkFaceNodes;
        for (auto const&[cell, nodes] : _facesNodes) {
            vtkFaceNodes.emplace_back();
            for (auto const &node : nodes)
                vtkFaceNodes.back().push_back(node);
        }

        unstructuredGrid->Allocate(_cellsN);
        for (auto const &faceNodes : vtkFaceNodes)
            unstructuredGrid->InsertNextCell(VTK_QUAD, 4, faceNodes.data());

        for (auto &ent : _facesArrays) {
            auto array = vtkSmartPointer<vtkDoubleArray>::New();
            array->SetNumberOfComponents(1);
            array->SetNumberOfTuples(_facesN);
            array->SetName(ent.first.c_str());
            array->SetArray(_facesArrays.at(ent.first).data(), _facesN, 1);
            unstructuredGrid->GetCellData()->AddArray(array);
        }
    }


    for (auto &ent : _pointsArrays) {
        auto array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(_pointsN);
        array->SetName(ent.first.c_str());
        array->SetArray(_pointsArrays.at(ent.first).data(), _pointsN, 1);
        unstructuredGrid->GetPointData()->AddArray(array);
    }

    auto unstructuredGridWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    unstructuredGridWriter->SetFileName(fileName.c_str());
    unstructuredGridWriter->SetInputData(unstructuredGrid);
    unstructuredGridWriter->SetDataModeToAscii();
    unstructuredGridWriter->Write();


}

void Sgrid::saveCells(const std::string &fileName) {
    save(fileName, "cells");
}

void Sgrid::saveFaces(const std::string &fileName) {
    save(fileName, "faces");
}


void saveFilesCollectionToFile(const std::string &fileName,
                               const std::vector<std::string> &filesNames,
                               const std::vector<std::string> &filesDescriptions) {

    pugi::xml_document doc;

    pugi::xml_node vtkFile = doc.append_child("VTKFile");
    vtkFile.append_attribute("type") = "Collection";
    vtkFile.append_attribute("version") = "1.0";
    vtkFile.append_attribute("byte_order") = "LittleEndian";
    vtkFile.append_attribute("header_type") = "UInt64";

    pugi::xml_node collection = vtkFile.append_child("Collection");

    for (int i = 0; i < filesNames.size(); i++) {

        collection.append_child("DataSet");

        collection.last_child().append_attribute("timestep").set_value(
                filesDescriptions[i].c_str());

        collection.last_child().append_attribute("file").set_value(
                filesNames[i].c_str());

    }

    doc.save_file(fileName.c_str());

}

