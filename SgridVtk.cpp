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
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkThreshold.h>

#include <pugixml.hpp>

void Sgrid::save(const std::string &fileName,
                 const bool &isStructured=false) {

    // Create VtkStructuredGrid

    auto structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(_pointsDims[0], _pointsDims[1], _pointsDims[2]);

    // Create point and incorporate it into VtkStructuredGrid

    auto xyzPtr = new double[_pointsN * 3];
    for (int k = 0; k < _pointsDims[2]; k++)
        for (int j = 0; j < _pointsDims[1]; j++)
            for (int i = 0; i < _pointsDims[0]; i++) {
                int scalarInd =
                        i + j * _pointsDims[0] + k * _pointsDims[0] * _pointsDims[1];
                xyzPtr[scalarInd * 3] = _pointsOrigin[0] + _spacing[0] * i;
                xyzPtr[scalarInd * 3 + 1] = _pointsOrigin[1] + _spacing[1] * j;
                xyzPtr[scalarInd * 3 + 2] = _pointsOrigin[2] + _spacing[2] * k;
            }


    auto xyz = vtkSmartPointer<vtkDoubleArray>::New();
    xyz->SetNumberOfComponents(3);
    xyz->SetNumberOfTuples(_pointsN);
    xyz->SetName("xyz");
    xyz->SetComponentName(0, "x");
    xyz->SetComponentName(1, "y");
    xyz->SetComponentName(2, "z");
    xyz->SetArray(xyzPtr, xyz->GetDataSize(), 1);


    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    points->SetData(xyz);
    structuredGrid->SetPoints(points);


    // Add arrays to points and cells

    for (auto &ent : _pointsArrays) {
        auto array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(_pointsN);
        array->SetName(ent.first.c_str());
        array->SetArray(_pointsArrays.at(ent.first).data(), _pointsN, 1);
        structuredGrid->GetPointData()->AddArray(array);
    }

    for (auto &ent : _cellsArrays) {
        auto array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(_cellsN);
        array->SetName(ent.first.c_str());
        array->SetArray(_cellsArrays.at(ent.first).data(), _cellsN, 1);
        structuredGrid->GetCellData()->AddArray(array);
    }


    if (isStructured) {

        // Write VtkStructuredGrid to a file

        auto structuredGridWriter = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        structuredGridWriter->SetFileName(fileName.c_str());
        structuredGridWriter->SetInputData(structuredGrid);
        structuredGridWriter->SetDataModeToAscii();
        structuredGridWriter->Write();

    } else {

        // Converts VtkStructuredGrid into VtkUnstructuredGrid

        // Add just scalar to the points
        auto scalarsPointsPtr = new int[_pointsN];
        for (int i = 0; i < _pointsN; i++)
            scalarsPointsPtr[i] = 0;
        auto scalarsPoints = vtkSmartPointer<vtkIntArray>::New();
        scalarsPoints->SetNumberOfComponents(1);
        scalarsPoints->SetNumberOfTuples(_pointsN);
        scalarsPoints->SetName("scalar");
        scalarsPoints->SetArray(scalarsPointsPtr, _pointsN, 1);
        structuredGrid->GetPointData()->SetScalars(scalarsPoints);

        auto threshold = vtkSmartPointer<vtkThreshold>::New();
        threshold->SetInputData(structuredGrid);
        threshold->ThresholdByUpper(-1);
        threshold->Update();

        // Write vtkUnstructuredGrid to a file

        auto unstructuredGridWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        unstructuredGridWriter->SetFileName(fileName.c_str());
        unstructuredGridWriter->SetDataModeToAscii();
        unstructuredGridWriter->SetInputConnection(threshold->GetOutputPort());
        unstructuredGridWriter->Write();

    }

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

