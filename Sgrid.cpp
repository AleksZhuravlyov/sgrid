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

#include<string>
#include<vector>
#include<map>

#include <pugixml.hpp>


Sgrid::Sgrid(Eigen::Ref<Eigen::Vector3i> pointsDims,
             Eigen::Ref<Eigen::Vector3d> pointsOrigin,
             Eigen::Ref<Eigen::Vector3d> spacing) :

    _pointsDims(pointsDims.data(), pointsDims.size()),
    _pointsOrigin(pointsOrigin.data(), pointsOrigin.size()),
    _spacing(spacing.data(), spacing.size()),
    _cellsDims(Eigen::Map<Eigen::Vector3i>(new int[3])),
    _faceS(Eigen::Map<Eigen::Vector3d>(new double[3])) {

  calculateGridProps();

}


void Sgrid::calculatePointsN() {
  _pointsN = _pointsDims.prod();
}


void Sgrid::calculateCellsDims() {
  _cellsDims = _pointsDims - Eigen::Vector3i::Ones(3);
}

void Sgrid::calculateCellsN() {
  _cellsN = _cellsDims.prod();
}


void Sgrid::calculateFacesDims() {
  for (int i = 0; i < 3; ++i) {
    auto dims = Eigen::Map<Eigen::Vector3i>(new int[3]);
    dims = _pointsDims;
    dims((i + 2) % 3) -= 1;
    dims((i + 1) % 3) -= 1;
    _facesDims.insert(std::pair<int, Eigen::Map<Eigen::Vector3i>>(i, dims));
  }
}

void Sgrid::calculateFacesN() {
  _facesN = _facesDims.at(0).prod() + _facesDims.at(1).prod() + _facesDims.at(2).prod();
}


void Sgrid::calculateCellV() {
  _cellV = _spacing.prod();
}

void Sgrid::calculateFaceS() {
  _faceS = Eigen::Vector3d(_spacing(1) * _spacing(2),
                           _spacing(0) * _spacing(2),
                           _spacing(0) * _spacing(1));
}


void Sgrid::calculateNeighborsFaces() {

  for (int iCell = 0; iCell < _cellsN; iCell++)
    _neighborsFaces.insert(std::pair<int, Eigen::Map<Eigen::VectorXi>>(
        iCell, Eigen::Map<Eigen::VectorXi>(new int[6], 6)));


  int offset0 = 0;
  int offset1 = _facesDims.at(0).prod();
  int offset2 = _facesDims.at(0).prod() + _facesDims.at(1).prod();

  for (int k = 0; k < _cellsDims(2); ++k) {
    for (int j = 0; j < _cellsDims(1); ++j) {
      for (int i = 0; i < _cellsDims(0); ++i) {

        auto iCell = i + j * _cellsDims(0) + k * _cellsDims(0) * _cellsDims(1);

        _neighborsFaces.at(iCell)(0) = iCell + offset0;
        _neighborsFaces.at(iCell)(1) = iCell + 1 + offset0;
        _neighborsFaces.at(iCell)(2) = iCell + offset1;
        _neighborsFaces.at(iCell)(3) = iCell + _cellsDims(0) + offset1;
        _neighborsFaces.at(iCell)(4) = iCell + offset2;
        _neighborsFaces.at(iCell)(5) = iCell + _cellsDims(0) * _cellsDims(1) + offset2;

      }
      offset0++;
    }
    offset1 += _cellsDims(0);
  }
}

void Sgrid::calculateNeighborsCells() {

  std::map<int, std::vector<int>> neighborsCellsStd;
  for (int iFace = 0; iFace < _facesN; iFace++)
    neighborsCellsStd[iFace] = std::vector<int>();

  for (const auto &ent : _neighborsFaces)
    for (int i = 0; i < _neighborsFaces.at(ent.first).size(); i++)
      neighborsCellsStd[_neighborsFaces.at(ent.first)(i)].push_back(ent.first);

  for (int iFace = 0; iFace < _facesN; iFace++) {

    auto size = neighborsCellsStd[iFace].size();
    int *array = new int[size];
    for (int i = 0; i < size; i++)
      array[i] = neighborsCellsStd[iFace][i];

    _neighborsCells.insert(std::pair<int, Eigen::Map<Eigen::VectorXi>>(
        iFace, Eigen::Map<Eigen::VectorXi>(array, size)));
  }

}


void Sgrid::calculateNormalsNeighborsCells() {

  for (auto &ent : _neighborsCells) {
    auto size = _neighborsCells.at(ent.first).size();
    _normalsNeighborsCells.insert(std::pair<int, Eigen::Map<Eigen::VectorXi>>(
        ent.first, Eigen::Map<Eigen::VectorXi>(new int[size], size)));
  }

  for (auto &ent : _normalsNeighborsCells)
    for (int i = 0; i < _normalsNeighborsCells.at(ent.first).size(); i++)
      _normalsNeighborsCells.at(ent.first)(i) = 2 * i - 1;

}

void Sgrid::calculateNormalsNeighborsFaces() {

  for (auto &ent : _neighborsFaces)
    _normalsNeighborsFaces.insert(std::pair<int, Eigen::Map<Eigen::VectorXi>>(
        ent.first, Eigen::Map<Eigen::VectorXi>(new int[6], 6)));


  for (auto const &ent : _neighborsFaces)
    for (int i = 0; i < _neighborsFaces.at(ent.first).size(); i++) {

      int j;
      for (j = 0; j < _neighborsCells.at(_neighborsFaces.at(ent.first)(i)).size(); j++)

        if (_neighborsCells.at(_neighborsFaces.at(ent.first)(i))(j) == ent.first)
          break;

      if (_neighborsCells.at(_neighborsFaces.at(ent.first)(i))(j) != ent.first) {
        std::cerr << "Inconsistent stuff with normalsNeighborsFaces." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      _normalsNeighborsFaces.at(ent.first)(i) =
          _normalsNeighborsCells.at(_neighborsFaces.at(ent.first)(i))(j);

    }


}

void Sgrid::calculateGridProps() {

  calculatePointsN();

  calculateCellsDims();
  calculateCellsN();

  calculateFacesDims();
  calculateFacesN();

  calculateCellV();
  calculateFaceS();

  calculateNeighborsFaces();
  calculateNeighborsCells();

  calculateNormalsNeighborsCells();
  calculateNormalsNeighborsFaces();

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
