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
#include<set>
#include <numeric>


Sgrid::Sgrid(const std::vector<uint16_t> &pointsDims,
             const std::vector<double> &pointsOrigin,
             const std::vector<double> &spacing) :

        _pointsDims(pointsDims),
        _pointsN(pointsDims[0] * pointsDims[1] * pointsDims[2]),
        _pointsOrigin(pointsOrigin),
        _spacing(spacing),

        _cellsDims({uint16_t(pointsDims[0] - 1),
                    uint16_t(pointsDims[1] - 1),
                    uint16_t(pointsDims[2] - 1)}),
        _cellsN((pointsDims[0] - 1) * (pointsDims[1] - 1) * (pointsDims[2] - 1)),
        _cellV(spacing[0] * spacing[1] * spacing[2]) {

    calculateGridProps();

}


// Users methods

uint64_t Sgrid::calculateICell(const uint16_t &iX,
                               const uint16_t &iY,
                               const uint16_t &iZ) {
    return iX + iY * _cellsDims[0] + iZ * _cellsDims[0] * _cellsDims[1];
}

uint16_t Sgrid::calculateIXCell(const uint64_t &iCell) {
    return (iCell % (_cellsDims[0] * _cellsDims[1])) % _cellsDims[0];
}

uint16_t Sgrid::calculateIYCell(const uint64_t &iCell) {
    return (iCell % (_cellsDims[0] * _cellsDims[1])) / _cellsDims[0];
}

uint16_t Sgrid::calculateIZCell(const uint64_t &iCell) {
    return iCell / (_cellsDims[0] * _cellsDims[1]);
}


uint8_t Sgrid::calculateAxisFace(const uint64_t &iFace) {
    uint8_t axis = 0;
    if (iFace >= _facesNs[0] + _facesNs[1])
        axis = 2;
    else if (iFace >= _facesNs[0])
        axis = 1;
    return axis;
}

uint64_t Sgrid::calculateOffsetFace(const uint8_t &axis) {
    uint64_t offset = 0;
    if (axis == 1)
        offset = _facesNs[0];
    else if (axis == 2)
        offset = _facesNs[0] + _facesNs[1];
    return offset;
}

uint64_t Sgrid::calculateIFace(const uint8_t &axis,
                               const uint16_t &iX, const uint16_t &iY,
                               const uint16_t &iZ) {
    return calculateOffsetFace(axis) + iX + iY * _facesDimss[axis][0] +
           iZ * _facesDimss[axis][0] * _facesDimss[axis][1];
}


uint16_t Sgrid::calculateIXFace(const uint64_t &iFace) {
    auto axis = calculateAxisFace(iFace);
    auto offset = calculateOffsetFace(axis);
    return ((iFace - offset) % (_facesDimss[axis][0] * _facesDimss[axis][1])) %
           _facesDimss[axis][0];
}

uint64_t Sgrid::calculateIYFace(const uint16_t &iFace) {
    auto axis = calculateAxisFace(iFace);
    auto offset = calculateOffsetFace(axis);
    return ((iFace - offset) % (_facesDimss[axis][0] * _facesDimss[axis][1])) /
           _facesDimss[axis][0];
}

uint64_t Sgrid::calculateIZFace(const uint16_t &iFace) {
    auto axis = calculateAxisFace(iFace);
    auto offset = calculateOffsetFace(axis);
    return (iFace - offset) / (_facesDimss[axis][0] * _facesDimss[axis][1]);
}


void Sgrid::setCellsType(const std::string &name, Eigen::Ref<Eigen::VectorXui64> cells) {
    _typesCells.erase(name);
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            name, Eigen::Map<Eigen::VectorXui64>(cells.data(), cells.size())));
}

void Sgrid::processTypesByCellsType(const std::string &name) {

    std::set<uint64_t> facesBound;
    std::set<uint64_t> facesNonbound;
    for (uint64_t i = 0; i < _typesCells.at(name).size(); i++)
        for (auto &face : _neighborsFaces[_typesCells.at(name)(i)])
            if (facesBound.find(face) == facesBound.end())
                facesBound.insert(face);
            else {
                facesBound.erase(face);
                facesNonbound.insert(face);
            }


    _typesFaces.erase(name + "_bound");

    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            name + "_bound", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    copyStdSetToEigenVector<uint64_t>(facesBound, _typesFaces.at(name + "_bound"));

    _typesFaces.erase(name + "_nonbound");
    if (!facesNonbound.empty()) {
        _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
                name + "_nonbound", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
        copyStdSetToEigenVector<uint64_t>(facesNonbound,
                                          _typesFaces.at(name + "_nonbound"));
    }


    std::set<uint64_t> cellsNeighborBound;
    for (auto &face : facesBound)
        for (auto &cell : _neighborsCells[face])
            cellsNeighborBound.insert(cell);

    std::set<uint64_t> cells;
    for (int i = 0; i < _typesCells.at(name).size(); i++)
        cells.insert(_typesCells.at(name)(i));

    std::set<uint64_t> cellsBound;
    std::set_intersection(cells.begin(), cells.end(),
                          cellsNeighborBound.begin(), cellsNeighborBound.end(),
                          std::inserter(cellsBound, cellsBound.begin()));

    std::set<uint64_t> cellsNonbound;
    std::set_difference(cells.begin(), cells.end(),
                        cellsBound.begin(), cellsBound.end(),
                        std::inserter(cellsNonbound, cellsNonbound.begin()));


    _typesCells.erase(name + "_bound");
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            name + "_bound", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    copyStdSetToEigenVector<uint64_t>(cellsBound, _typesCells.at(name + "_bound"));

    _typesCells.erase(name + "_nonbound");
    if (!cellsNonbound.empty()) {
        _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
                name + "_nonbound", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
        copyStdSetToEigenVector<uint64_t>(cellsNonbound,
                                          _typesCells.at(name + "_nonbound"));
    }


}


// Accessory shitty constructor methods

void Sgrid::calculateFacessDims() {
    for (int8_t i = 0; i < 3; ++i) {
        auto dims = _pointsDims;
        dims[(i + 2) % 3] -= 1;
        dims[(i + 1) % 3] -= 1;
        _facesDimss.push_back(dims);
    }
}

void Sgrid::calculateFacesNs() {
    _facesNs = {uint64_t(_facesDimss[0][0] * _facesDimss[0][1] * _facesDimss[0][2]),
                uint64_t(_facesDimss[1][0] * _facesDimss[1][1] * _facesDimss[1][2]),
                uint64_t(_facesDimss[2][0] * _facesDimss[2][1] * _facesDimss[2][2])};
}

void Sgrid::calculateFacesN() {
    _facesN = _facesNs[0] + _facesNs[1] + _facesNs[2];
}

void Sgrid::calculateFacesAxes() {
    _facesAxes = std::vector<uint8_t>(_facesN);
    for (uint64_t i = 0; i < _facesN; i++)
        if (i < _facesNs[0])
            _facesAxes[i] = 0;
        else if (i >= _facesNs[0] + _facesNs[1])
            _facesAxes[i] = 2;
        else
            _facesAxes[i] = 1;
}


void Sgrid::calculateFaceSs() {
    _facesSs = {_spacing[1] * _spacing[2],
                _spacing[0] * _spacing[2],
                _spacing[0] * _spacing[1]};
}

void Sgrid::calculateNeighborsFaces() {

    for (uint64_t iCell = 0; iCell < _cellsN; iCell++)
        _neighborsFaces[iCell] = std::vector<uint64_t>(6);

    uint64_t offset0 = 0;
    uint64_t offset1 = _facesNs[0];
    uint64_t offset2 = _facesNs[0] + _facesNs[1];

    for (uint16_t k = 0; k < _cellsDims[2]; ++k) {
        for (uint16_t j = 0; j < _cellsDims[1]; ++j) {
            for (uint16_t i = 0; i < _cellsDims[0]; ++i) {
                auto iCell = calculateICell(i, j, k);
                _neighborsFaces[iCell][0] = iCell + offset0;
                _neighborsFaces[iCell][1] = iCell + 1 + offset0;
                _neighborsFaces[iCell][2] = iCell + offset1;
                _neighborsFaces[iCell][3] = iCell + _cellsDims[0] + offset1;
                _neighborsFaces[iCell][4] = iCell + offset2;
                _neighborsFaces[iCell][5] =
                        iCell + _cellsDims[0] * _cellsDims[1] + offset2;
            }
            offset0++;
        }
        offset1 += _cellsDims[0];
    }

}

void Sgrid::calculateNeighborsCells() {

    for (uint64_t iFace = 0; iFace < _facesN; iFace++)
        _neighborsCells[iFace] = std::vector<uint64_t>();

    for (auto &[cell, faces] : _neighborsFaces)
        for (auto &face : faces)
            _neighborsCells[face].push_back(cell);

}


void Sgrid::calculateNormalsNeighborsCells() {

    for (uint64_t iCell = 0; iCell < _cellsN; iCell++)
        _normalsNeighborsCells[iCell] = std::vector<int8_t>();

    for (auto &[face, cells] : _neighborsCells)
        for (uint8_t i = 0; i < cells.size(); i++)
            _normalsNeighborsCells[face].push_back(2 * i - 1);

}

void Sgrid::calculateNormalsNeighborsFaces() {

    for (uint64_t iCell = 0; iCell < _cellsN; iCell++)
        _normalsNeighborsFaces[iCell] = std::vector<int8_t>(6);

    for (auto &[cell, faces] : _neighborsFaces)
        for (uint8_t i = 0; i < faces.size(); i++) {

            uint8_t j;
            for (j = 0; j < _neighborsCells[faces[i]].size(); j++)
                if (_neighborsCells[faces[i]][j] == cell)
                    break;

            if (_neighborsCells[faces[i]][j] != cell) {
                std::cerr << "Inconsistent stuff with normalsNeighborsFaces."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }

            _normalsNeighborsFaces[cell][i] = _normalsNeighborsCells[faces[i]][j];

        }

}

void Sgrid::calculateMainTypesCells() {

    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "left", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "right", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "front", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "back", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "bottom", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "top", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));

    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "nonbound", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));


    bool isX1D = (_cellsDims[0] == 1);
    bool isY1D = (_cellsDims[1] == 1);
    bool isZ1D = (_cellsDims[2] == 1);

    if (isX1D or isY1D or isZ1D)
        _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
                "nonbound_1D", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));


    std::vector<uint64_t> left;
    std::vector<uint64_t> right;
    std::vector<uint64_t> front;
    std::vector<uint64_t> back;
    std::vector<uint64_t> bottom;
    std::vector<uint64_t> top;

    std::vector<uint64_t> nonbound;
    std::vector<uint64_t> nonbound1D;


    for (uint64_t iCell = 0; iCell < _cellsN; iCell++) {

        bool isIXCellNonbound = true;
        auto iX = calculateIXCell(iCell);
        if (iX == 0) {
            left.push_back(iCell);
            isIXCellNonbound = false;
        }
        if (iX == _cellsDims[0] - 1) {
            right.push_back(iCell);
            isIXCellNonbound = false;
        }

        bool isIYCellNonbound = true;
        auto iY = calculateIYCell(iCell);
        if (iY == 0) {
            front.push_back(iCell);
            isIYCellNonbound = false;
        }
        if (iY == _cellsDims[1] - 1) {
            back.push_back(iCell);
            isIYCellNonbound = false;
        }

        bool isIZCellNonbound = true;
        auto iZ = calculateIZCell(iCell);
        if (iZ == 0) {
            bottom.push_back(iCell);
            isIZCellNonbound = false;
        }
        if (iZ == _cellsDims[2] - 1) {
            top.push_back(iCell);
            isIZCellNonbound = false;
        }

        if (isIXCellNonbound and
            isIYCellNonbound and isIZCellNonbound)
            nonbound.push_back(iCell);

        if ((isIXCellNonbound or isX1D) and
            (isIYCellNonbound or isY1D) and (isIZCellNonbound or isZ1D))
            nonbound1D.push_back(iCell);

    }

    copyStdVectotToEigenVector<uint64_t>(left, _typesCells.at("left"));
    copyStdVectotToEigenVector<uint64_t>(right, _typesCells.at("right"));
    copyStdVectotToEigenVector<uint64_t>(front, _typesCells.at("front"));
    copyStdVectotToEigenVector<uint64_t>(back, _typesCells.at("back"));
    copyStdVectotToEigenVector<uint64_t>(bottom, _typesCells.at("bottom"));
    copyStdVectotToEigenVector<uint64_t>(top, _typesCells.at("top"));

    copyStdVectotToEigenVector<uint64_t>(nonbound, _typesCells.at("nonbound"));
    if (isX1D or isY1D or isZ1D)
        copyStdVectotToEigenVector<uint64_t>(nonbound1D, _typesCells.at("nonbound_1D"));

}

void Sgrid::calculateMainTypesFaces() {

    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "left", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "right", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "front", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "back", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "bottom", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "top", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));

    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXui64>>(
            "nonbound", Eigen::Map<Eigen::VectorXui64>(new uint64_t[1], 1)));


    std::vector<uint64_t> left;
    std::vector<uint64_t> right;
    std::vector<uint64_t> front;
    std::vector<uint64_t> back;
    std::vector<uint64_t> bottom;
    std::vector<uint64_t> top;

    std::set<uint64_t> nonboundSet;

    for (uint64_t i = 0; i < _typesCells.at("left").size(); i++) {
        auto iCell = _typesCells.at("left")(i);
        auto &faces = _neighborsFaces.at(iCell);
        auto buffer = std::vector<uint64_t>();
        for (auto &face : faces)
            if (_neighborsCells[face].size() == 1 and _facesAxes[face] == 0)
                buffer.push_back(face);
            else if (_neighborsCells[face].size() != 1)
                nonboundSet.insert(face);
        if (buffer.size() == 2)
            if (calculateIXFace(buffer.front()) > calculateIXFace(buffer.back()))
                left.push_back(buffer.back());
            else
                left.push_back(buffer.front());
        else
            left.push_back(buffer.front());
        buffer.clear();
    }

    for (uint64_t i = 0; i < _typesCells.at("right").size(); i++) {
        auto iCell = _typesCells.at("right")(i);
        auto &faces = _neighborsFaces.at(iCell);
        auto buffer = std::vector<uint64_t>();
        for (auto &face : faces)
            if (_neighborsCells[face].size() == 1 and _facesAxes[face] == 0)
                buffer.push_back(face);
            else if (_neighborsCells[face].size() != 1)
                nonboundSet.insert(face);
        if (buffer.size() == 2)
            if (calculateIXFace(buffer.front()) < calculateIXFace(buffer.back()))
                right.push_back(buffer.back());
            else
                right.push_back(buffer.front());
        else
            right.push_back(buffer.front());
        buffer.clear();
    }

    for (uint64_t i = 0; i < _typesCells.at("front").size(); i++) {
        auto iCell = _typesCells.at("front")(i);
        auto &faces = _neighborsFaces.at(iCell);
        auto buffer = std::vector<uint64_t>();
        for (auto &face : faces)
            if (_neighborsCells[face].size() == 1 and _facesAxes[face] == 1)
                buffer.push_back(face);
            else if (_neighborsCells[face].size() != 1)
                nonboundSet.insert(face);
        if (buffer.size() == 2)
            if (calculateIYFace(buffer.front()) > calculateIYFace(buffer.back()))
                front.push_back(buffer.back());
            else
                front.push_back(buffer.front());
        else
            front.push_back(buffer.front());
        buffer.clear();
    }

    for (uint64_t i = 0; i < _typesCells.at("back").size(); i++) {
        auto iCell = _typesCells.at("back")(i);
        auto &faces = _neighborsFaces.at(iCell);
        auto buffer = std::vector<uint64_t>();
        for (auto &face : faces)
            if (_neighborsCells[face].size() == 1 and _facesAxes[face] == 1)
                buffer.push_back(face);
            else if (_neighborsCells[face].size() != 1)
                nonboundSet.insert(face);
        if (buffer.size() == 2)
            if (calculateIYFace(buffer.front()) < calculateIYFace(buffer.back()))
                back.push_back(buffer.back());
            else
                back.push_back(buffer.front());
        else
            back.push_back(buffer.front());
        buffer.clear();
    }

    for (uint64_t i = 0; i < _typesCells.at("bottom").size(); i++) {
        auto iCell = _typesCells.at("bottom")(i);
        auto &faces = _neighborsFaces.at(iCell);
        auto buffer = std::vector<uint64_t>();
        for (auto &face : faces)
            if (_neighborsCells[face].size() == 1 and _facesAxes[face] == 2)
                buffer.push_back(face);
            else if (_neighborsCells[face].size() != 1)
                nonboundSet.insert(face);
        if (buffer.size() == 2)
            if (calculateIZFace(buffer.front()) > calculateIZFace(buffer.back()))
                bottom.push_back(buffer.back());
            else
                bottom.push_back(buffer.front());
        else
            bottom.push_back(buffer.front());
        buffer.clear();
    }

    for (uint64_t i = 0; i < _typesCells.at("top").size(); i++) {
        auto iCell = _typesCells.at("top")(i);
        auto &faces = _neighborsFaces.at(iCell);
        auto buffer = std::vector<uint64_t>();
        for (auto &face : faces)
            if (_neighborsCells[face].size() == 1 and _facesAxes[face] == 2)
                buffer.push_back(face);
            else if (_neighborsCells[face].size() != 1)
                nonboundSet.insert(face);
        if (buffer.size() == 2)
            if (calculateIZFace(buffer.front()) < calculateIZFace(buffer.back()))
                top.push_back(buffer.back());
            else
                top.push_back(buffer.front());
        else
            top.push_back(buffer.front());
        buffer.clear();
    }


    for (int i = 0; i < _typesCells.at("nonbound").size(); i++) {
        auto iCell = _typesCells.at("nonbound")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (auto &face : faces)
            nonboundSet.insert(face);
    }


    copyStdVectotToEigenVector<uint64_t>(left, _typesFaces.at("left"));
    copyStdVectotToEigenVector<uint64_t>(right, _typesFaces.at("right"));
    copyStdVectotToEigenVector<uint64_t>(front, _typesFaces.at("front"));
    copyStdVectotToEigenVector<uint64_t>(back, _typesFaces.at("back"));
    copyStdVectotToEigenVector<uint64_t>(bottom, _typesFaces.at("bottom"));
    copyStdVectotToEigenVector<uint64_t>(top, _typesFaces.at("top"));

    if (!nonboundSet.empty())
        copyStdSetToEigenVector<uint64_t>(nonboundSet, _typesFaces.at("nonbound"));

}


void Sgrid::calculateGridProps() {

    calculateFacessDims();
    calculateFacesNs();
    calculateFacesN();
    calculateFacesAxes();
    calculateFaceSs();

    calculateNeighborsFaces();
    calculateNeighborsCells();

    calculateNormalsNeighborsCells();
    calculateNormalsNeighborsFaces();

    calculateMainTypesCells();
    calculateMainTypesFaces();

}
