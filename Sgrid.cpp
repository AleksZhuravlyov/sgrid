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


int Sgrid::calculateICell(const int &iX, const int &iY, const int &iZ) {
    return iX + iY * _cellsDims(0) + iZ * _cellsDims(0) * _cellsDims(1);
}

int Sgrid::calculateIXCell(const int &iCell) {
    return (iCell % (_cellsDims(0) * _cellsDims(1))) % _cellsDims(0);
}

int Sgrid::calculateIYCell(const int &iCell) {
    return (iCell % (_cellsDims(0) * _cellsDims(1))) / _cellsDims(0);
}

int Sgrid::calculateIZCell(const int &iCell) {
    return iCell / (_cellsDims(0) * _cellsDims(1));
}


void Sgrid::setCellsType(const std::string &name, Eigen::Ref<Eigen::VectorXi> cells) {

    _typesCells.erase(name);

    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            name, Eigen::Map<Eigen::VectorXi>(cells.data(), cells.size())));

}

void Sgrid::processFacesTypesByCellsType(const std::string &name) {

    std::set<int> setFacesBound;
    std::set<int> setFacesNonbound;
    auto &cells = _typesCells.at(name);
    for (int i = 0; i < cells.size(); i++) {
        auto &faces = _neighborsFaces.at(cells(i));
        for (int j = 0; j < faces.size(); j++) {
            auto face = faces(j);
            if (setFacesBound.find(face) == setFacesBound.end())
                setFacesBound.insert(face);
            else {
                setFacesBound.erase(face);
                setFacesNonbound.insert(face);
            }
        }
    }

    _typesFaces.erase(name + "_bound");
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            name + "_bound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    std::vector<int> vectorFacesBound(setFacesBound.size());
    std::copy(setFacesBound.begin(), setFacesBound.end(), vectorFacesBound.begin());
    copyStdVectotToEigenVector<int>(vectorFacesBound, _typesFaces.at(name + "_bound"));

    _typesFaces.erase(name + "_nonbound");
    if (!setFacesNonbound.empty()) {
        _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
                name + "_nonbound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
        std::vector<int> vectorFacesNonbound(setFacesNonbound.size());
        std::copy(setFacesNonbound.begin(), setFacesNonbound.end(),
                  vectorFacesNonbound.begin());
        copyStdVectotToEigenVector<int>(vectorFacesNonbound,
                                        _typesFaces.at(name + "_nonbound"));
    }

    std::set<int> setCellsNeighborBound;
    for (auto &face : setFacesBound)
        for (int i = 0; i < _neighborsCells.at(face).size(); i++)
            setCellsNeighborBound.insert(_neighborsCells.at(face)(i));
    std::set<int> setCells;
    for (int i = 0; i < cells.size(); i++)
        setCells.insert(_typesCells.at(name)(i));
    std::set<int> setCellsBound;
    std::set_intersection(setCells.begin(), setCells.end(),
                          setCellsNeighborBound.begin(), setCellsNeighborBound.end(),
                          std::inserter(setCellsBound, setCellsBound.begin()));
    std::set<int> setCellsNonbound;
    std::set_difference(setCells.begin(), setCells.end(),
                        setCellsBound.begin(), setCellsBound.end(),
                        std::inserter(setCellsNonbound, setCellsNonbound.begin()));


    _typesCells.erase(name + "_bound");
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            name + "_bound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    std::vector<int> vectorCellsBound(setCellsBound.size());
    std::copy(setCellsBound.begin(), setCellsBound.end(), vectorCellsBound.begin());
    copyStdVectotToEigenVector<int>(vectorCellsBound, _typesCells.at(name + "_bound"));

    _typesCells.erase(name + "_nonbound");
    if (!setCellsNonbound.empty()) {
        _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
                name + "_nonbound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
        std::vector<int> vectorCellsNonbound(setCellsNonbound.size());
        std::copy(setCellsNonbound.begin(), setCellsNonbound.end(),
                  vectorCellsNonbound.begin());
        copyStdVectotToEigenVector<int>(vectorCellsNonbound,
                                        _typesCells.at(name + "_nonbound"));
    }


    /*std::set<int> setCellsNonbound;
    std::set<int> setCellsBound;
    for (auto &face : setFacesNonbound)
        for (int i = 0; i < _neighborsCells.at(face).size(); i++)
            setCellsNonbound.insert(_neighborsCells.at(face)(i));
    std::set<int> setCells;
    for (int i = 0; i < cells.size(); i++)
        setCells.insert(_typesCells.at(name)(i));
    std::set_difference(setCells.begin(), setCells.end(),
                        setCellsNonbound.begin(), setCellsNonbound.end(),
                        std::inserter(setCellsBound, setCellsBound.begin()));

    _typesCells.erase(name + "_bound");
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            name + "_bound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    std::vector<int> vectorCellsBound(setCellsBound.size());
    std::copy(setCellsBound.begin(), setCellsBound.end(), vectorCellsBound.begin());
    copyStdVectotToEigenVector<int>(vectorCellsBound, _typesCells.at(name + "_bound"));

    _typesCells.erase(name + "_nonbound");
    if (!setCellsNonbound.empty()) {
        _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
                name + "_nonbound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
        std::vector<int> vectorCellsNonbound(setCellsNonbound.size());
        std::copy(setCellsNonbound.begin(), setCellsNonbound.end(),
                  vectorCellsNonbound.begin());
        copyStdVectotToEigenVector<int>(vectorCellsNonbound,
                                        _typesCells.at(name + "_nonbound"));
    }*/


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

void Sgrid::calculateFacesAxis() {
    for (int i = 0; i < _facesN; i++)
        if (i < _facesDims.at(0).prod())
            _facesAxis[i] = 0;
        else if (i >= _facesDims.at(0).prod() + _facesDims.at(1).prod())
            _facesAxis[i] = 2;
        else
            _facesAxis[i] = 1;
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

                auto iCell = calculateICell(i, j, k);

                _neighborsFaces.at(iCell)(0) = iCell + offset0;
                _neighborsFaces.at(iCell)(1) = iCell + 1 + offset0;
                _neighborsFaces.at(iCell)(2) = iCell + offset1;
                _neighborsFaces.at(iCell)(3) = iCell + _cellsDims(0) + offset1;
                _neighborsFaces.at(iCell)(4) = iCell + offset2;
                _neighborsFaces.at(iCell)(5) =
                        iCell + _cellsDims(0) * _cellsDims(1) + offset2;

            }
            offset0++;
        }
        offset1 += _cellsDims(0);
    }
}

void Sgrid::calculateNeighborsCells() {

    std::map<int, std::vector<int>> neighbors;
    for (int iFace = 0; iFace < _facesN; iFace++)
        neighbors[iFace] = std::vector<int>();

    for (const auto &ent : _neighborsFaces)
        for (int i = 0; i < _neighborsFaces.at(ent.first).size(); i++)
            neighbors[_neighborsFaces.at(ent.first)(i)].push_back(ent.first);

    for (int iFace = 0; iFace < _facesN; iFace++) {

        _neighborsCells.insert(std::pair<int, Eigen::Map<Eigen::VectorXi>>(
                iFace, Eigen::Map<Eigen::VectorXi>(new int[1], 1)));

        copyStdVectotToEigenVector<int>(neighbors[iFace], _neighborsCells.at(iFace));
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
            for (j = 0;
                 j < _neighborsCells.at(_neighborsFaces.at(ent.first)(i)).size(); j++)

                if (_neighborsCells.at(_neighborsFaces.at(ent.first)(i))(j) == ent.first)
                    break;

            if (_neighborsCells.at(_neighborsFaces.at(ent.first)(i))(j) != ent.first) {
                std::cerr << "Inconsistent stuff with normalsNeighborsFaces."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }

            _normalsNeighborsFaces.at(ent.first)(i) =
                    _normalsNeighborsCells.at(_neighborsFaces.at(ent.first)(i))(j);

        }


}

void Sgrid::calculateMainTypesCells() {

    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "left", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "right", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "front", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "back", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "bottom", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "top", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));


    _typesCells.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "nonbound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));

    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> front;
    std::vector<int> back;
    std::vector<int> bottom;
    std::vector<int> top;

    std::vector<int> nonbound;

    for (int iCell = 0; iCell < _cellsN; iCell++) {

        bool isICellNonbound = true;

        auto iX = calculateIXCell(iCell);
        if (iX == 0) {
            left.push_back(iCell);
            isICellNonbound = false;
        }
        if (iX == _cellsDims[0] - 1) {
            right.push_back(iCell);
            isICellNonbound = false;
        }

        auto iY = calculateIYCell(iCell);
        if (iY == 0) {
            front.push_back(iCell);
            isICellNonbound = false;
        }
        if (iY == _cellsDims[1] - 1) {
            back.push_back(iCell);
            isICellNonbound = false;
        }

        auto iZ = calculateIZCell(iCell);
        if (iZ == 0) {
            bottom.push_back(iCell);
            isICellNonbound = false;
        }
        if (iZ == _cellsDims[2] - 1) {
            top.push_back(iCell);
            isICellNonbound = false;
        }

        if (isICellNonbound)
            nonbound.push_back(iCell);

    }

    copyStdVectotToEigenVector<int>(left, _typesCells.at("left"));
    copyStdVectotToEigenVector<int>(right, _typesCells.at("right"));
    copyStdVectotToEigenVector<int>(front, _typesCells.at("front"));
    copyStdVectotToEigenVector<int>(back, _typesCells.at("back"));
    copyStdVectotToEigenVector<int>(bottom, _typesCells.at("bottom"));
    copyStdVectotToEigenVector<int>(top, _typesCells.at("top"));

    copyStdVectotToEigenVector<int>(nonbound, _typesCells.at("nonbound"));

}

void Sgrid::calculateMainTypesFaces() {

    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "left", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "right", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "front", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "back", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "bottom", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));
    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "top", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));

    _typesFaces.insert(std::pair<std::string, Eigen::Map<Eigen::VectorXi>>(
            "nonbound", Eigen::Map<Eigen::VectorXi>(new int[1], 1)));

    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> front;
    std::vector<int> back;
    std::vector<int> bottom;
    std::vector<int> top;

    std::set<int> nonboundSet;

    for (int i = 0; i < _typesCells.at("left").size(); i++) {
        auto iCell = _typesCells.at("left")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            if (_neighborsCells.at(faces(j)).size() == 1 and _facesAxis[faces(j)] == 0)
                left.push_back(faces(j));
    }

    for (int i = 0; i < _typesCells.at("right").size(); i++) {
        auto iCell = _typesCells.at("right")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            if (_neighborsCells.at(faces(j)).size() == 1 and _facesAxis[faces(j)] == 0)
                right.push_back(faces(j));
    }

    for (int i = 0; i < _typesCells.at("front").size(); i++) {
        auto iCell = _typesCells.at("front")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            if (_neighborsCells.at(faces(j)).size() == 1 and _facesAxis[faces(j)] == 1)
                front.push_back(faces(j));
    }

    for (int i = 0; i < _typesCells.at("back").size(); i++) {
        auto iCell = _typesCells.at("back")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            if (_neighborsCells.at(faces(j)).size() == 1 and _facesAxis[faces(j)] == 1)
                back.push_back(faces(j));
    }

    for (int i = 0; i < _typesCells.at("bottom").size(); i++) {
        auto iCell = _typesCells.at("bottom")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            if (_neighborsCells.at(faces(j)).size() == 1 and _facesAxis[faces(j)] == 2)
                bottom.push_back(faces(j));
    }

    for (int i = 0; i < _typesCells.at("top").size(); i++) {
        auto iCell = _typesCells.at("top")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            if (_neighborsCells.at(faces(j)).size() == 1 and _facesAxis[faces(j)] == 2)
                top.push_back(faces(j));
    }


    for (int i = 0; i < _typesCells.at("nonbound").size(); i++) {
        auto iCell = _typesCells.at("nonbound")(i);
        auto &faces = _neighborsFaces.at(iCell);
        for (int j = 0; j < faces.size(); j++)
            nonboundSet.insert(faces(j));
    }


    copyStdVectotToEigenVector<int>(left, _typesFaces.at("left"));
    copyStdVectotToEigenVector<int>(right, _typesFaces.at("right"));
    copyStdVectotToEigenVector<int>(front, _typesFaces.at("front"));
    copyStdVectotToEigenVector<int>(back, _typesFaces.at("back"));
    copyStdVectotToEigenVector<int>(bottom, _typesFaces.at("bottom"));
    copyStdVectotToEigenVector<int>(top, _typesFaces.at("top"));

    std::vector<int> nonbound(nonboundSet.size());
    std::copy(nonboundSet.begin(), nonboundSet.end(), nonbound.begin());
    copyStdVectotToEigenVector<int>(nonbound, _typesFaces.at("nonbound"));

}


void Sgrid::calculateGridProps() {

    calculatePointsN();

    calculateCellsDims();
    calculateCellsN();

    calculateFacesDims();
    calculateFacesN();
    calculateFacesAxis();

    calculateCellV();
    calculateFaceS();

    calculateNeighborsFaces();
    calculateNeighborsCells();

    calculateNormalsNeighborsCells();
    calculateNormalsNeighborsFaces();

    calculateMainTypesCells();
    calculateMainTypesFaces();

}
