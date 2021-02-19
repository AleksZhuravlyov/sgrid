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
#ifndef SGRID_H
#define SGRID_H

#include<string>
#include<vector>
#include<map>
#include<set>
#include <iostream>

#include <Eigen/Dense>


namespace Eigen {
    typedef Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> VectorXui64;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;
}

class Sgrid {

public:

    Sgrid(const std::vector<uint16_t> &pointsDims,
          const std::vector<double> &pointsOrigin,
          const std::vector<double> &spacing);

    virtual ~Sgrid() = default;


    // Users methods

    void save(const std::string &fileName, const std::string &type);

    void saveCells(const std::string &fileName);

    void saveFaces(const std::string &fileName);


    uint64_t calculateINode(const uint16_t &iX, const uint16_t &iY, const uint16_t &iZ);


    uint64_t calculateICell(const uint16_t &iX, const uint16_t &iY, const uint16_t &iZ);

    uint16_t calculateIXCell(const uint64_t &iCell);

    uint16_t calculateIYCell(const uint64_t &iCell);

    uint16_t calculateIZCell(const uint64_t &iCell);


    uint8_t calculateAxisFace(const uint64_t &iFace);

    uint64_t calculateOffsetFace(const uint8_t &axis);

    uint64_t calculateIFace(const uint8_t &axis, const uint16_t &iX, const uint16_t &iY,
                            const uint16_t &iZ);

    uint16_t calculateIXFace(const uint64_t &iFace);

    uint64_t calculateIYFace(const uint16_t &iFace);

    uint64_t calculateIZFace(const uint16_t &iFace);


    void setCellsType(const std::string &name, Eigen::Ref<Eigen::VectorXui64> cells);

    void processTypesByCellsType(const std::string &name);


    // Accessory shitty constructor methods

    void calculateNodesCoordinates();

    void calculatePointsNodes();

    void calculateCellsNodes();

    void calculateFacesDimss();

    void calculateFacesNs();

    void calculateFacesN();

    void calculateFacesAxes();

    void calculateFaceSs();

    void calculateFacesNodes();

    void calculateNeighborsFaces();

    void calculateNeighborsCells();

    void calculateNormalsNeighborsCells();

    void calculateNormalsNeighborsFaces();

    void calculateMainTypesCells();

    void calculateMainTypesFaces();


    void calculateGridProps();


    // Eigens accessors and mutators

    std::map<std::string, Eigen::Ref<Eigen::VectorXui64>> getTypesCells();

    void setTypesCells(std::map<std::string, Eigen::Ref<Eigen::VectorXui64>> typesCells);


    std::map<std::string, Eigen::Ref<Eigen::VectorXui64>> getTypesFaces();

    void setTypesFaces(std::map<std::string, Eigen::Ref<Eigen::VectorXui64>> typesFaces);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getPointsArrays();

    void setPointsArrays(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> pointsArrays);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getCellsArrays();

    void setCellsArrays(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> cellsArrays);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getFacesArrays();

    void setFacesArrays(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> facesArrays);


    std::map<std::string, Eigen::Ref<Eigen::VectorXb>> getCellsConditions();

    void setCellsConditions(std::map<std::string, Eigen::Ref<Eigen::VectorXb>> cellsConditions);


    /// Attributes

    std::vector<uint16_t> _pointsDims;
    uint64_t _pointsN;
    std::vector<double> _pointsOrigin;
    std::vector<double> _spacing;

    std::vector<double> _nodesCoordinates;

    std::map<uint64_t, std::vector<uint64_t>> _pointsNodes;

    std::vector<uint16_t> _cellsDims;
    uint64_t _cellsN;
    std::map<uint64_t, std::vector<uint64_t>> _cellsNodes;
    double _cellV;
    std::map<std::string, Eigen::Map<Eigen::VectorXui64>> _typesCells;

    std::vector<std::vector<uint16_t>> _facesDimss;
    std::vector<uint64_t> _facesNs;
    uint64_t _facesN;
    std::map<uint64_t, std::vector<uint64_t>> _facesNodes;
    std::vector<double> _facesSs;
    std::vector<uint8_t> _facesAxes;
    std::map<std::string, Eigen::Map<Eigen::VectorXui64>> _typesFaces;

    std::map<uint64_t, std::vector<uint64_t>> _neighborsFaces;
    std::map<uint64_t, std::vector<uint64_t>> _neighborsCells;

    std::map<uint64_t, std::vector<int8_t>> _normalsNeighborsCells;
    std::map<uint64_t, std::vector<int8_t>> _normalsNeighborsFaces;

    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _pointsArrays;
    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _cellsArrays;
    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _facesArrays;

    std::map<std::string, Eigen::Map<Eigen::VectorXb>> _cellsConditions;

};


void saveFilesCollectionToFile(const std::string &fileName,
                               const std::vector<std::string> &filesNames,
                               const std::vector<std::string> &filesDescriptions);

template<class T>
void copyStdVectotToEigenVector(
        std::vector<T> &stdVector,
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> &eigenVector) {
    auto size = stdVector.size();
    delete[] eigenVector.data();
    new(&eigenVector) Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(new T[size], size);
    for (int i = 0; i < size; i++)
        eigenVector(i) = stdVector[i];
}

template<class T>
void copyStdSetToEigenVector(
        std::set<T> &stdSet,
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> &eigenVector) {
    auto size = stdSet.size();
    delete[] eigenVector.data();
    new(&eigenVector) Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(new T[size], size);
    int i = 0;
    for (auto &value: stdSet) {
        eigenVector(i) = value;
        i++;
    }
}


#endif // SGRID_H