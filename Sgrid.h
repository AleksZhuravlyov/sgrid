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

class Sgrid {

public:

    Sgrid(const std::vector<int> &pointsDims,
          const std::vector<double> &pointsOrigin,
          const std::vector<double> &spacing);

    Sgrid(const std::string &fileName);

    virtual ~Sgrid() = default;


    // Users methods

    void save(const std::string &fileName,
              const bool &isStructured);


    int calculateICell(const int &iX, const int &iY, const int &iZ);

    int calculateIXCell(const int &iCell);

    int calculateIYCell(const int &iCell);

    int calculateIZCell(const int &iCell);


    int calculateAxisFace(const int &iFace);

    int calculateOffsetFace(const int &axis);

    int calculateIFace(const int &axis, const int &iX, const int &iY, const int &iZ);

    int calculateIXFace(const int &iFace);

    int calculateIYFace(const int &iFace);

    int calculateIZFace(const int &iFace);


    void setCellsType(const std::string &name, Eigen::Ref<Eigen::VectorXi> cells);

    void processTypesByCellsType(const std::string &name);


    // Accessory shitty constructor methods

    void calculateFacesDims();

    void calculateFacesNs();

    void calculateFacesN();

    void calculateFacesAxes();

    void calculateFaceSs();

    void calculateNeighborsFaces();

    void calculateNeighborsCells();

    void calculateNormalsNeighborsCells();

    void calculateNormalsNeighborsFaces();

    void calculateMainTypesCells();

    void calculateMainTypesFaces();


    void calculateGridProps();


    // Eigens accessors and mutators

    std::map<std::string, Eigen::Ref<Eigen::VectorXi>> getTypesCells();

    void setTypesCells(std::map<std::string, Eigen::Ref<Eigen::VectorXi>> typesCells);


    std::map<std::string, Eigen::Ref<Eigen::VectorXi>> getTypesFaces();

    void setTypesFaces(std::map<std::string, Eigen::Ref<Eigen::VectorXi>> typesFaces);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getPointsArrays();

    void setPointsArrays(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> pointsArrays);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getCellsArrays();

    void setCellsArrays(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> cellsArrays);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getFacesArrays();

    void setFacesArrays(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> facesArrays);


    /// Attributes

    std::vector<int> _pointsDims;
    int _pointsN;
    std::vector<double> _pointsOrigin;
    std::vector<double> _spacing;

    std::vector<int> _cellsDims;
    int _cellsN;
    double _cellV;
    std::map<std::string, Eigen::Map<Eigen::VectorXi>> _typesCells;

    std::vector<std::vector<int>> _facesDimss;
    std::vector<int> _facesNs;
    int _facesN;
    std::vector<double> _facesSs;
    std::vector<int> _facesAxes;
    std::map<std::string, Eigen::Map<Eigen::VectorXi>> _typesFaces;

    std::map<int, std::vector<int>> _neighborsFaces;
    std::map<int, std::vector<int>> _neighborsCells;

    std::map<int, std::vector<int>> _normalsNeighborsCells;
    std::map<int, std::vector<int>> _normalsNeighborsFaces;

    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _pointsArrays;
    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _cellsArrays;
    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _facesArrays;


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