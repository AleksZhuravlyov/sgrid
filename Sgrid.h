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
#include <iostream>

#include <Eigen/Dense>

class Sgrid {

public:

    Sgrid(Eigen::Ref<Eigen::Vector3i> pointsDims,
          Eigen::Ref<Eigen::Vector3d> pointsOrigin,
          Eigen::Ref<Eigen::Vector3d> spacing);

    Sgrid(const std::string &fileName);

    virtual ~Sgrid() = default;

    /// Main methods

    void save(const std::string &fileName);


    /// Accessors and mutators

    Eigen::Ref<Eigen::Vector3i> getPointsDims();

    void setPointsDims(Eigen::Ref<Eigen::Vector3i> pointsDims);


    Eigen::Ref<Eigen::Vector3d> getPointsOrigin();

    void setPointsOrigin(Eigen::Ref<Eigen::Vector3d> pointsOrigin);


    Eigen::Ref<Eigen::Vector3d> getSpacing();

    void setSpacing(Eigen::Ref<Eigen::Vector3d> spacing);


    Eigen::Ref<Eigen::Vector3i> getCellsDims();

    void setCellsDims(Eigen::Ref<Eigen::Vector3i> cellsDims);


    std::map<std::string, Eigen::Ref<Eigen::VectorXi>> getActiveCells();

    void setActiveCells(
            std::map<std::string, Eigen::Ref<Eigen::VectorXi>> activeCells);


    std::map<int, Eigen::Ref<Eigen::Vector3i>> getFacesDims();

    void setFacesDims(std::map<int, Eigen::Ref<Eigen::Vector3i>> facesDims);


    std::map<std::string, Eigen::Ref<Eigen::VectorXi>> getBoundFaces();

    void setBoundFaces(
            std::map<std::string, Eigen::Ref<Eigen::VectorXi>> boundFaces);


    Eigen::Ref<Eigen::Vector3d> getFaceS();

    void setFaceS(Eigen::Ref<Eigen::Vector3d> faceS);


    std::map<int, Eigen::Ref<Eigen::VectorXi>> getNeighborsFaces();

    void setNeighborsFaces(
            std::map<int, Eigen::Ref<Eigen::VectorXi>> neighborsFaces);

    std::map<int, Eigen::Ref<Eigen::VectorXi>> getNeighborsCells();

    void setNeighborsCells(
            std::map<int, Eigen::Ref<Eigen::VectorXi>> neighborsCells);


    std::map<int, Eigen::Ref<Eigen::VectorXi>> getNormalsNeighborsCells();

    void setNormalsNeighborsCells(
            std::map<int, Eigen::Ref<Eigen::VectorXi>> normalsNeighborsCells);

    std::map<int, Eigen::Ref<Eigen::VectorXi>> getNormalsNeighborsFaces();

    void setNormalsNeighborsFaces(
            std::map<int, Eigen::Ref<Eigen::VectorXi>> normalsNeighborsFaces);


    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getPointsArrays();

    void setPointsArrays(
            std::map<std::string, Eigen::Ref<Eigen::VectorXd>> pointsArrays);

    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getCellsArrays();

    void setCellsArrays(
            std::map<std::string, Eigen::Ref<Eigen::VectorXd>> cellsArrays);

    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getFacesArrays();

    void setFacesArrays(
            std::map<std::string, Eigen::Ref<Eigen::VectorXd>> facesArrays);


    Eigen::Map<Eigen::Vector3i> _pointsDims;
    Eigen::Map<Eigen::Vector3d> _spacing;
    Eigen::Map<Eigen::Vector3d> _pointsOrigin;

    int _pointsN;

    Eigen::Map<Eigen::Vector3i> _cellsDims;
    int _cellsN;
    std::map<std::string, Eigen::Map<Eigen::VectorXi>> _activeCells;


    std::map<int, Eigen::Map<Eigen::Vector3i>> _facesDims;
    int _facesN;
    std::map<int, int> _facesAxis;
    std::map<std::string, Eigen::Map<Eigen::VectorXi>> _boundFaces;

    double _cellV;
    Eigen::Map<Eigen::Vector3d> _faceS;

    std::map<int, Eigen::Map<Eigen::VectorXi>> _neighborsFaces;
    std::map<int, Eigen::Map<Eigen::VectorXi>> _neighborsCells;

    std::map<int, Eigen::Map<Eigen::VectorXi>> _normalsNeighborsCells;
    std::map<int, Eigen::Map<Eigen::VectorXi>> _normalsNeighborsFaces;


    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _pointsArrays;
    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _cellsArrays;
    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _facesArrays;


    /// Accessory shitty constructor methods

    void calculatePointsN();

    void calculateCellsDims();

    void calculateCellsN();

    void calculateFacesDims();

    void calculateFacesN();

    void calculateFacesAxis();

    void calculateCellV();

    void calculateFaceS();

    void calculateNeighborsFaces();

    void calculateNeighborsCells();

    void calculateNormalsNeighborsCells();

    void calculateNormalsNeighborsFaces();


    void calculateGridProps();


};


void saveFilesCollectionToFile(const std::string &fileName,
                               const std::vector<std::string> &filesNames,
                               const std::vector<std::string> &filesDescriptions);

#endif // SGRID_H