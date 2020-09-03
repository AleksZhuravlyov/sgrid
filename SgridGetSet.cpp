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

namespace Eigen {

    template<class T>
    void setter(Eigen::Ref<T> &source, Eigen::Map<T> &sink) {
        new(&sink) Eigen::Map<T>(source.data(), source.size());
    };

    template<class K, class V>
    std::map<K, Eigen::Ref<V>> mapGetter(std::map<K, Eigen::Map<V>> &source) {

        std::map<K, Eigen::Ref<V>> sink;
        for (auto const &ent : source)
            sink.insert(
                    std::pair<K, Eigen::Ref<V>>(ent.first,
                                                source.at(ent.first)));

        return sink;
    }

    template<class K, class V>
    void mapSetter(std::map<K, Eigen::Ref<V>> &source,
                   std::map<K, Eigen::Map<V>> &sink) {
        sink.clear();
        for (auto const &ent : source)
            sink.insert(
                    std::pair<K, Eigen::Map<V>>(ent.first, Eigen::Map<V>(
                            source.at(ent.first).data(),
                            source.at(ent.first).size())));
    }
}


Eigen::Ref<Eigen::Vector3i> Sgrid::getPointsDims() {
    return _pointsDims;
}

void Sgrid::setPointsDims(Eigen::Ref<Eigen::Vector3i> pointsDims) {
    Eigen::setter<Eigen::Vector3i>(pointsDims, _pointsDims);
}

Eigen::Ref<Eigen::Vector3d> Sgrid::getPointsOrigin() {
    return _pointsOrigin;
}

void Sgrid::setPointsOrigin(Eigen::Ref<Eigen::Vector3d> pointsOrigin) {
    Eigen::setter<Eigen::Vector3d>(pointsOrigin, _pointsOrigin);
}


Eigen::Ref<Eigen::Vector3d> Sgrid::getSpacing() {
    return _spacing;
}

void Sgrid::setSpacing(Eigen::Ref<Eigen::Vector3d> spacing) {
    Eigen::setter<Eigen::Vector3d>(spacing, _spacing);
}


Eigen::Ref<Eigen::Vector3i> Sgrid::getCellsDims() {
    return _cellsDims;
}

void Sgrid::setCellsDims(Eigen::Ref<Eigen::Vector3i> cellsDims) {
    Eigen::setter<Eigen::Vector3i>(cellsDims, _cellsDims);
}


std::map<int, Eigen::Ref<Eigen::Vector3i>> Sgrid::getFacesDims() {
    return Eigen::mapGetter<int, Eigen::Vector3i>(_facesDims);
}

void Sgrid::setFacesDims(std::map<int, Eigen::Ref<Eigen::Vector3i>> facesDims) {
    Eigen::mapSetter<int, Eigen::Vector3i>(facesDims, _facesDims);
}


Eigen::Ref<Eigen::Vector3d> Sgrid::getFaceS() {
    return _faceS;
}

void Sgrid::setFaceS(Eigen::Ref<Eigen::Vector3d> faceS) {
    Eigen::setter<Eigen::Vector3d>(faceS, _faceS);
}


std::map<int, Eigen::Ref<Eigen::VectorXi>> Sgrid::getNeighborsFaces() {
    return Eigen::mapGetter<int, Eigen::VectorXi>(_neighborsFaces);
}

void Sgrid::setNeighborsFaces(
        std::map<int, Eigen::Ref<Eigen::VectorXi>> neighborsFaces) {
    Eigen::mapSetter<int, Eigen::VectorXi>(neighborsFaces, _neighborsFaces);
}

std::map<int, Eigen::Ref<Eigen::VectorXi>> Sgrid::getNeighborsCells() {
    return Eigen::mapGetter<int, Eigen::VectorXi>(_neighborsCells);
}

void Sgrid::setNeighborsCells(
        std::map<int, Eigen::Ref<Eigen::VectorXi>> neighborsCells) {
    Eigen::mapSetter<int, Eigen::VectorXi>(neighborsCells, _neighborsCells);
}


std::map<int, Eigen::Ref<Eigen::VectorXi>> Sgrid::getNormalsNeighborsCells() {
    return Eigen::mapGetter<int, Eigen::VectorXi>(_normalsNeighborsCells);
}

void Sgrid::setNormalsNeighborsCells(
        std::map<int, Eigen::Ref<Eigen::VectorXi>> normalsNeighborsCells) {
    Eigen::mapSetter<int, Eigen::VectorXi>(normalsNeighborsCells,
                                           _normalsNeighborsCells);
}

std::map<int, Eigen::Ref<Eigen::VectorXi>> Sgrid::getNormalsNeighborsFaces() {
    return Eigen::mapGetter<int, Eigen::VectorXi>(_normalsNeighborsFaces);
}

void Sgrid::setNormalsNeighborsFaces(
        std::map<int, Eigen::Ref<Eigen::VectorXi>> normalsNeighborsFaces) {
    Eigen::mapSetter<int, Eigen::VectorXi>(normalsNeighborsFaces,
                                           _normalsNeighborsFaces);
}


std::map<std::string, Eigen::Ref<Eigen::VectorXd>> Sgrid::getPointsArrays() {
    return Eigen::mapGetter<std::string, Eigen::VectorXd>(_pointsArrays);
}

void Sgrid::setPointsArrays(
        std::map<std::string, Eigen::Ref<Eigen::VectorXd>> pointsArrays) {
    Eigen::mapSetter<std::string, Eigen::VectorXd>(pointsArrays, _pointsArrays);
}

std::map<std::string, Eigen::Ref<Eigen::VectorXd>> Sgrid::getCellsArrays() {
    return Eigen::mapGetter<std::string, Eigen::VectorXd>(_cellsArrays);
}

void Sgrid::setCellsArrays(
        std::map<std::string, Eigen::Ref<Eigen::VectorXd>> cellsArrays) {
    Eigen::mapSetter<std::string, Eigen::VectorXd>(cellsArrays, _cellsArrays);
}

std::map<std::string, Eigen::Ref<Eigen::VectorXd>> Sgrid::getFacesArrays() {
    return Eigen::mapGetter<std::string, Eigen::VectorXd>(_facesArrays);
}

void Sgrid::setFacesArrays(
        std::map<std::string, Eigen::Ref<Eigen::VectorXd>> facesArrays) {
    Eigen::mapSetter<std::string, Eigen::VectorXd>(facesArrays, _facesArrays);
}
