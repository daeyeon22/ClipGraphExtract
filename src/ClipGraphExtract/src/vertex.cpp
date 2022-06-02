#include <vector>
#include <string>
#include "opendb/db.h"
#include "instGraph.h"

namespace ClipGraphExtract {

using namespace odb;
using namespace std;


odb::dbInst* Vertex::inst() const {
  return inst_;
}

const std::vector<Edge*> & Vertex::inEdges() const {
  return inEdges_;
}

const std::vector<Edge*> & Vertex::outEdges() const {
  return outEdges_;
}

float Vertex::weight() const {
  return weight_;
}

void Vertex::setId(int id) {
  id_ = id;
}

int Vertex::id() const {
  return id_;
}

Vertex::Vertex() : 
  inst_(nullptr), 
  weight_(0),
  id_(0) {};

Vertex::Vertex(odb::dbInst* inst, 
    int id, float weight) 
  : Vertex() {
  inst_ = inst;
  id_ = id;
  weight_ = weight;
    absSlack_ = 0.0;
    relSlack_ = 0.0;
    numAccPoints_ = 0;
    numBlkPoints_ = 0;
    numBndPoints_ = 0;
    whiteSpaceL_ = 0;
    whiteSpaceR_ = 0;
    whiteSpaceD_ = 0;
    whiteSpaceT_ = 0;
    relPosX_=0;
    relPosY_=0;
    isCrit_=false;
    cutEdges_=0;
    bboxSize_=0.0;
    sWireOverlap_=0.0;
    numDrvs_=0;
} 


void Vertex::setInst(odb::dbInst* inst) {
  inst_ = inst;
}

void Vertex::setWeight(float weight) {
  weight_ = weight;
}

void Vertex::addInEdge(Edge* edge) {
  inEdges_.push_back(edge);
}

void Vertex::addOutEdge(Edge* edge) { 
  outEdges_.push_back(edge);
}

void Vertex::setNumDrvs(int numDrvs) {
    numDrvs_ = numDrvs;
}


void Vertex::setBBoxSize(double bboxSize) {
    bboxSize_ = bboxSize;
}

void Vertex::setCellType(int cellType) {
    cellType_ = cellType;
}

void Vertex::setIsCrit(bool isCrit) {
    isCrit_ = isCrit;
}

void Vertex::setSWireOverlap(double overlap) {
    sWireOverlap_ = overlap;
}

void Vertex::setCutEdges(int cutEdges) {
    cutEdges_= cutEdges;
}

void Vertex::setSlack(double clockPeriod, double slack) {

    absSlack_ = slack;
    relSlack_ = slack / clockPeriod;
}

void Vertex::setNumAccPoints(int numPoints) {
    numAccPoints_ = numPoints;
}

void Vertex::setNumBlkPoints(int numPoints) {
    numBlkPoints_ = numPoints;
}

void Vertex::setNumBndPoints(int numPoints) {
    numBndPoints_ = numPoints;
}

void Vertex::setViaFeature(double powerViaDistance) {
    powerViaDistance_ = powerViaDistance;
}

void Vertex::setWhiteSpaceL(double space) {
    whiteSpaceL_ = space;
}

void Vertex::setWhiteSpaceR(double space) {
    whiteSpaceR_ = space;
}

void Vertex::setWhiteSpaceT(double space) {
    whiteSpaceT_ = space;
}

void Vertex::setWhiteSpaceD(double space) {
    whiteSpaceD_ = space;
}


void Vertex::setRelPos(double posX, double posY) {
    relPosX_ = posX, relPosY_ = posY;
}

double Vertex::getRelPosX() {
    return relPosX_;
}

double Vertex::getRelPosY() {
    return relPosY_;
}

int Vertex::getNumDrvs() {
    return numDrvs_;
}



int Vertex::getSize() {
    dbBox* box = inst_->getBBox();
    int width = box->xMax() - box->xMin();
    int height = box->yMax() - box->yMin();
    return width*height;
}


int Vertex::getCutEdges() {
    return cutEdges_;
}

int Vertex::getCellType() {
    return cellType_;
}

double Vertex::getSWireOverlap() {
    return sWireOverlap_;
}

double Vertex::getBBoxSize() {
    return bboxSize_;
}

bool Vertex::getIsCrit() {
    return isCrit_;
}


int Vertex::getDegree() {
    return getNumInEdges() + getNumOutEdges();
}

int Vertex::getNumInEdges() {
    return inEdges_.size();
}

int Vertex::getNumOutEdges() {
    return outEdges_.size();
}

bool Vertex::isClocked() {
    if(inst_->getClockedTerm() != NULL) {
        return true;
    } else {
        return false;
    }
}

double Vertex::getWhiteSpaceL() {
    return whiteSpaceL_;
}

double Vertex::getWhiteSpaceR() {
    return whiteSpaceR_;
}

double Vertex::getWhiteSpaceT() {
    return whiteSpaceT_;
}

double Vertex::getWhiteSpaceD() {
    return whiteSpaceD_;
}


double Vertex::getAbsSlack() {
    return absSlack_;
}

double Vertex::getRelSlack() {
    return relSlack_;
}

int Vertex::getNumAccPoints() {
    return numAccPoints_;
}

int Vertex::getNumBlkPoints() {
    return numBlkPoints_;
}

int Vertex::getNumBndPoints() {
    return numBndPoints_;
}



};

