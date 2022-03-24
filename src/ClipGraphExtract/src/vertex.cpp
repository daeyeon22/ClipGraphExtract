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


void Vertex::setMinSlack(double slack) {
    minSlack_ = slack;
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

void Vertex::setWhiteSpaceL(int space) {
    whiteSpaceL_ = space;
}

void Vertex::setWhiteSpaceR(int space) {
    whiteSpaceR_ = space;
}

int Vertex::getSize() {
    dbBox* box = inst_->getBBox();
    int width = box->xMax() - box->xMin();
    int height = box->yMax() - box->yMin();
    return width*height;
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

int Vertex::getWhiteSpaceL() {
    return whiteSpaceL_;
}

int Vertex::getWhiteSpaceR() {
    return whiteSpaceR_;
}

double Vertex::getMinSlack() {
    return minSlack_;
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

