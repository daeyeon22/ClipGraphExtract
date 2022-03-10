#include "grid.h"


namespace feature_extractor {

using namespace odb;
using namespace std;


Marker::Marker():
    type_(""), rule_(""), tag_(Marker::Tag::N2N), 
    fromNet_(nullptr), toNet_(nullptr), toInst_(nullptr),
    bbox_(Rect(0,0,0,0)) {}



void Marker::setType(string type) {
    type_ = type;
}

void Marker::setRule(string rule) {
    rule_ = rule;
}

void Marker::setBoundary(Rect rect) {
    bbox_ = rect;
}

void Marker::setTag(Marker::Tag tag) {
    tag_ = tag;
}

void Marker::setFromNet(dbNet* net) {
    fromNet_ = net;    
}

void Marker::setToNet(dbNet* net) {
    toNet_ = net;
}

void Marker::setToInst(dbInst* inst) {
    toInst_ = inst;
}

odb::dbNet* Marker::getFromNet() {
    return fromNet_;
}

odb::dbNet* Marker::getToNet() {
    return toNet_;
}

odb::dbInst* Marker::getToInst() {
    return toInst_;
}

bgBox Marker::getQueryBox() {
    return bgBox(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
}


Rect Marker::getBBox() {
    return bbox_;
}

Point Marker::getCentor() {
    int cx = bbox_.xMin() + bbox_.dx() / 2;
    int cy = bbox_.yMin() + bbox_.dy() / 2;
    return Point(cx, cy);
}









};




