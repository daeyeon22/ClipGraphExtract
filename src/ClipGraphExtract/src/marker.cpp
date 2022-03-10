#include "grid.h"
#include <iostream>


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

void Marker::setFromNet(RSMT* rsmt) {
    fromNet_ = rsmt;    
}

void Marker::setToNet(RSMT* rsmt) {
    toNet_ = rsmt;
}

void Marker::setToInst(dbInst* inst) {
    toInst_ = inst;
}

Marker::Tag Marker::getTag() {
    if(fromNet_ != NULL && toNet_ != NULL)
        return Tag::N2N;
    else if(fromNet_ != NULL && toInst_ != NULL)
        return Tag::N2I;
    else
        return Tag::N2I;
}


string Marker::getType() {
    return type_;
}

string Marker::getRule() {
    return rule_;
}





RSMT* Marker::getFromNet() {
    return fromNet_;
}

RSMT* Marker::getToNet() {
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



Marker::Category Marker::getCategory() {
    if(tag_ == Tag::N2N) {
        if(fromNet_->isLocalNet() && toNet_->isLocalNet()) {
            return Category::L2L;
        } else if(!fromNet_->isLocalNet() && toNet_->isLocalNet()) {
            return Category::L2G;
        } else if(fromNet_->isLocalNet() && !toNet_->isLocalNet()) {
            return Category::L2G;
        } else {
            return Category::G2G;
        } 
    } else {
        //if (tag_ == Tag::N2I) {
        if(fromNet_->isLocalNet()) {
            return Category::L2I;
        } else {
            return Category::G2I;
        }
    }
}




void Marker::print() {
    Category ctgy = getCategory();
    switch(ctgy) {
        case Category::L2L:
            cout << "Marker from Local to Local" << endl; break;
        case Category::L2G:
            cout << "Marker from Local to Global" << endl; break;
        case Category::L2I:
            cout << "Marker from Local to Instance" << endl; break;
        case Category::G2G:
            cout << "Marker from Global to Global" << endl; break;
        case Category::G2I:
            cout << "Marker from Global to Instance" << endl; break;
        default:
            cout << "Exception case..." << endl; break;
    }
}







};




