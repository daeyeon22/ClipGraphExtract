#include "grid.h"
#include <iostream>


namespace ClipGraphExtract {

using namespace odb;
using namespace std;


Marker::Marker():
    type_(""), rule_(""), 
    fromTag_(Marker::Tag::NONE), 
    toTag_(Marker::Tag::NONE), 
    fromNet_(nullptr), toNet_(nullptr), 
    fromInst_(nullptr), toInst_(nullptr),
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

void Marker::setFromTag(Marker::Tag tag) {
    fromTag_ = tag;
}

void Marker::setToTag(Marker::Tag tag) {
    toTag_ = tag;
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

void Marker::setFromInst(dbInst* inst) {
    fromInst_ = inst;
}

void Marker::setLayer(dbTechLayer* layer) {
    layer_ = layer;
}



Marker::Tag Marker::getFromTag() {
    return fromTag_;
}

Marker::Tag Marker::getToTag() {
    return toTag_;
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

odb::dbInst* Marker::getFromInst() {
    return fromInst_;
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



bool Marker::isFromNet() {
    return (fromTag_ == Tag::RWoN)? true :  false;
}

bool Marker::isToNet() {
    return (toTag_ == Tag::RWoN) ? true : false;
}

bool Marker::isToInst() {
    return (toTag_ == Tag::BoC || toTag_ == Tag::PoC) ? true : false;
}

bool Marker::isFromInst() {
    return (fromTag_ == Tag::BoC || fromTag_ == Tag::PoC) ? true : false;
}





Marker::Category Marker::getCategory() {
    if(isFromNet() && isToNet()) {
        //if(toNet_ == NULL)
        //    cout << "tonet is nullptr" << endl;
        //if(fromNet_ == NULL) 
        //    cout << "fromnet is nullptr" << endl;
       
        //cout << fromTag_ << " " << toTag_ << endl;

        if(toNet_->isLocalNet() && fromNet_->isLocalNet()) {
            return Category::L2L;
        } else if (!toNet_->isLocalNet() && fromNet_->isLocalNet()) {
            return Category::L2G;
        } else if (toNet_->isLocalNet() && !fromNet_->isLocalNet()) {
            return Category::L2G;
        } else if (!toNet_->isLocalNet() && !fromNet_->isLocalNet()) {
            return Category::G2G;
        } else {
            return Category::ERR;
        }
    } else if (isFromInst() && isToNet()) {
		if(toNet_->isLocalNet())
			return Category::L2I;
		else
			return Category::G2I;
    } else if (isFromNet() && isToInst()) {
        if(fromNet_->isLocalNet())
            return Category::L2I;
        else
            return Category::G2I;
    } else if (isFromInst() && isToInst()) {
        return Category::I2I;
    } else {
        if(toTag_ == Tag::NONE) {
            return Category::SELF;
        } else {
            return Category::ERR;
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
        case Category::I2I:
            cout << "Marker from Instance to Instance" << endl; break;
        case Category::SELF:
            cout << "Marker from Itself" << endl; break;
        default:
            cout << "Exception case..." << endl; break;
    }
}







};




