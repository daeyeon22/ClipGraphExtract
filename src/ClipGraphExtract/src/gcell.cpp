#include "grid.h"
#include <iostream>


namespace feature_extractor {

using namespace std;
using namespace odb;

Gcell::Gcell() :
    numInstances_(0), numTerminals_(0), numLocalNets_(0),
    numGlobalNets_(0), totalCellArea_(0), totalPinArea_(0),
    cellDensity_(0), pinDensity_(0), RUDY_(0) {}


void Gcell::setBoundary(Rect rect) { 
    bbox_ = rect;
}

void Gcell::setTrackSupply(int tSup) {
    rmEGR.setTrackSupply(tSup);
    rmRSMT.setTrackSupply(tSup);
    rmDR.setTrackSupply(tSup);
}

void Gcell::setWireCapacity(int wCap) {
    rmEGR.setWireCapacity(wCap);
    rmRSMT.setTrackSupply(wCap);
    rmDR.setTrackSupply(wCap);
}



bgBox Gcell::getQueryBox() {
    return bgBox(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
}

uint Gcell::getArea() {
    return bbox_.area();    
}

double Gcell::getCellDensity() {
    return cellDensity_;
}

double Gcell::getPinDensity() {
    return pinDensity_;
}

double Gcell::getRUDY() {
    return RUDY_;
}

void
Gcell::extractFeaturePL(BoxRtree<odb::dbInst*> &rtree) {
    vector<pair<bgBox, dbInst*>> queryResults;
    // Query
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));
    for(auto& val : queryResults) {
        dbInst* inst = val.second;
        insts_.push_back(inst);
        odb::Rect instBBox;
        inst->getBBox()->getBox(instBBox);

        // Cell Area
        odb::Rect overlap = instBBox.intersect(bbox_);
        int cellArea = overlap.dx() * overlap.dy();

        // Termianls
        for(dbITerm* iterm : inst->getITerms()) {

            int pinArea = 0;
            dbMTerm* master = iterm->getMTerm();
            dbSet<dbMPin> mPins = master->getMPins();
            // Pin area 
            for(dbMPin* mPin : master->getMPins()) {
                for(dbBox* pBox : mPin->getGeometry()) {
                    Rect pinBBox;
                    pBox->getBox(pinBBox);
                    odb::Rect pinOverlap = pinBBox.intersect(bbox_);
                    pinArea += pinOverlap.dx() * pinOverlap.dy();
                }
            }
            
            if(iterm->getNet() == NULL) {
                // unconnected pin
            } else {
                // connected pin
            }
            numTerminals_++;
            totalPinArea_ += pinArea;
        }
        numInstances_++;
        totalCellArea_ += cellArea; 
    }
    //
cellDensity_ = 1.0* totalCellArea_ / getArea();
    pinDensity_ = 1.0 * totalPinArea_ / getArea();

}


void
Gcell::extractFeatureEGR(SegRtree<odb::dbNet*> &rtree) {

    vector<pair<bgSeg, dbNet*>> queryResults;

    // Query
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));

    bgSeg lb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMin(), bbox_.yMax()));
    bgSeg rb(bgPoint(bbox_.xMax(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg tb(bgPoint(bbox_.xMin(), bbox_.yMax()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg bb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMin()));

    for(auto &val : queryResults) {
        bgSeg wire_seg = val.first;
        dbNet* net = val.second;
        
        if(bg::intersects(lb, wire_seg)) {
            rmEGR.trackDemand[Orient::LEFT]++;
        }
        if(bg::intersects(rb, wire_seg)) {
            rmEGR.trackDemand[Orient::RIGHT]++;
        }
        if(bg::intersects(tb, wire_seg)) {
            rmEGR.trackDemand[Orient::TOP]++;
        }
        if(bg::intersects(bb, wire_seg)) {
            rmEGR.trackDemand[Orient::BOTTOM]++;
        }

        int x0 = bg::get<0,0>(wire_seg);
        int y0 = bg::get<0,1>(wire_seg);
        int x1 = bg::get<1,0>(wire_seg);
        int y1 = bg::get<1,1>(wire_seg);
        x0 = max(x0, bbox_.xMin());
        y0 = max(y0, bbox_.yMin());
        x1 = min(x1, bbox_.xMax());
        y1 = min(y1, bbox_.yMax());
        // intersection wirelength
        rmEGR.wireLength += (x1-x0) + (y1-y0);
    }
}


void
Gcell::extractFeatureRSMT(SegRtree<RSMT*> &rtree) {
    vector<pair<bgSeg, RSMT*>> queryResults;

    // Query
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));
    bgSeg lb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMin(), bbox_.yMax()));
    bgSeg rb(bgPoint(bbox_.xMax(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg tb(bgPoint(bbox_.xMin(), bbox_.yMax()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg bb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMin()));

    set<RSMT*> RSMTs;

    for(auto &val : queryResults) {
        bgSeg wire_seg = val.first;
        RSMT* myRSMT = val.second;
        
        if(bg::intersects(lb, wire_seg)) {
            rmRSMT.trackDemand[Orient::LEFT]++;
        }
        if(bg::intersects(rb, wire_seg)) {
            rmRSMT.trackDemand[Orient::RIGHT]++;
        }
        if(bg::intersects(tb, wire_seg)) {
            rmRSMT.trackDemand[Orient::TOP]++;
        }
        if(bg::intersects(bb, wire_seg)) {
            rmRSMT.trackDemand[Orient::BOTTOM]++;
        }

        int x0 = bg::get<0,0>(wire_seg);
        int y0 = bg::get<0,1>(wire_seg);
        int x1 = bg::get<1,0>(wire_seg);
        int y1 = bg::get<1,1>(wire_seg);
        x0 = max(x0, bbox_.xMin());
        y0 = max(y0, bbox_.yMin());
        x1 = min(x1, bbox_.xMax());
        y1 = min(y1, bbox_.yMax());
        // intersection wirelength
        rmRSMT.wireLength += (x1-x0) + (y1-y0);
        RSMTs.insert(myRSMT);
    }

    cout << "GCELL (" << bbox_.xMin() << " " << bbox_.yMin() << ") (" << bbox_.xMax() << " " << bbox_.yMax() << ")" << endl;

    for(RSMT* myRSMT : RSMTs) {
        // update RUDY
        odb::Rect r1 = bbox_;
        odb::Rect r2 = myRSMT->getBBox();
        uint area1 = r2.intersect(r1).area();
        uint area2 = r1.area();
        double dn = myRSMT->getWireUniformDensity();
        double R = 1.0* area1 / area2;
        //double R = 1.0 * r2.intersect(r1).area() / r1.area();
        double partialRUDY = dn * R;   
       
        cout << "   - "<< myRSMT->getNet()->getName() << " RUDY = " << partialRUDY << " (" << dn << " " << R << ")" << endl;
        RUDY_ += partialRUDY;
    }
}



}
