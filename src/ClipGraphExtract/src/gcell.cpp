#include "grid.h"
#include "opendb/dbTransform.h"
#include <iostream>


namespace feature_extractor {

using namespace std;
using namespace odb;

Gcell::Gcell() :
    numInstances_(0), numTerminals_(0), 
    //numLocalNets_(0), numGlobalNets_(0), 
    totalCellArea_(0), totalPinArea_(0),
    cellDensity_(0), pinDensity_(0), RUDY_(0) {}


void Gcell::setBoundary(Rect rect) { 
    bbox_ = rect;
}

void Gcell::setTrackSupply(int tSup) {
    rmEGR_.setTrackSupply((uint)tSup);
    rmPL_.setTrackSupply((uint)tSup);
    rmDR_.setTrackSupply((uint)tSup);
}

void Gcell::setWireCapacity(int wCap) {
    rmEGR_.setWireCapacity((uint)wCap);
    rmPL_.setWireCapacity((uint)wCap);
    rmDR_.setWireCapacity((uint)wCap);
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

uint Gcell::getTrackDemand(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::DR: 
            return rmDR_.getTrackDemand(orient);
        case ModelType::EGR:
            return rmEGR_.getTrackDemand(orient);
        case ModelType::PL:
            return rmPL_.getTrackDemand(orient);
        default:
            return 0; 
    }
}

uint Gcell::getTrackSupply(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::DR: 
            return rmDR_.getTrackSupply(orient);
        case ModelType::EGR:
            return rmEGR_.getTrackSupply(orient);
        case ModelType::PL:
            return rmPL_.getTrackSupply(orient);
        default:
            return 0; 
    }
}

uint Gcell::getWireCapacity(ModelType type) {
    switch(type) {
        case ModelType::DR: 
            return rmDR_.getWireCapacity();
        case ModelType::EGR:
            return rmEGR_.getWireCapacity();
        case ModelType::PL:
            return rmPL_.getWireCapacity();
        default:
            return 0; 
    }
}

uint Gcell::getNumInstances() {
    return numInstances_;
}

uint Gcell::getNumTerminals() {
    return numTerminals_;
}

uint Gcell::getNumNets() {
    return (uint)(rsmts_.size());
}


uint Gcell::getNumLocalNets() {
    uint count =0;
    for(RSMT* rsmt : rsmts_) {
        if(rsmt->isLocalNet())
            count++;
    }

    return count;
}

uint Gcell::getNumGlobalNets() {
    uint count = getNumNets() - getNumLocalNets();
    return count;

}

uint Gcell::getNumMarkers() {
    return (uint)(markers_.size());
}


double Gcell::getWireDensity(ModelType type) {
     switch(type) {
        case ModelType::DR: 
            return rmDR_.getWireDensity();
        case ModelType::EGR:
            return rmEGR_.getWireDensity();
        case ModelType::PL:
            return rmPL_.getWireDensity();
        default:
            return 0; 
    }   
}

double Gcell::getLNetDensity() {
    uint wireCap = rmPL_.getWireCapacity();
    uint wireLen = 0;
    for(RSMT* rsmt : rsmts_) {
        if(rsmt->isLocalNet()) {
            wireLen += rsmt->getWireLengthRSMT();
        }
    }
    return 1.0 * wireLen / wireCap;
}

double Gcell::getGNetDensity() {
    //uint wireCap = rmPL_.getWireCapacity();
    //uint wireLen = rmPL_.getWireLength();
    return getWireDensity(ModelType::PL) - getLNetDensity();
}




double Gcell::getChannelDensity(ModelType type) {
    double chanDenL = getChannelDensity(Orient::LEFT, type);
    double chanDenU = getChannelDensity(Orient::TOP, type);
    double chanDenR = getChannelDensity(Orient::RIGHT, type);
    double chanDenB = getChannelDensity(Orient::BOTTOM, type);
    double chanDen = (chanDenL + chanDenU + chanDenR + chanDenB) / 4;
    return chanDen;
}
double Gcell::getChannelDensity(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::DR: 
            return rmDR_.getChannelDensity(orient);
        case ModelType::EGR:
            return rmEGR_.getChannelDensity(orient);
        case ModelType::PL:
            return rmPL_.getChannelDensity(orient);
        default:
            return 0.0; 
    }
}





void
Gcell::extractFeaturePL(BoxRtree<odb::dbInst*> &rtree) {
    vector<pair<bgBox, dbInst*>> queryResults;
    // Query
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));

    //DEBUG

    set<dbInst*> visit;

    for(auto& val : queryResults) {
        dbInst* inst = val.second;
        insts_.push_back(inst);
        odb::Rect instBBox;
        inst->getBBox()->getBox(instBBox);

        if(visit.find(inst) != visit.end()){
            cout << "????" << endl;
            exit(0);
        }
        visit.insert(inst);


        int xOrigin, yOrigin;
        inst->getOrigin(xOrigin,yOrigin);
        dbTransform transform;
        transform.setOrient( inst->getOrient() );
        transform.setOffset( Point(xOrigin, yOrigin) );

        // Cell Area
        odb::Rect overlap = instBBox.intersect(bbox_);
        uint cellArea = overlap.dx() * overlap.dy();
        //cout << inst->getName() << " (" << instBBox.xMin() << " " << instBBox.yMin() 
        //    << ") (" << instBBox.xMax() << " " << instBBox.yMax() << ")" << endl;
        // Termianls
        //cout << inst->getName() << " (" << overlap.xMin() << " " << overlap.yMin() 
        //    << ") (" << overlap.xMax() << " " << overlap.yMax() << ")" << endl;
        for(dbITerm* iterm : inst->getITerms()) {
            uint pinArea = 0;
            dbMTerm* master = iterm->getMTerm();
            dbSet<dbMPin> mPins = master->getMPins();
            // Pin area 
            for(dbMPin* mPin : master->getMPins()) {
                for(dbBox* pBox : mPin->getGeometry()) {
                    if(pBox->getTechLayer()->getType() == dbTechLayerType::ROUTING) {
                        Rect pinBBox;
                        pBox->getBox(pinBBox);
                        transform.apply( pinBBox );
                        //cout << "   - (" << pinBBox.xMin() << " " << pinBBox.yMin() 
                        //<< ") (" << pinBBox.xMax() << " " << pinBBox.yMax() << ")" << endl;
                        if(pinBBox.intersects(overlap)) {
                            odb::Rect pinOverlap = pinBBox.intersect(overlap);
                            pinArea += pinOverlap.area(); //pinOverlap.dx() * pinOverlap.dy();
                            //cout << "   - " << pBox->getTechLayer()->getName() 
                            //    << " (" << pinOverlap.xMin() << " " << pinOverlap.yMin() 
                            //<< ") (" << pinOverlap.xMax() << " " << pinOverlap.yMax() 
                            //<< ") "  << pinOverlap.dx() << " " << pinOverlap.dy() << " " << pinOverlap.area() << endl;
                        }
                        //cout << "   - (" << pBox->xMin() << " " << pBox->yMin() 
                        //    << ") (" << pBox->xMax() << " " << pBox->yMax() << ")" << endl;
                        //cout << pBox->getTechLayer()->getName() << endl;
                    }
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

        //
        //if(totalPinArea_ > getArea()) {
        //    cout << "HERE" << endl;
        //}
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
            rmEGR_.incrTrackDemand(Orient::LEFT);
        }
        if(bg::intersects(rb, wire_seg)) {
            rmEGR_.incrTrackDemand(Orient::RIGHT);
        }
        if(bg::intersects(tb, wire_seg)) {
            rmEGR_.incrTrackDemand(Orient::TOP);
        }
        if(bg::intersects(bb, wire_seg)) {
            rmEGR_.incrTrackDemand(Orient::BOTTOM);
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
        uint wl = (x1-x0) + (y1-y0);
        rmEGR_.addWireLength(wl);
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

    uint wirelength=0;
    for(auto &val : queryResults) {
        bgSeg wire_seg = val.first;
        RSMT* myRSMT = val.second;
        
        if(bg::intersects(lb, wire_seg)) {
            rmPL_.incrTrackDemand(Orient::LEFT);
        }
        if(bg::intersects(rb, wire_seg)) {
            rmPL_.incrTrackDemand(Orient::RIGHT);
        }
        if(bg::intersects(tb, wire_seg)) {
            rmPL_.incrTrackDemand(Orient::TOP);
        }
        if(bg::intersects(bb, wire_seg)) {
            rmPL_.incrTrackDemand(Orient::BOTTOM);
        }

        int x0 = bg::get<0,0>(wire_seg);
        int y0 = bg::get<0,1>(wire_seg);
        int x1 = bg::get<1,0>(wire_seg);
        int y1 = bg::get<1,1>(wire_seg);

        //if(x0 > x1 || y0 > y1) {
        //    cout << "??????????" << endl;
        //    exit(0);
        //}

        
        x0 = max(x0, bbox_.xMin());
        y0 = max(y0, bbox_.yMin());
        x1 = min(x1, bbox_.xMax());
        y1 = min(y1, bbox_.yMax());
        // intersection wirelength
        uint wl = (x1-x0) + (y1-y0);
        rmPL_.addWireLength(wl);
        RSMTs.insert(myRSMT);
    }



    //cout << "GCELL (" << bbox_.xMin() << " " << bbox_.yMin() << ") (" << bbox_.xMax() << " " << bbox_.yMax() << ")" << endl;

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
       
        //cout << "   - "<< myRSMT->getNet()->getName() << " RUDY = " << partialRUDY << " (" << dn << " " << R << ")" << endl;
        RUDY_ += partialRUDY;

        //
        rsmts_.push_back(myRSMT);
        //if(myRSMT->isLocalNet())
        //    numLocalNets_++;
        //else
        //    numGlobalNets_++;

    }
}







void Gcell::print() {

    cout << "GCELL (" << bbox_.xMin() << " " << bbox_.yMin() << ") (" << bbox_.xMax() << " " << bbox_.yMax() << ")" << endl;
    cout << "   - CellDen   : " << getCellDensity() << endl;
    cout << "   - PinDen    : " << getPinDensity() << endl;
    cout << "   - WireDen   : " << getWireDensity(ModelType::PL) 
         << "(" << getWireCapacity(ModelType::PL) << ")" <<  endl;
    cout << "   - ChaDen (u): " << getChannelDensity(Orient::TOP, ModelType::PL) 
         << " (" << getTrackDemand(Orient::TOP, ModelType::PL) << " " << getTrackSupply(Orient::TOP, ModelType::PL) << ")" << endl;
    cout << "   - ChaDen (r): " << getChannelDensity(Orient::RIGHT, ModelType::PL)
         << " (" << getTrackDemand(Orient::RIGHT, ModelType::PL) << " " << getTrackSupply(Orient::RIGHT, ModelType::PL) << ")" << endl;
    cout << "   - ChaDen (b): " << getChannelDensity(Orient::BOTTOM, ModelType::PL)
         << " (" << getTrackDemand(Orient::BOTTOM, ModelType::PL) << " " << getTrackSupply(Orient::BOTTOM, ModelType::PL) << ")" << endl;
    cout << "   - ChaDen (l): " << getChannelDensity(Orient::LEFT, ModelType::PL) 
         << " (" << getTrackDemand(Orient::LEFT, ModelType::PL) << " " << getTrackSupply(Orient::LEFT, ModelType::PL) << ")" << endl;
    cout << "   - RUDY      : " << getRUDY() << endl;
    cout << "   - #Insts    : " << getNumInstances() << endl;
    cout << "   - #Terms    : " << getNumTerminals() << endl;
    cout << "   - #GNets    : " << getNumGlobalNets() << endl;
    cout << "   - #LNets    : " << getNumLocalNets() << endl;
    cout << "   - LNetDen   : " << getLNetDensity() << endl;
    cout << "   - GNetDen   : " << getGNetDensity() << endl;
    cout << "   - #DRVs     : " << getNumMarkers() << endl;


    
    cout << endl;

}


}
