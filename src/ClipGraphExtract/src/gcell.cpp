#include "grid.h"
#include "opendb/dbTransform.h"
#include <iostream>


namespace ClipGraphExtract {

using namespace std;
using namespace odb;

Gcell::Gcell() :
    numInstances_(0), numTerminals_(0), 
    //numLocalNets_(0), numGlobalNets_(0), 
    numLayers_(1), graph_(nullptr),
    totalCellArea_(0), totalPinArea_(0),
    cellDensity_(0), pinDensity_(0), RUDY_(0),
    lnetRUDY_(0), gnetRUDY_(0), snetRUDY_(0) {}


void Gcell::setBoundary(Rect rect) { 
    bbox_ = rect;
}

void Gcell::setTrackSupply(int tSup) {
    //rmEGR_.setTrackSupply((uint)tSup);
    rmPL_.setTrackSupply((uint)tSup);
    rmDR_.setTrackSupply((uint)tSup);
}

void Gcell::setWireCapacity(int wCap) {
    //rmEGR_.setWireCapacity((uint)wCap);
    rmPL_.setWireCapacity((uint)wCap);
    rmDR_.setWireCapacity((uint)wCap);
}

void Gcell::setNumLayers(int nLyr) {
    numLayers_ = nLyr;
}


set<dbInst*> Gcell::getInstSet() {
    return set<dbInst*>(insts_.begin(), insts_.end());
}

void Gcell::setGraph(Graph* graph) {
    graph_ = graph;
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

double Gcell::getLNetRUDY() {
    return lnetRUDY_;
}

double Gcell::getGNetRUDY() {
    return gnetRUDY_;
}

double Gcell::getSNetRUDY() {
    return snetRUDY_;
}

uint Gcell::getTrackDemand(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::DR: 
            return rmDR_.getTrackDemand(orient);
        //case ModelType::EGR:
        //    return rmEGR_.getTrackDemand(orient);
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
        //case ModelType::EGR:
        //    return rmEGR_.getTrackSupply(orient);
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
        //case ModelType::EGR:
        //    return rmEGR_.getWireCapacity();
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

void Gcell::getNumMarkers(uint &lnet, uint &gnet, uint &inst) {
    lnet =0;
    gnet =0;
    inst =0;
    for(Marker* mark : markers_) {
        switch(mark->getCategory()) {
            case Marker::Category::L2L:
                lnet++; break;
            case Marker::Category::L2I:
                lnet++; inst++; break;
            case Marker::Category::L2G:
                lnet++; gnet++; break;
            case Marker::Category::G2I:
                gnet++; inst++; break;
            case Marker::Category::I2I:
                inst++; break;
            case Marker::Category::SELF: {
                if(mark->getFromTag() == Marker::Tag::BoC || mark->getFromTag() == Marker::Tag::PoC) {
                    inst++;
                } else if (mark->getFromTag() == Marker::Tag::RWoN) {
                    if(mark->getFromNet()->isLocalNet())
                        lnet++;
                    else
                        gnet++;
                }
                break;
            }
            case Marker::Category::ERR:
                break;
            default:
                break;
        }
    }
}



double Gcell::getWireDensity(ModelType type) {
     switch(type) {
        case ModelType::DR: 
            return rmDR_.getWireDensity();
        //case ModelType::EGR:
        //    return rmEGR_.getWireDensity();
        case ModelType::PL:
            return rmPL_.getWireDensity();
        default:
            return 0; 
    }   
}

double Gcell::getLNetDensity(ModelType type) {
    uint wireCap = 0;
    uint wireLen = 0;
    
    if(type == ModelType::PL) {
        wireCap = rmPL_.getWireCapacity();
        for(RSMT* rsmt : rsmts_) {
            if(rsmt->isLocalNet()) {
                wireLen += rsmt->getWireLengthRSMT();
            }
        }
    } else if(type == ModelType::DR) {
        wireCap = rmDR_.getWireCapacity();
        for(RSMT* rsmt : rsmts_) {
            if(rsmt->isLocalNet()) {
                dbWire* wire = rsmt->getNet()->getWire();
                if(wire != NULL) {
                    wireLen += wire->getLength();
                }
            }
        }
    } else {
        cout << "???" << endl;
        exit(0);
    }

    double c = min(1.0, 3.0 / numLayers_);
    return 1.0 * wireLen / (c*wireCap);

}

double Gcell::getGNetDensity(ModelType type) {
    //uint wireCap = rmPL_.getWireCapacity();
    //uint wireLen = rmPL_.getWireLength();
    return getWireDensity(type) - getLNetDensity(type);
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
        //case ModelType::EGR:
        //    return rmEGR_.getChannelDensity(orient);
        case ModelType::PL:
            return rmPL_.getChannelDensity(orient);
        default:
            return 0.0; 
    }
}

double Gcell::getChannelDensityV(ModelType type) {
    double chanDenU = getChannelDensity(Orient::TOP, type);
    double chanDenB = getChannelDensity(Orient::BOTTOM, type);
    double chanDen = (chanDenU + chanDenB) / 2;
    return chanDen;
}


double Gcell::getChannelDensityH(ModelType type) {
    double chanDenL = getChannelDensity(Orient::LEFT, type);
    double chanDenR = getChannelDensity(Orient::RIGHT, type);
    double chanDen = (chanDenL + chanDenR) / 2;
    return chanDen;
}





void
Gcell::extractFeaturePL(BoxRtree<odb::dbInst*> &rtree) {
    vector<pair<bgBox, dbInst*>> queryResults;

    //bbox_.print();
   
    //bgBox qBox = getQueryBox();
    //cout << bg::get<0,0>(qBox) << " " << bg::get<0,1>(qBox) << " " << bg::get<1,0>(qBox) << " " << bg::get<1,1>(qBox) << endl;


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

    //cout << "extractPL done" << endl;
}


void
Gcell::extractFeatureDR(SegRtree<odb::dbNet*> &rtree) {

    vector<pair<bgSeg, dbNet*>> queryResults;

    // Query
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));

    bgSeg lb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMin(), bbox_.yMax()));
    bgSeg rb(bgPoint(bbox_.xMax(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg tb(bgPoint(bbox_.xMin(), bbox_.yMax()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg bb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMin()));


    //cout << "# of intersections : " << queryResults.size() << endl;

    for(auto &val : queryResults) {
        bgSeg wire_seg = val.first;
        dbNet* net = val.second;

        try {
            if(bg::intersects(lb, wire_seg)) {
                rmDR_.incrTrackDemand(Orient::LEFT);
            }
            if(bg::intersects(rb, wire_seg)) {
                rmDR_.incrTrackDemand(Orient::RIGHT);
            }
            if(bg::intersects(tb, wire_seg)) {
                rmDR_.incrTrackDemand(Orient::TOP);
            }
            if(bg::intersects(bb, wire_seg)) {
                rmDR_.incrTrackDemand(Orient::BOTTOM);
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
            rmDR_.addWireLength(wl);

        } catch(boost::numeric::negative_overflow &e) {
            cout << e.what() << endl;
            cout << bg::wkt(rb) << " " << bg::wkt(wire_seg) << endl;
        } catch(boost::numeric::positive_overflow &e) {
            cout << e.what() << endl;
        }
    
    }
}


void
Gcell::extractFeatureRSMT(SegRtree<RSMT*> &rtree) {
    //cout << "extract feature rsmt start" << endl;
    vector<pair<bgSeg, RSMT*>> queryResults;

    // Query
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));
    bgSeg lb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMin(), bbox_.yMax()));
    bgSeg rb(bgPoint(bbox_.xMax(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg tb(bgPoint(bbox_.xMin(), bbox_.yMax()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg bb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMin()));

	//cout << bbox_.xMin() << " " << bbox_.yMin() << " " << bbox_.xMax() << " " << bbox_.yMax() << endl;
   
	set<RSMT*> RSMTs;
    uint wirelength=0;
    for(auto &val : queryResults) {
        bgSeg wire_seg = val.first;
        RSMT* myRSMT = val.second;
        try {
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
            
            x0 = max(x0, bbox_.xMin());
            y0 = max(y0, bbox_.yMin());
            x1 = min(x1, bbox_.xMax());
            y1 = min(y1, bbox_.yMax());
            // intersection wirelength
            uint wl = (x1-x0) + (y1-y0);
            rmPL_.addWireLength(wl);
            RSMTs.insert(myRSMT);

        } catch(boost::numeric::negative_overflow &e) {
            cout << e.what() << endl;
            cout << bg::wkt(rb) << " " << bg::wkt(wire_seg) << endl;
        } catch(boost::numeric::positive_overflow &e) {
            cout << e.what() << endl;
        }


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

        if(myRSMT->getNet()->isSpecial()) {
            snetRUDY_ += partialRUDY;
        } else {
            if(myRSMT->isLocalNet()) {
                lnetRUDY_ += partialRUDY;
            } else {
                gnetRUDY_ += partialRUDY;
            }
        }


        //if(myRSMT->isLocalNet())
        //    numLocalNets_++;
        //else
        //    numGlobalNets_++;

    }

    RUDY_ = RUDY_ / numLayers_;
    //cout << "extract feature rsmt done" << endl;
}







void Gcell::print() {

    uint lnet, gnet, inst;
    getNumMarkers(lnet, gnet, inst);

    cout << "GCELL (" << bbox_.xMin() << " " << bbox_.yMin() << ") (" << bbox_.xMax() << " " << bbox_.yMax() << ")" << endl;
    cout << "   - CellDen   : " << getCellDensity() << endl;
    cout << "   - PinDen    : " << getPinDensity() << endl;
    cout << "   - #Insts    : " << getNumInstances() << endl;
    cout << "   - #Terms    : " << getNumTerminals() << endl;
    cout << "   - #GNets    : " << getNumGlobalNets() << endl;
    cout << "   - #LNets    : " << getNumLocalNets() << endl;
    cout << "   - RUDY      : " << getRUDY() << endl;
    cout << "   Measured using RSMT " << endl;
    cout << "   - LNetDen   : " << getLNetDensity(ModelType::PL) << endl;
    cout << "   - GNetDen   : " << getGNetDensity(ModelType::PL) << endl;
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
     cout << "   Measured using DR " << endl;
    cout << "   - LNetDen   : " << getLNetDensity(ModelType::DR) << endl;
    cout << "   - GNetDen   : " << getGNetDensity(ModelType::DR) << endl;
    cout << "   - WireDen   : " << getWireDensity(ModelType::DR) 
         << "(" << getWireCapacity(ModelType::DR) << ")" <<  endl;
    cout << "   - ChaDen (u): " << getChannelDensity(Orient::TOP, ModelType::DR) 
         << " (" << getTrackDemand(Orient::TOP, ModelType::DR) << " " << getTrackSupply(Orient::TOP, ModelType::DR) << ")" << endl;
    cout << "   - ChaDen (r): " << getChannelDensity(Orient::RIGHT, ModelType::DR)
         << " (" << getTrackDemand(Orient::RIGHT, ModelType::DR) << " " << getTrackSupply(Orient::RIGHT, ModelType::DR) << ")" << endl;
    cout << "   - ChaDen (b): " << getChannelDensity(Orient::BOTTOM, ModelType::DR)
         << " (" << getTrackDemand(Orient::BOTTOM, ModelType::DR) << " " << getTrackSupply(Orient::BOTTOM, ModelType::DR) << ")" << endl;
    cout << "   - ChaDen (l): " << getChannelDensity(Orient::LEFT, ModelType::DR) 
         << " (" << getTrackDemand(Orient::LEFT, ModelType::DR) << " " << getTrackSupply(Orient::LEFT, ModelType::DR) << ")" << endl;
 
    cout << "   - #DRVs     : " << getNumMarkers() << endl;
    cout << "   -   due to local net = " << lnet << endl;
    cout << "   -   due to global net = " << gnet << endl;
    cout << "   -   due to instance = " << inst << endl;


    
    cout << endl;

}


}
