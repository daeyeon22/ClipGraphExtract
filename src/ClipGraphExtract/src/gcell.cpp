#include "grid.h"
#include "opendb/dbTransform.h"
#include "sta/Sta.hh"
#include "sta/Network.hh"
#include "sta/DelayFloat.hh"
#include "sta/Graph.hh"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"

#include <iostream>


namespace ClipGraphExtract {

using namespace std;
using namespace odb;

Gcell::Gcell(int col, int row) :
    col_(col), row_(row),
    numInsts_(0), numTerms_(0),
    //numLNets_(0), numGNets_(0), 
    wns_(0), tns_(0),
    numLayers_(1), //graph_(nullptr),
    totalCellArea_(0), totalPinArea_(0),
    cellUtil_(0), pinUtil_(0), RUDY_(0),
    lNetRUDY_(0), gNetRUDY_(0), sNetRUDY_(0) {}


void Gcell::setBoundary(Rect rect) { 
    bbox_ = rect;
}

void Gcell::setTotalTrackSupply(int tSup) {
    //rmEGR_.setTrackSupply((int)tSup);
    rmPL_.setTrackSupply((int)tSup);
}

void Gcell::setTrackSupply(int tSup, int layerNum) {
    //rmEGR_.setTrackSupply((int)tSup);
    rmDR_[layerNum].setTrackSupply((int)tSup);
}

void Gcell::setTotalWireCapacity(int wCap) {
    //rmEGR_.setWireCapacity((int)wCap);
    rmPL_.setWireCapacity((int)wCap);
}

void Gcell::setWireCapacity(int wCap, int layerNum) {
    //rmEGR_.setWireCapacity((int)wCap);
    rmDR_[layerNum].setWireCapacity((int)wCap);
}

void Gcell::setNumLayers(int nLyr) {
    numLayers_ = nLyr;
}

set<dbInst*> Gcell::getInstSet() {
    set<dbInst*> instSet;
    copy(insts_.begin(), insts_.end(), inserter(instSet, instSet.end()));
    return instSet;
}

//void Gcell::setGraph(Graph* graph) {
//    graph_ = graph;
//}

bgBox Gcell::getQueryBox() {
    return bgBox(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
}

uint Gcell::getArea() {
    return bbox_.area();    
}

double Gcell::getCellUtil() {
    return cellUtil_;
}

double Gcell::getPinUtil() {
    return pinUtil_;
}

double Gcell::getRUDY() {
    return RUDY_;
}

double Gcell::getLNetRUDY() {
    return lNetRUDY_;
}

double Gcell::getGNetRUDY() {
    return gNetRUDY_;
}

double Gcell::getSNetRUDY() {
    return sNetRUDY_;
}

int Gcell::getCol() {
    return col_;
}

int Gcell::getRow() {
    return row_;
}

int Gcell::getTrackDemand(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: {
            int sumTrack = 0;
            for(int layer = 1; layer <= numLayers_; layer++)
                sumTrack += rmDR_[layer].getTrackDemand(orient);
            return sumTrack;
        }
        //case ModelType::EGR:
        //    return rmEGR_.getTrackDemand(orient);
        case ModelType::TREE:
            return rmPL_.getTrackDemand(orient);
        default:
            return 0; 
    }
}

int Gcell::getTrackDemand(Orient orient, int layerNum, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: 
            return rmDR_[layerNum].getTrackDemand(orient);
        //case ModelType::EGR:
        //    return rmEGR_.getTrackDemand(orient);
        case ModelType::TREE:
            return rmPL_.getTrackDemand(orient);
        default:
            return 0; 
    }
}

int Gcell::getTrackSupply(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: {
            int sumTrack = 0;
            for(int layer = 1; layer <= numLayers_; layer++)
                sumTrack += rmDR_[layer].getTrackSupply(orient);
            return sumTrack;
        }
        //case ModelType::EGR:
        //    return rmEGR_.getTrackSupply(orient);
        case ModelType::TREE:
            return rmPL_.getTrackSupply(orient);
        default:
            return 0; 
    }
}

int Gcell::getTrackSupply(Orient orient, int layerNum, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: 
            return rmDR_[layerNum].getTrackSupply(orient);
        //case ModelType::EGR:
        //    return rmEGR_.getTrackSupply(orient);
        case ModelType::TREE:
            return rmPL_.getTrackSupply(orient);
        default:
            return 0; 
    }
}

int Gcell::getWireCapacity(ModelType type) {
    switch(type) {
        case ModelType::ROUTE: {
            double sumWireCap = 0;
            for(int layer = 1; layer <= numLayers_; layer++)
                sumWireCap += rmDR_[layer].getWireCapacity();
            return sumWireCap;
        }
        //case ModelType::EGR:
        //    return rmEGR_.getWireCapacity();
        case ModelType::TREE:
            return rmPL_.getWireCapacity();
        default:
            return 0; 
    }
}

int Gcell::getWireCapacity(int layerNum, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: 
            return rmDR_[layerNum].getWireCapacity();
        //case ModelType::EGR:
        //    return rmEGR_.getWireCapacity();
        case ModelType::TREE:
            return rmPL_.getWireCapacity();
        default:
            return 0; 
    }
}

int Gcell::getNumInsts() {
    return numInsts_;
}

double Gcell::getViaUtil() {
    double totalViaUtil = 0;
    for(int layer = 1; layer < numLayers_; layer++)
        totalViaUtil += rViaUtil_[layer] + sViaUtil_[layer] + pViaUtil_[layer];
    
    return totalViaUtil/(numLayers_-1);
}

double Gcell::getViaUtil(int layerNum) {
    double totalViaUtil = rViaUtil_[layerNum] + sViaUtil_[layerNum] + pViaUtil_[layerNum];
    if(totalViaUtil > 1) totalViaUtil = 1;
    
    return totalViaUtil;
}

int Gcell::getNumTerms() {
    return numTerms_;
}

int Gcell::getNumNets() {
    return (int)(rsmts_.size());
}

int Gcell::getNumLNets() {
    int count =0;
    for(RSMT* rsmt : rsmts_) {
        if(rsmt->isLocalNet())
            count++;
    }
    return count;
}

int Gcell::getNumGNets() {
    int count = getNumNets() - getNumLNets();
    return count;
}

int Gcell::getNumMarkers() {
    return (int)(markers_.size());
}

void Gcell::getNumMarkers(int &lnet, int &gnet, int &inst) {
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
                } else if (mark->getFromTag() == Marker::Tag::RWoN || mark->getFromTag() == Marker::Tag::SWoN) {
                    if(mark->getFromNet() != NULL) {
                        if(mark->getFromNet()->isLocalNet())
                            lnet++;
                        else
                            gnet++;
                    } else {
                        lnet++;   
                    }
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

double Gcell::getWireUtil(ModelType type) {
     switch(type) {
        case ModelType::ROUTE:  {
            double sumWireUtil = 0;
            for(int layer = 1; layer <= numLayers_; layer++)
                sumWireUtil += rmDR_[layer].getWireUtil();
    
            return sumWireUtil/numLayers_;
        }
        //case ModelType::EGR:
        //    return rmEGR_.getWireUtil();
        case ModelType::TREE:
            return rmPL_.getWireUtil();
        default:
            return 0; 
    }   
}


double Gcell::getWireUtil(int layerNum, ModelType type) {
     switch(type) {
        case ModelType::ROUTE:{
            if(layerNum <= numLayers_)
                return rmDR_[layerNum].getWireUtil(); 
            else 
                return 1;
        }
        //case ModelType::EGR:
        //    return rmEGR_.getWireUtil();
        case ModelType::TREE:
            return rmPL_.getWireUtil();
        default:
            return 0; 
    }   
}

double Gcell::getLNetUtil(ModelType type) {
    int wireCap = 0;
    int wireLen = 0;
    
    if(type == ModelType::TREE) {
        wireCap = rmPL_.getWireCapacity();
        for(RSMT* rsmt : rsmts_) {
            if(rsmt->isLocalNet()) {
                wireLen += rsmt->getWireLengthRSMT();
            }
        }
    } else if(type == ModelType::ROUTE) {
        for(int layer = 1; layer <= numLayers_; layer++)
            wireCap += rmDR_[layer].getWireCapacity();
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

double Gcell::getGNetUtil(ModelType type) {
    //int wireCap = rmPL_.getWireCapacity();
    //int wireLen = rmPL_.getWireLength();
    return getWireUtil(type) - getLNetUtil(type);
}

double Gcell::getChanUtil(ModelType type) {
    double chanDenL = getChanUtil(Orient::LEFT, type);
    double chanDenU = getChanUtil(Orient::TOP, type);
    double chanDenR = getChanUtil(Orient::RIGHT, type);
    double chanDenB = getChanUtil(Orient::BOTTOM, type);
    double chanDen = (chanDenL + chanDenU + chanDenR + chanDenB) / 4;
    return chanDen;
}

double Gcell::getChanUtil(int layerNum, ModelType type) {
    if(layerNum <= numLayers_) {
        double chanDenL = getChanUtil(Orient::LEFT, layerNum, type);
        double chanDenU = getChanUtil(Orient::TOP, layerNum, type);
        double chanDenR = getChanUtil(Orient::RIGHT, layerNum, type);
        double chanDenB = getChanUtil(Orient::BOTTOM, layerNum, type);
        double chanDen = (chanDenL + chanDenU + chanDenR + chanDenB) / 4;
        return chanDen;
    } else 
        return 1;
}

double Gcell::getChanUtil(Orient orient, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: {
            double sumChanUtil = 0;
            for(int layer = 1; layer <= numLayers_; layer++)
                sumChanUtil += rmDR_[layer].getChanUtil(orient);
            return sumChanUtil/numLayers_;
        }
        //case ModelType::EGR:
        //    return rmEGR_.getChanUtil(orient);
        case ModelType::TREE:
            return rmPL_.getChanUtil(orient);
        default:
            return 0.0; 
    }
}

double Gcell::getChanUtil(Orient orient, int layerNum, ModelType type) {
    switch(type) {
        case ModelType::ROUTE: 
            return rmDR_[layerNum].getChanUtil(orient);
        //case ModelType::EGR:
        //    return rmEGR_.getChanUtil(orient);
        case ModelType::TREE:
            return rmPL_.getChanUtil(orient);
        default:
            return 0.0; 
    }
}

double Gcell::getChanUtilV(ModelType type) {
    double chanDenU = getChanUtil(Orient::TOP, type);
    double chanDenB = getChanUtil(Orient::BOTTOM, type);
    double chanDen = (chanDenU + chanDenB) / 2;
    return chanDen;
}


double Gcell::getChanUtilH(ModelType type) {
    double chanDenL = getChanUtil(Orient::LEFT, type);
    double chanDenR = getChanUtil(Orient::RIGHT, type);
    double chanDen = (chanDenL + chanDenR) / 2;
    return chanDen;
}

double Gcell::getAvgTerms() {
    if(getNumNets() > 0) 
        return 1.0* getNumTerms() / getNumNets();
    else
        return 0.0;
}

double Gcell::getBufferUtil() {
    // TODO
    return 0.0;
}


double Gcell::getTNS() {
    return tns_;    
}

double Gcell::getWNS() {
    return wns_;
}

//sta::dbSta* db_sta) {
//
//    // TODO
//    sta::Graph* sta_graph = db_sta->ensureGraph();
//    double tns = 0;
//    for(dbInst* tarInst : insts_) {
//        for(dbITerm* tarITerm : tarInst->getITerms()) {
//            sta::Vertex* tarVertex = sta_graph->vertex(tarITerm->staVertexId());
//            if(tarVertex == NULL)
//                continue;
//            sta::Slack slk = db_sta->vertexSlack(tarVertex, sta::MinMax::max());
//            double slkValue = (double)sta::delayAsFloat(slk);
//
//            if(slkValue < 0) {
//                tns += slkValue;
//            }
//        }
//    }
//    return tns;
//}
//
//double Gcell::getWNS(sta::dbSta* db_sta) {
//    sta::Graph* sta_graph = db_sta->ensureGraph();
//    double wns = 0;
//    for(dbInst* tarInst : insts_) {
//        for(dbITerm* tarITerm : tarInst->getITerms()) {
//            sta::Vertex* tarVertex = sta_graph->vertex(tarITerm->staVertexId());
//            if(tarVertex == NULL)
//                continue;
//            sta::Slack slk = db_sta->vertexSlack(tarVertex, sta::MinMax::max());
//            double slkValue = (double)sta::delayAsFloat(slk);
//
//            if(slkValue < 0) {
//                wns = min(wns, slkValue);
//            }
//        }
//    }
//    return wns;
//}



double Gcell::getClkRatio() {
    // TODO
    int numInsts = insts_.size();
    if(numInsts>0) {
        int numClked = 0;
        for(dbInst* tarInst : insts_) {
            dbMaster* tarMaster = tarInst->getMaster();
            if(tarMaster->isSequential()) {
                numClked++;
            }
        }
        return 1.0 * numClked / numInsts;
    } else {
        return 0.0;
    }
}

void Gcell::updateTimingInfo(unordered_map<dbInst*, double> &slack) {


    //cout << "# of insts = " << insts_.size() << endl;
    for(dbInst* tarInst : insts_) {
        if(slack.find(tarInst) == slack.end()) {
            cout << "?" << endl;
            continue;
        }
        wns_ = max(wns_, slack[tarInst]);
        tns_ += slack[tarInst];
    }
}



//void Gcell::saveGraph(string dirPath, string fileName) {
//    string vertFeaFile = dirPath + "/" + fileName + ".x";
//    string edgeIdxFile = dirPath + "/" + fileName + ".edge_index";
//    string edgeAttFile = dirPath + "/" + fileName + ".edge_attr";
//
//    graph_->saveNodeFeaFile(vertFeaFile);
//    graph_->saveEdgeIdxFile(edgeIdxFile);
//    graph_->saveEdgeAttFile(edgeAttFile);
//}

void
Gcell::extractPlaceFeature(BoxRtree<odb::dbInst*> *rtree) {
    vector<pair<bgBox, dbInst*>> queryResults;

    // Query
    rtree->query(bgi::intersects(getQueryBox()), back_inserter(queryResults));

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
            numTerms_++;
            totalPinArea_ += pinArea;
        }

        //
        //if(totalPinArea_ > getArea()) {
        //    cout << "HERE" << endl;
        //}
        numInsts_++;
        totalCellArea_ += cellArea; 
    }
    //
    cellUtil_ = 1.0* totalCellArea_ / getArea();
    pinUtil_ = 1.0 * totalPinArea_ / getArea();

    //cout << "extractPL done" << endl;
}


void
Gcell::extractRouteFeature(unordered_map<int, SegRtree<odb::dbNet*>> *rRtree, 
                           int maxTechLayer, int maxRouteLayer) {
    bgSeg lb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMin(), bbox_.yMax()));
    bgSeg rb(bgPoint(bbox_.xMax(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg tb(bgPoint(bbox_.xMin(), bbox_.yMax()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg bb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMin()));



    for(int layerNum = 1; layerNum < maxTechLayer; layerNum++) {
        vector<pair<bgSeg, dbNet*>> rQueryResults;
        
        if(layerNum <= maxRouteLayer) {
            if(rRtree->find(layerNum) != rRtree->end()) {
                (*rRtree)[layerNum].query(bgi::intersects(getQueryBox()), back_inserter(rQueryResults));

                for(auto &val : rQueryResults) {
                    bgSeg wire_seg = val.first;
                    dbNet* net = val.second;

                    try {
                        if(bg::intersects(lb, wire_seg)) {
                            rmDR_[layerNum].incrTrackDemand(Orient::LEFT);
                        }
                        if(bg::intersects(rb, wire_seg)) {
                            rmDR_[layerNum].incrTrackDemand(Orient::RIGHT);
                        }
                        if(bg::intersects(tb, wire_seg)) {
                            rmDR_[layerNum].incrTrackDemand(Orient::TOP);
                        }
                        if(bg::intersects(bb, wire_seg)) {
                            rmDR_[layerNum].incrTrackDemand(Orient::BOTTOM);
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
                        int wl = (x1-x0) + (y1-y0);
                        rmDR_[layerNum].addWireLength(wl);

                    } catch(boost::numeric::negative_overflow &e) {
                        cout << e.what() << endl;
                        cout << bg::wkt(rb) << " " << bg::wkt(wire_seg) << endl;
                    } catch(boost::numeric::positive_overflow &e) {
                        cout << e.what() << endl;
                    }
                }
                
                /*
                cout << "layerNum: m" << layerNum << " ";
                cout << rmDR_[layerNum].getTrackSupply(Orient::LEFT) << " "
                     << rmDR_[layerNum].getTrackDemand(Orient::LEFT) << " "
                     << rmDR_[layerNum].getTrackSupply(Orient::RIGHT) << " "
                     << rmDR_[layerNum].getTrackDemand(Orient::RIGHT) << " "
                     << rmDR_[layerNum].getTrackSupply(Orient::TOP) << " "
                     << rmDR_[layerNum].getTrackDemand(Orient::TOP) << " "
                     << rmDR_[layerNum].getTrackSupply(Orient::BOTTOM) << " "
                     << rmDR_[layerNum].getTrackDemand(Orient::BOTTOM) << endl;
                */

            }
        }
    }
}

void
Gcell::extractViaFeature(std::unordered_map<int, BoxRtree<odb::dbTechVia*>> *rViaRtreeMap, 
                         std::unordered_map<int, BoxRtree<odb::dbTechVia*>> *sViaRtreeMap,
                         std::unordered_map<int, BoxRtree<odb::dbTechVia*>> *pViaRtreeMap,
                         int maxTechLayer, int maxRouteLayer) {
    
   
    bgPoint minCorner = getQueryBox().min_corner();
    bgPoint maxCorner = getQueryBox().max_corner();

    // cout << minCorner.get<0>()/2000.0 << " ";
    // cout << minCorner.get<1>()/2000.0 << " ";
    // cout << maxCorner.get<0>()/2000.0 << " ";
    // cout << maxCorner.get<1>()/2000.0 << endl;

    for(int layerNum = 1; layerNum < maxTechLayer; layerNum++) {
        vector<pair<bgBox, dbTechVia*>> rQueryResults;
        vector<pair<bgBox, dbTechVia*>> sQueryResults;
        vector<pair<bgBox, dbTechVia*>> pQueryResults;
   
        // cout << "layerNum : " << layerNum << endl;

        // cout << "Signal Vias" << endl;
        // Layers that cannot be used are blocked.
        if(maxRouteLayer <= layerNum) {
            rViaUtil_[layerNum] = 1;
            sViaUtil_[layerNum] = 1;
            pViaUtil_[layerNum] = 1;
        } else {
            if(rViaRtreeMap->find(layerNum) != rViaRtreeMap->end()) {
                (*rViaRtreeMap)[layerNum].query(bgi::intersects(getQueryBox()), back_inserter(rQueryResults));
                for(auto &val : rQueryResults) {
                    dbTechVia* via = val.second;
                    int bottomPitch = via->getBottomLayer()->getPitch();
                    int topPitch = via->getTopLayer()->getPitch();

                    // The bbox is "square".
                    int numViaGrid = (bbox_.dx()/bottomPitch) * (bbox_.dy()/topPitch);

                    rViaUtil_[layerNum] += 1.0/numViaGrid;
                    
                    // bgBox box = val.first;
            
                    // bgPoint minCorner = box.min_corner();
                    // bgPoint maxCorner = box.max_corner();
            
                    // cout << (minCorner.get<0>()+maxCorner.get<0>())/4000.0 << " ";
                    // cout << (minCorner.get<1>()+maxCorner.get<1>())/4000.0 << " ";
            
                    // cout << via->getBottomLayer()->getRoutingLevel() << endl;
                }
            } else rViaUtil_[layerNum] = 0;
            
            // cout << endl;
            // cout << "Special Vias" << endl;
            if(sViaRtreeMap->find(layerNum) != sViaRtreeMap->end()) {
                (*sViaRtreeMap)[layerNum].query(bgi::intersects(getQueryBox()), back_inserter(sQueryResults));
                for(auto &val : sQueryResults) {
                    dbTechVia* via = val.second;
                    int bottomPitch = via->getBottomLayer()->getPitch();
                    int topPitch = via->getTopLayer()->getPitch();
                    
                    int numViaGrid = (bbox_.dx()/bottomPitch) * (bbox_.dy()/topPitch);

                    sViaUtil_[layerNum] += 1.0/numViaGrid;
                }
            } else sViaUtil_[layerNum] = 0;
            // cout << endl;
            // cout << "Power Vias" << endl;
            if(pViaRtreeMap->find(layerNum) != pViaRtreeMap->end()) {
                (*pViaRtreeMap)[layerNum].query(bgi::intersects(getQueryBox()), back_inserter(pQueryResults));
                for(auto &val : pQueryResults) {
                    dbTechVia* via = val.second;
                    int bottomPitch = via->getBottomLayer()->getPitch();
                    int topPitch = via->getTopLayer()->getPitch();
                    
                    // The bbox is "square".
                    int numViaGrid = (bbox_.dx()/bottomPitch) * (bbox_.dy()/topPitch);

                    pViaUtil_[layerNum] += 1.0/numViaGrid;
                }
            } else pViaUtil_[layerNum] = 0;
        }
        
//        cout << "Routing Via : " << rViaUtil_[layerNum] << endl;
//        cout << "Special Via : " << sViaUtil_[layerNum] << endl;
//        cout << "Power Via : " << pViaUtil_[layerNum] << endl;
//        cout << endl;
    }
}



void
Gcell::extractPlaceFeature(SegRtree<RSMT*> *rtree) {
    //cout << "extract feature rsmt start" << endl;
    vector<pair<bgSeg, RSMT*>> queryResults;
    clock_t start, end;
    double runtime;
    int c =0;
    runtime = 0;

    start = clock();
    // Query
    rtree->query(bgi::intersects(getQueryBox()), back_inserter(queryResults));
    bgSeg lb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMin(), bbox_.yMax()));
    bgSeg rb(bgPoint(bbox_.xMax(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg tb(bgPoint(bbox_.xMin(), bbox_.yMax()), bgPoint(bbox_.xMax(), bbox_.yMax()));
    bgSeg bb(bgPoint(bbox_.xMin(), bbox_.yMin()), bgPoint(bbox_.xMax(), bbox_.yMin()));

	//cout << bbox_.xMin() << " " << bbox_.yMin() << " " << bbox_.xMax() << " " << bbox_.yMax() << endl;

  
	set<RSMT*> RSMTs;
    int wirelength=0;
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
            int wl = (x1-x0) + (y1-y0);
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
        int area1 = r2.intersect(r1).area();
        int area2 = r1.area();
        double dn = myRSMT->getWireUniformUtil();
        double R = 1.0* area1 / area2;
        double partialRUDY = dn * R;   
        RUDY_ += partialRUDY;
        //
        rsmts_.push_back(myRSMT);

        if(myRSMT->getNet()->isSpecial()) {
            sNetRUDY_ += partialRUDY;
        } else {
            if(myRSMT->isLocalNet()) {
                lNetRUDY_ += partialRUDY;
            } else {
                gNetRUDY_ += partialRUDY;
            }
        }
    }
    //end = clock();
    //runtime = double(end-start); 
    //cout << "(" << c++ << ") " << runtime << endl;
    //start = end;


    RUDY_ /=  numLayers_;
    lNetRUDY_ /= numLayers_;
    gNetRUDY_ /= numLayers_;
    sNetRUDY_ /= numLayers_;
    //cout << "extract feature rsmt done" << endl;
}







void Gcell::print() {

    int lnet, gnet, inst;
    getNumMarkers(lnet, gnet, inst);

    cout << "GCELL (" << bbox_.xMin() << " " << bbox_.yMin() << ") (" << bbox_.xMax() << " " << bbox_.yMax() << ")" << endl;
    cout << "   - CellDen   : " << getCellUtil() << endl;
    cout << "   - PinDen    : " << getPinUtil() << endl;
    cout << "   - #Insts    : " << getNumInsts() << endl;
    cout << "   - #Terms    : " << getNumTerms() << endl;
    cout << "   - #GNets    : " << getNumGNets() << endl;
    cout << "   - #LNets    : " << getNumLNets() << endl;
    cout << "   - RUDY      : " << getRUDY() << endl;
    cout << "   Measured using RSMT " << endl;
    cout << "   - LNetDen   : " << getLNetUtil(ModelType::TREE) << endl;
    cout << "   - GNetDen   : " << getGNetUtil(ModelType::TREE) << endl;
    cout << "   - WireDen   : " << getWireUtil(ModelType::TREE) 
         << "(" << getWireCapacity(ModelType::TREE) << ")" <<  endl;
    cout << "   - ChaDen (u): " << getChanUtil(Orient::TOP, ModelType::TREE) 
         << " (" << getTrackDemand(Orient::TOP, ModelType::TREE) << " " << getTrackSupply(Orient::TOP, ModelType::TREE) << ")" << endl;
    cout << "   - ChaDen (r): " << getChanUtil(Orient::RIGHT, ModelType::TREE)
         << " (" << getTrackDemand(Orient::RIGHT, ModelType::TREE) << " " << getTrackSupply(Orient::RIGHT, ModelType::TREE) << ")" << endl;
    cout << "   - ChaDen (b): " << getChanUtil(Orient::BOTTOM, ModelType::TREE)
         << " (" << getTrackDemand(Orient::BOTTOM, ModelType::TREE) << " " << getTrackSupply(Orient::BOTTOM, ModelType::TREE) << ")" << endl;
    cout << "   - ChaDen (l): " << getChanUtil(Orient::LEFT, ModelType::TREE) 
         << " (" << getTrackDemand(Orient::LEFT, ModelType::TREE) << " " << getTrackSupply(Orient::LEFT, ModelType::TREE) << ")" << endl;
     cout << "   Measured using DR " << endl;
    cout << "   - LNetDen   : " << getLNetUtil(ModelType::ROUTE) << endl;
    cout << "   - GNetDen   : " << getGNetUtil(ModelType::ROUTE) << endl;
    cout << "   - WireDen   : " << getWireUtil(ModelType::ROUTE) 
         << "(" << getWireCapacity(ModelType::ROUTE) << ")" <<  endl;
    cout << "   - ChaDen (u): " << getChanUtil(Orient::TOP, ModelType::ROUTE) 
         << " (" << getTrackDemand(Orient::TOP, ModelType::ROUTE) << " " << getTrackSupply(Orient::TOP, ModelType::ROUTE) << ")" << endl;
    cout << "   - ChaDen (r): " << getChanUtil(Orient::RIGHT, ModelType::ROUTE)
         << " (" << getTrackDemand(Orient::RIGHT, ModelType::ROUTE) << " " << getTrackSupply(Orient::RIGHT, ModelType::ROUTE) << ")" << endl;
    cout << "   - ChaDen (b): " << getChanUtil(Orient::BOTTOM, ModelType::ROUTE)
         << " (" << getTrackDemand(Orient::BOTTOM, ModelType::ROUTE) << " " << getTrackSupply(Orient::BOTTOM, ModelType::ROUTE) << ")" << endl;
    cout << "   - ChaDen (l): " << getChanUtil(Orient::LEFT, ModelType::ROUTE) 
         << " (" << getTrackDemand(Orient::LEFT, ModelType::ROUTE) << " " << getTrackSupply(Orient::LEFT, ModelType::ROUTE) << ")" << endl;
 
    cout << "   - #DRVs     : " << getNumMarkers() << endl;
    cout << "   -   due to local net = " << lnet << endl;
    cout << "   -   due to global net = " << gnet << endl;
    cout << "   -   due to instance = " << inst << endl;


    
    cout << endl;

}


}
