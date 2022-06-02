#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include "opendb/geom.h"
#include "sta/Sta.hh"
#include "sta/Network.hh"
#include "sta/DelayFloat.hh"
#include "sta/Graph.hh"
#include "sta/Clock.hh"


#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "instGraph.h"
#include "grid.h"
#include "opendb/dbTransform.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <regex>
#include <fstream>
#include <vector>
#include <cassert>




namespace ClipGraphExtract {

using namespace odb;
using namespace std;


void initInstRtree(odb::dbDatabase* db, BoxRtree<dbInst*> &instRtree) {
    dbBlock* block = db->getChip()->getBlock();
    // make instRtree
    for( dbInst* inst : block->getInsts() ) {
        dbBox* bBox = inst->getBBox();
        bgBox b (bgPoint(bBox->xMin(), bBox->yMin()),
                    bgPoint(bBox->xMax(), bBox->yMax()));
        instRtree.insert( make_pair(b, inst) );
    }
}



void initWireRtree(odb::dbDatabase* db, SegRtree<dbNet*> &rwireRtree, BoxRtree<dbNet*> &swireRtree, 
        BoxRtree<dbTechVia*> &rviaRtree, BoxRtree<dbTechVia*> &sviaRtree) {

    // make wireRtree
    int x, y, ext;
    dbBlock* block = db->getChip()->getBlock();
    for(dbNet* net : block->getNets()) {
        dbWire* wire = net->getWire();
        if( wire && wire->length() ) {
            int wl = wire->getLength();
            int wl_ = 0;
            dbWireDecoder decoder;
            decoder.begin(wire);
            vector<odb::Point> points;

            dbWireDecoder::OpCode opcode = decoder.peek();
            while(opcode != dbWireDecoder::END_DECODE) {
                bool hasPath=false;
                bool isVia = false;
                switch(opcode) {
                    case dbWireDecoder::PATH:
                        points.clear(); break;
                    case dbWireDecoder::JUNCTION: break;
                    case dbWireDecoder::SHORT: break;

                    case dbWireDecoder::TECH_VIA: {
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        isVia = true;
                        hasPath=true;
                        break;
                    }
                    case dbWireDecoder::VIA: break;
                    case dbWireDecoder::VWIRE: break;
                    case dbWireDecoder::POINT_EXT:
                        decoder.getPoint(x,y,ext);
                        points.push_back(odb::Point(x,y));
                        hasPath=true;
                        break;
                    case dbWireDecoder::BTERM: break;
                    case dbWireDecoder::ITERM: break;
                    case dbWireDecoder::RULE: break;
                    case dbWireDecoder::POINT: {
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        hasPath=true;
                        break;
                    }
                    default: break;
                }
         
                // Wire segment
                if(hasPath && points.size() > 1) {
                    odb::Point pt1 = points[points.size()-2];
                    odb::Point pt2 = points[points.size()-1];
                    
//                    cout << net->getName() << " Type Check : " << net->getSigType().getString() << endl;

                    int xMin = min(pt1.getX(), pt2.getX());
                    int xMax = max(pt1.getX(), pt2.getX());
                    int yMin = min(pt1.getY(), pt2.getY());
                    int yMax = max(pt1.getY(), pt2.getY());
                    int dist = (xMax-xMin) + (yMax-yMin);
                    
                    if(dist != 0) {
                        wl_ += dist;
                        bgSeg wireSeg( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                        bgBox wireBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                        if(net->isSpecial()) swireRtree.insert( make_pair( wireBox, net ) );
                        else rwireRtree.insert( make_pair( wireSeg, net ) );
                    }
                }
                
                // Via point
                if(hasPath && isVia) {
                    dbTechVia* via = decoder.getTechVia();
                    odb::Point pt = points[points.size()-1];

                    int xMin = pt.getX() + via->getBBox()->xMin();
                    int xMax = pt.getX() + via->getBBox()->xMax();
                    int yMin = pt.getY() + via->getBBox()->yMin();
                    int yMax = pt.getY() + via->getBBox()->yMax();
                
//                    cout << xMin << " ";
//                    cout << xMax << " ";
//                    cout << yMin << " ";
//                    cout << yMax << endl;
                    
                    bgBox viaBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                    if(net->isSpecial()) sviaRtree.insert( make_pair( viaBox, via ) );
                    else rviaRtree.insert( make_pair( viaBox, via ) );
                }
                opcode = decoder.next();
            }
        }
    } 
}

void initPWireRtree(odb::dbDatabase* db,
    BoxRtree<dbNet*> &pwireRtree, BoxRtree<dbTechVia*> &pviaRtree) {
    
    dbBlock* block = db->getChip()->getBlock();
    set<dbNet*> findNet;
    dbSet<dbITerm> iterms = block->getITerms();

    for(dbITerm* iterm : iterms){
        string type = iterm->getSigType().getString();
        if (type == "POWER" || type == "GROUND"){
            dbNet* net = iterm->getNet();
            if(findNet.find(net) == findNet.end()) {
                findNet.insert(net);
                dbSet<dbSWire> swires = net->getSWires();
                for(dbSWire* swire : swires) {
                    dbSet<dbSBox> sboxes = swire->getWires();
                    for(dbSBox* sbox : sboxes){
                        int xMin = sbox->xMin();
                        int xMax = sbox->xMax();
                        int yMin = sbox->yMin();
                        int yMax = sbox->yMax();
                        if (sbox->isVia()){
                            dbTechVia* via = sbox->getTechVia();
                            bgBox viaBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                            pviaRtree.insert( make_pair( viaBox, via ) );
                        }
                        else {
                            bgBox wireBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                            pwireRtree.insert( make_pair( wireBox, net ) );
                        }
                    }
                }
            }
        }
    }
}


void initGcellRtree(Grid* grid, BoxRtree<Gcell*> &gcellRtree) {
    // make gcellRtree
    for( Gcell* gcell : grid->getGcells() ) {
        bgBox gcellBox = gcell->getQueryBox();
        gcellRtree.insert( make_pair(gcellBox, gcell) );
    }
}


void initRsmtRtree(Grid* grid, SegRtree<RSMT*> &rsmtRtree) {

    // make rsmtRtree
    for( RSMT* rsmt : grid->getRSMTs() ) {       
        vector<Rect> segments = rsmt->getSegments();
        // insert segments into rtree
        for(Rect& seg : segments) {
            // update (1) #cut-nets (2) wire utilization
            bgSeg bgseg(bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()));
            rsmtRtree.insert( make_pair( bgseg, rsmt ) );
        }
    }

}


void getTimingInfo(dbDatabase* db_, sta::dbSta* sta_, 
        double &clockPeriod_,
        unordered_map<dbInst*, double> &absSlack_, 
        unordered_map<dbInst*, double> &relSlack_, 
        unordered_map<dbInst*, bool> &isCritical_) {

    sta_->updateTiming(false);
    sta::dbNetwork* network = sta_->getDbNetwork();
    sta::Graph* graph = sta_->ensureGraph();

    sta::Vertex* wstVertex;
    sta::Slack wstSlack, totNegSlack;
    sta_->worstSlack(sta::MinMax::max(), wstSlack, wstVertex);
    totNegSlack = sta_->totalNegativeSlack(sta::MinMax::max());
    cout << "Worst Slack    : " << wstSlack << endl;
    cout << "Total NegSlack : " << totNegSlack << endl;

    double wstSlk = abs((double)sta::delayAsFloat(wstSlack));


    string clockName = "clk";
    const sta::Clock* staClock = sta_->findClock(clockName.c_str());
    if(staClock == NULL) {
        cout << "Cannot find clock (" << clockName << ")" << endl;
        clockPeriod_ = 1.0;
    } else {
        cout << staClock->name() << " -period " << staClock->period() << endl;
        clockPeriod_ = staClock->period();
    }

    dbBlock* block = db_->getChip()->getBlock();
    // Get slack info.
    for(dbInst* inst : block->getInsts()) {
        absSlack_[inst] = 0.0;
        for(dbITerm* iterm : inst->getITerms()) {
            if(iterm->getIoType() == dbSigType::POWER || iterm->getIoType() == dbSigType::GROUND)
                continue;
            uint32_t vertexId = iterm->staVertexId();
            sta::Vertex* vertex = graph->vertex(vertexId);
            if(vertex == NULL)
                continue;
            sta::Slack staSlk = sta_->vertexSlack(vertex, sta::MinMax::max());
            double slack = (double)sta::delayAsFloat(staSlk);
         
            absSlack_[inst] = min(absSlack_[inst], slack);
            //cout << slack << " " << absSlack_[inst] << endl;
        }


        absSlack_[inst] = abs(absSlack_[inst]);
        relSlack_[inst] = absSlack_[inst] / clockPeriod_;
        isCritical_[inst] = absSlack_[inst] > 0.0 ?  true : false;
       

        absSlack_[inst] *= 1e+9;
        if(isCritical_[inst]) {
            cout << inst->getName() << " is in a critical path (" << absSlack_[inst] << ") " << wstSlk  << endl;
        }   


    }
}

void ClipGraphExtractor::init() {
    Flute::readLUT();

    grid_ = (void*) new Grid();
    Grid* grid = (Grid*) grid_;
    dbBlock* block = getDb()->getChip()->getBlock();

   // Calculate grid size = rowHeight * numRows_
    dbSite* site = block->getRows().begin()->getSite();
    int gcellWidth = site->getHeight() * numRows_;
    int gcellHeight = site->getHeight() * numRows_;
    // get routing supply / capacity for each gcell
    dbSet<dbTechLayer> techLayers = getDb()->getTech()->getLayers();
    int numLayer=0;
    int trackSupply=0;
    int wireCapacity=0;
    int minWidth =INT_MAX;
    for(dbTechLayer* layer : techLayers) {
        if(layer->getType() == dbTechLayerType::ROUTING) {
            numLayer++;
            minWidth = min(minWidth, (int)layer->getWidth());
            int minPitch = layer->getPitch();
            int minSpacing = layer->getSpacing();
            int capacity=0;
            int supply=0;
            if(layer->getDirection() == dbTechLayerDir::HORIZONTAL) {
                supply = gcellHeight / minPitch;
                capacity = supply * gcellWidth;
            }else if(layer->getDirection() == dbTechLayerDir::VERTICAL) {
                supply = gcellWidth / minPitch;
                capacity = supply * gcellHeight;
            }

            trackSupply += supply;
            wireCapacity += capacity;

            if(numLayer == maxRouteLayer_)
                break;
        }
    }

    Rect blockArea;
    block->getDieArea(blockArea);
    cout << "Die area"
            << " (" << blockArea.xMin() << " " << blockArea.yMin() << ")"
            << " (" << blockArea.xMax() << " " << blockArea.yMax() << ")" << endl;
    cout << "TrackSupply    : " << trackSupply << endl;
    cout << "WireCapacity   : " << wireCapacity << endl;
    // Get core area


    // Initialize Gcell Grid
    grid->setDb(getDb());
    grid->setSta(getSta());
    grid->setBoundary(blockArea);
    grid->setGcellWidth(gcellWidth);
    grid->setGcellHeight(gcellHeight);
    grid->setWireCapacity(wireCapacity);
    grid->setTrackSupply(trackSupply);
    grid->setNumLayers(numLayer);
    grid->setWireMinWidth(minWidth);
    grid->init();
    //  


}


void ClipGraphExtractor::extract() {

    Grid* grid = (Grid*)grid_;

    dbBlock* block = db_->getChip()->getBlock();
    dbSet<dbTechLayer> techLayers = db_->getTech()->getLayers();
    double dbu = block->getDbUnitsPerMicron();

    set<dbNet*> findNet;
    dbSet<dbITerm> iterms = block->getITerms();

    // Init Rtrees
    SegRtree<dbNet*> rWireRtree;
    BoxRtree<dbTechVia*> rViaRtree;

    BoxRtree<dbNet*> sWireRtree;
    BoxRtree<dbTechVia*> sViaRtree;

    BoxRtree<dbNet*> pWireRtree;
    BoxRtree<dbTechVia*> pViaRtree;
    
    BoxRtree<dbInst*> instRtree;
    BoxRtree<Gcell*> gcellRtree;
    SegRtree<RSMT*> rsmtRtree;
    
    initWireRtree(db_, rWireRtree, sWireRtree, rViaRtree, sViaRtree);
    initPWireRtree(db_, pWireRtree, pViaRtree);
    initInstRtree(db_, instRtree);
    initGcellRtree((Grid*)grid_, gcellRtree);
    initRsmtRtree((Grid*)grid_, rsmtRtree);
    // 
    //
    //unordered_map<dbInst*, bool> isCritical_;
    //unordered_map<dbInst*, double> absSlack_;
    double clockPeriod_;
    // Get timing info.
    getTimingInfo(db_, sta_, clockPeriod_, absSlack_, relSlack_, isCritical_);

    for( dbNet* net : block->getNets() ) {
        // search overlapping gcells
        RSMT* rsmt = grid->getRSMT(net);
        rsmt->searchOverlaps(&gcellRtree);
    }


    // add Instance in gcell 
    for( Gcell* gcell : grid->getGcells() ) {
        gcell->extractPlaceFeature(&instRtree);
        gcell->extractPlaceFeature(&rsmtRtree);
        gcell->extractRouteFeature(&rWireRtree);
        gcell->extractViaFeature(&rViaRtree, &sViaRtree, &pViaRtree);
        gcell->updateTimingInfo(absSlack_);
    }


    int gcellWidth = grid->getGcellWidth();
    int gcellHeight = grid->getGcellHeight();


    // Encode dbMaster to int
    int count=0;
    unordered_map<dbMaster*, int> typeEncoder;
    for(dbLib* tarLib : db_->getLibs()) {
        for(dbMaster* tarMaster : tarLib->getMasters()) {
            typeEncoder[tarMaster]=count++;
        }
    }

    // Get # of access points info.
    //unordered_map<dbInst*, int> instAccPoints_;
    //unordered_map<dbInst*, int> instBlkPoints_;
    //unordered_map<dbInst*, int> instBndPoints_;
    unordered_map<dbITerm*, int> termAccPoints;
    // Get left/right white space of inst
    //unordered_map<dbInst*, double> whiteSpaceL_;
    //unordered_map<dbInst*, double> whiteSpaceR_;
    //unordered_map<dbInst*, double> whiteSpaceT_;
    //unordered_map<dbInst*, double> whiteSpaceD_;
    // Get overlaps by PDN stripes
    //unordered_map<dbInst*, double> sWireOverlap_;
    //unordered_map<dbInst*, double> stnBBox_;
    //unordered_map<dbInst*, int> numCutEdges;
    //unordered_map<dbInst*, int> cellType_;

    ////// FOR DEBUG
    double maxWireDenDR = 0;
    double maxWireDenPL = 0;
    double avgRUDY = 0;
    for(Gcell* gcell : grid->getGcells()) {
       avgRUDY += gcell->getRUDY();
        maxWireDenPL = max(maxWireDenPL, gcell->getWireUtil(ModelType::TREE));
        maxWireDenDR = max(maxWireDenDR, gcell->getWireUtil(ModelType::ROUTE));
    }
    cout << "[REP] Max RUDY : " << grid->getMaxRUDY() << endl;
    avgRUDY /= grid->getGcells().size();
    cout << "[REP] Avg RUDY : " << avgRUDY << endl;
    cout << "[REP] Max WireDen (PL) : " << maxWireDenPL << endl;
    cout << "[REP] Max WireDen (DR) : " << maxWireDenDR << endl;
    cout << "Feature extraction finished" << endl;

    // Rtree for M1-M2 via access point
    unordered_map<dbTechLayer*, vector<int>> xGrid;
    unordered_map<dbTechLayer*, vector<int>> yGrid;
    dbTech* tech = getDb()->getTech();

    cout << "Init x-y Grid" << endl;
    for(dbTechLayer* techLayer : techLayers) {
        if(techLayer->getType() == dbTechLayerType::ROUTING) {
            dbTrackGrid* trackGrid = block->findTrackGrid(techLayer);
            trackGrid->getGridX(xGrid[techLayer]);
            trackGrid->getGridY(yGrid[techLayer]);
        }
    }
    cout << "Init trackgrid is done" << endl;

    for(dbInst* inst : block->getInsts()) {

        dbMaster* tarMaster = inst->getMaster();
        cellType_[inst] = typeEncoder[tarMaster];
        // Calculate # of points
        instAccPoints_[inst] = 0;
        instBndPoints_[inst] = 0;
        instBlkPoints_[inst] = 0;

        Rect instBBox;
        dbBox* instBox = inst->getBBox();
        instBox->getBox(instBBox);
        dbTechLayer* techM1Layer = getDb()->getTech()->findRoutingLayer(1);
        vector<int>::iterator xMinIter = lower_bound(xGrid[techM1Layer].begin(), 
                                                    xGrid[techM1Layer].end(), instBBox.xMin());
        vector<int>::iterator yMinIter = lower_bound(yGrid[techM1Layer].begin(),
                                                    yGrid[techM1Layer].end(), instBBox.yMin());
        vector<int>::iterator xMaxIter = upper_bound(xGrid[techM1Layer].begin(),
                                                    xGrid[techM1Layer].end(), instBBox.xMax());
        vector<int>::iterator yMaxIter = upper_bound(yGrid[techM1Layer].begin(),
                                                    yGrid[techM1Layer].end(), instBBox.yMax());

        // cout << inst->getName() << " (" << instBBox.xMin() << " " << instBBox.yMin() << ") (" 
        //                                << instBBox.xMax() << " " << instBBox.yMax() << ")" << endl;
       
        // cout << xMaxIter - xMinIter << " " << yMaxIter - yMinIter << endl;
        // cout << *xMinIter << " " << *yMinIter << " " << *xMaxIter << " " << *yMaxIter << endl;
       

        instBndPoints_[inst] = (yMaxIter-yMinIter) * (xMaxIter-xMinIter);
        int xOrig, yOrig;
        inst->getOrigin(xOrig, yOrig);
        dbTransform transform;
        transform.setOrient(inst->getOrient());
        transform.setOffset( Point(xOrig, yOrig) );

        int xMin, xMax, yMin, yMax;
        xMin = instBBox.xMin();
        yMin = instBBox.yMin();
        xMax = instBBox.xMax();
        yMax = instBBox.yMax();
        set<odb::dbInst*> sourceInsts;
        set<odb::dbInst*> sinkInsts;
        
        cellSize_[inst] = 1.0 * (yMax-yMin)*(xMax-xMin) / (dbu*dbu);
        isClocked_[inst] = tarMaster->isSequential();

        for(dbITerm* iterm : inst->getITerms()) {
            dbMTerm* mTerm = iterm->getMTerm();
            set<Point> accPoints;
            for(dbMPin* mPin : mTerm->getMPins()) {
                for(dbBox* pBox : mPin->getGeometry()) {
                    if(pBox->getTechLayer()->getType() != dbTechLayerType::ROUTING)
                        continue;
                    dbTechLayer* techLayer = pBox->getTechLayer();
                    Rect pinBBox;
                    pBox->getBox(pinBBox);
                    transform.apply(pinBBox);
                    vector<int>::iterator xMinIter = lower_bound(xGrid[techLayer].begin(), 
                                                                xGrid[techLayer].end(), pinBBox.xMin());
                    vector<int>::iterator xMaxIter = upper_bound(xGrid[techLayer].begin(),
                                                                xGrid[techLayer].end(), pinBBox.xMax());
                    vector<int>::iterator yMinIter = lower_bound(yGrid[techLayer].begin(),
                                                                yGrid[techLayer].end(), pinBBox.yMin());
                    vector<int>::iterator yMaxIter = upper_bound(yGrid[techLayer].begin(),
                                                                yGrid[techLayer].end(), pinBBox.yMax());
                    for(; xMinIter != xMaxIter; xMinIter++) {
                        for(; yMinIter != yMaxIter; yMinIter++) {
                            int x = *xMinIter;
                            int y = *yMinIter;
                            Point point(x,y);
                            accPoints.insert(point);
                        }
                    }
                }
            }


            int numPoints = accPoints.size();
            termAccPoints[iterm] = numPoints;
            dbNet* tarNet = iterm->getNet();
            
            if(tarNet == NULL) {
                instBlkPoints_[inst] += numPoints;
            } else {
                if(iterm->getSigType() == dbSigType::SIGNAL) {
                    instAccPoints_[inst] += numPoints;
                } else {
                    instBlkPoints_[inst] += numPoints;
                }
                // Calculate NetBBox for tarInst
                if(iterm->getIoType() == dbIoType::OUTPUT) {
                    
                    for(dbITerm* sinkITerm : tarNet->getITerms()) {
                        dbInst* sinkInst = sinkITerm->getInst();
                        if(sinkInst != inst) {
                            sinkInsts.insert(sinkInst);
                            xMin = min(xMin, sinkInst->getBBox()->xMin());
                            yMin = min(yMin, sinkInst->getBBox()->yMin());
                            xMax = min(xMax, sinkInst->getBBox()->xMax());
                            yMax = min(yMax, sinkInst->getBBox()->yMax());
                        }
                    }
                } else if(iterm->getIoType() == dbIoType::INPUT) {

                    dbITerm* sourceITerm = tarNet->getFirstOutput();
                    if(sourceITerm!=NULL) {
                        dbInst* sourceInst = sourceITerm->getInst();
                        sourceInsts.insert(sourceInst);        
                        xMin = min(xMin, sourceInst->getBBox()->xMin());
                        yMin = min(yMin, sourceInst->getBBox()->yMin());
                        xMax = min(xMax, sourceInst->getBBox()->xMax());
                        yMax = min(yMax, sourceInst->getBBox()->yMax());
                    }
                }
            }
        }


        numInEdges_[inst] = sourceInsts.size();
        numOutEdges_[inst] = sinkInsts.size();
        numEdges_[inst] = sourceInsts.size() + sinkInsts.size();

        
        stnBBox_[inst] = 1.0/(gcellWidth)*(xMax-xMin)*(yMax-yMin);

        
        bgBox tarBox(bgPoint(instBBox.xMin(), instBBox.yMin()), bgPoint(instBBox.xMax(), instBBox.yMax()));
        vector<pair<bgBox, dbInst*>> queryResults;
        
        vector<pair<bgBox, dbTechVia*>> queryViaResults;
        
        pViaRtree.query(bgi::nearest(tarBox, 2), back_inserter(queryViaResults));

        cout << "Inst Box" << " ";
        cout << instBBox.xMin() << " "
             << instBBox.yMin() << " "
             << instBBox.xMax() << " "
             << instBBox.yMax() << endl;

        int xCenInst = (instBBox.xMin()+instBBox.xMax())/2;
        int yCenInst = (instBBox.yMin()+instBBox.yMax())/2;

        for(auto& val : queryViaResults) {
            dbTechVia* via = val.second;
            bgBox viaBox = val.first;

            cout << "Via Box" << " ";
            cout << viaBox.min_corner().get<0>() << " "
                 << viaBox.min_corner().get<1>() << " "
                 << viaBox.max_corner().get<0>() << " "
                 << viaBox.max_corner().get<1>() << " ";
            
            int xCenVia = (viaBox.min_corner().get<0>()+viaBox.max_corner().get<0>())/2;
            int yCenVia = (viaBox.min_corner().get<1>()+viaBox.max_corner().get<1>())/2;

            int xDis = xCenInst - xCenVia;
            if(xDis < 0) xDis = -xDis;

            int yDis = yCenInst - yCenVia;
            if(yDis < 0) yDis = -yDis;

            int distance = xDis + yDis;

            cout << "Distance : " << distance << endl;
        }
        
        // Calculate white space (horizontal)
        xMin = instBBox.xMin() - gcellWidth;
        xMax = instBBox.xMax() + gcellWidth;
        yMin = instBBox.yMin()+1;
        yMax = instBBox.yMax()-1;

        bgBox horQueryBox(bgPoint(xMin, yMin), bgPoint(xMax,yMax));
        instRtree.query(bgi::intersects(horQueryBox), back_inserter(queryResults));
        whiteSpaceL_[inst] = 1.0;
        whiteSpaceR_[inst] = 1.0;
        for(auto& val : queryResults) {
            dbInst* adjInst = val.second;
            bgBox adjBox = val.first;

            if(inst == adjInst)
                continue;

            double dist = 1.0*bg::distance(tarBox, adjBox)/gcellWidth;
            int adjMinX = adjBox.min_corner().get<0>();
            int adjMaxX = adjBox.max_corner().get<0>();

            double adjCenX = 1.0*(adjMinX+adjMaxX)/2;
            //cout << inst->getName() << " " << adjInst->getName() << " dist " << dist << endl;
            if (adjCenX < instBBox.xMin()) {
                // leftside
                whiteSpaceL_[inst] = min(whiteSpaceL_[inst], dist);
            } else if (adjCenX > instBBox.xMax()) {
                // rightside
                whiteSpaceR_[inst] = min(whiteSpaceR_[inst], dist);
            } else {
                cout << "1 Overlapped (" << inst->getName() << " " << adjInst->getName() << ")" <<  endl;
                cout << bg::get<0,0>(tarBox) << " " << bg::get<0,1>(tarBox) << " "
                     << bg::get<1,0>(tarBox) << " " << bg::get<1,1>(tarBox) << endl;
                cout << bg::get<0,0>(adjBox) << " " << bg::get<0,1>(adjBox) << " "
                     << bg::get<1,0>(adjBox) << " " << bg::get<1,1>(adjBox) << endl;
                cout << adjCenX << endl;
                //exit(0);

            }
        }
        queryResults.clear();
        // Calculate white space (vertical)
        xMin = instBBox.xMin()+1;
        xMax = instBBox.xMax()-1;
        yMin = instBBox.yMin()-gcellHeight;
        yMax = instBBox.yMax()+gcellHeight;

        bgBox verQueryBox(bgPoint(xMin, yMin), bgPoint(xMax,yMax));
        instRtree.query(bgi::intersects(verQueryBox), back_inserter(queryResults));
        whiteSpaceT_[inst] = 1.0;
        whiteSpaceD_[inst] = 1.0;
         for(auto& val : queryResults) {
            dbInst* adjInst = val.second;
            bgBox adjBox = val.first;

            if(inst == adjInst)
                continue;

            double dist = 1.0*bg::distance(tarBox, adjBox)/gcellHeight;
            int adjMinY = adjBox.min_corner().get<1>();
            int adjMaxY = adjBox.max_corner().get<1>();
            double adjCenY = 1.0*(adjMinY+adjMaxY)/2;

            //cout << inst->getName() << " " << adjInst->getName() << " dist " << dist << endl;

            if (adjCenY < instBBox.yMin()) {
                // bottomside
                whiteSpaceD_[inst] = min(whiteSpaceD_[inst], dist);
            } else if (adjCenY > instBBox.yMax()) {
                // topside
                whiteSpaceT_[inst] = min(whiteSpaceT_[inst], dist);
            } else {
                cout << "2 Overlapped (" << inst->getName() << " " << adjInst->getName() << ")" <<  endl;
                cout << bg::get<0,0>(tarBox) << " " << bg::get<0,1>(tarBox) << " "
                     << bg::get<1,0>(tarBox) << " " << bg::get<1,1>(tarBox) << endl;
                cout << bg::get<0,0>(adjBox) << " " << bg::get<0,1>(adjBox) << " "
                     << bg::get<1,0>(adjBox) << " " << bg::get<1,1>(adjBox) << endl;
                cout << adjCenY << endl;
                //exit(0);

            }
        }

        sWireOverlap_[inst] = 0.0;
        bgBox queryBox(bgPoint(instBBox.xMin(), instBBox.yMin()), bgPoint(instBBox.xMax(), instBBox.yMax()));
        for( auto it =sWireRtree.qbegin(bgi::intersects(queryBox)); it != sWireRtree.qend(); it++) {
            xMin = bg::get<0,0>(it->first);
            yMin = bg::get<0,1>(it->first);
            xMax = bg::get<1,0>(it->first);
            yMax = bg::get<1,1>(it->first);
            Rect wireBBox(xMin, yMin, xMax, yMax);


            Rect overlap = instBBox.intersect(wireBBox);
            double area1 = overlap.dx() * overlap.dy();
            double area2 = instBBox.dx() * instBBox.dy();
            double overlapRatio = area1 / area2;
            
            sWireOverlap_[inst] += overlapRatio;
        }
    }

    
    cout << "Start to extract clip graphs" << endl; 
    for(Gcell* gcell : grid->getGcells()) {
//        cout << "Axis : " << gcell->getRow() << " " << gcell->getCol() << endl;
        set<dbInst*> instSet = gcell->getInstSet();

        //unordered_map<dbInst*, double> relPosX_;
        //unordered_map<dbInst*, double> relPosY_;
        //unordered_map<dbInst*, int> numCutEdges_;        

        Rect boundBox = gcell->getBBox();

        for(dbInst* tarInst : instSet) {
            dbBox* instBox = tarInst->getBBox();
            double cenX = 1.0*(instBox->xMin() + instBox->xMax()) / 2;
            double cenY = 1.0*(instBox->yMin() + instBox->yMax()) / 2;

            cenX = (cenX - boundBox.xMin())/boundBox.dx();
            cenY = (cenY - boundBox.yMin())/boundBox.dy();
            
            relPosX_[tarInst] = cenX;
            relPosY_[tarInst] = cenY;

            col_[tarInst] = gcell->getCol();
            row_[tarInst] = gcell->getRow();

            numCutEdges_[tarInst] = 0;

            for(dbITerm* tarITerm : tarInst->getITerms()) {
                dbNet* tarNet = tarITerm->getNet();
                dbMTerm* tarMTerm = tarITerm->getMTerm();

                if(tarNet != NULL) {
                    // Calculate numCutEdges_
//                    cout << "pin net: " << tarInst->getName() << "/" << tarMTerm->getName() << " " << tarNet->getName() << endl;
//                    cout << "IO type: " << tarITerm->getIoType() << endl;
                    
                    if(tarITerm->getIoType() == dbIoType::OUTPUT) {
                        for(dbITerm* sinkITerm : tarNet->getITerms()) {
                            dbInst* sinkInst = sinkITerm->getInst();
                            if(instSet.find(sinkInst) == instSet.end()) {
                                numCutEdges_[tarInst]++;
                            }
                        }
                    } else if(tarITerm->getIoType() == dbIoType::INPUT) {
                        dbITerm* fstITerm = tarNet->getFirstOutput();

                        if(fstITerm != NULL) {
                            
                            if(fstITerm->getIoType() != dbIoType::OUTPUT) {
                                cout << fstITerm->getInst()->getName() << " ????" << endl;
                                exit(0);
                            }

                            
                            dbInst* sourceInst = fstITerm->getInst();
                            if(instSet.find(sourceInst) == instSet.end()) {
                                numCutEdges_[tarInst]++;    
                            }
                        }
                    }
                }
            }
//            cout << endl;
        }
        Graph* instGraph = new Graph;
        instGraph->setDb(db_);
        instGraph->setSta(sta_);
        instGraph->setGraphModel(graphModel_);
        instGraph->setEdgeWeightModel(edgeWeightModel_);
        instGraph->init(instSet);

        // for timing
        instGraph->setIsCrit(isCritical_);
        instGraph->setSlack(clockPeriod_, absSlack_);
        // for pin accessibility
        instGraph->setNumPoints(instAccPoints_,
                                instBlkPoints_,
                                instBndPoints_);
        instGraph->setWhiteSpace(whiteSpaceL_,
                                 whiteSpaceR_,
                                 whiteSpaceT_,
                                 whiteSpaceD_);
        instGraph->setSWireOverlap(sWireOverlap_);
        instGraph->setCutEdges(numCutEdges_);
        instGraph->setBBoxSize(stnBBox_);
        instGraph->setCellType(cellType_);
        instGraph->setIsCrit(isCritical_); 
        instGraph->setRelPos(relPosX_, relPosY_);
        gcell->setGraph(instGraph);

   
        //instGraph->print();

        //for(dbInst* tarInst: instSet) {
        //    // For debug
        //    cout << tarInst->getName() << endl;
        //    cout << "   - relpos (x,y) : " << relPosX_[tarInst] << " " << relPosY_[tarInst] << endl;
        //    cout << "   - # acc points : " << instAccPoints_[tarInst] << endl;
        //    cout << "   - # blk points : " << instBlkPoints_[tarInst] << endl;
        //    cout << "   - # bnd points : " << instBndPoints_[tarInst] << endl;
        //    cout << "   - # cut edges  : " << numCutEdges_[tarInst] << endl;
        //    cout << "   - wspace (L/R) : " << whiteSpaceL_[tarInst] << " " << whiteSpaceR_[tarInst] << endl;
        //    cout << "   - wspace (T/D) : " << whiteSpaceT_[tarInst] << " " << whiteSpaceD_[tarInst] << endl;
        //    cout << "   - encoded type : " << cellType_[tarInst] << endl;
        //    cout << "   - size of bbox : " << stnBBox_[tarInst] << endl;
        //    cout << "   - is critical  : " << isCritical_[tarInst] << endl;
        //    cout << "   - wire overlap : " << sWireOverlap_[tarInst] << endl;
        //    cout << endl; 
        //}


    
    }
    cout << "Done!" << endl;
}








}
