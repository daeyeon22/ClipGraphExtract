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



void initWireRtree(odb::dbDatabase* db, SegRtree<dbNet*> &rwireRtree, BoxRtree<dbNet*> &swireRtree) {

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
                switch(opcode) {
                    case dbWireDecoder::PATH:
                        points.clear(); break;
                    case dbWireDecoder::JUNCTION: break;
                    case dbWireDecoder::SHORT: break;

                    case dbWireDecoder::TECH_VIA: {
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
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
         
                if(hasPath && points.size() > 1) {
                    odb::Point pt1 = points[points.size()-2];
                    odb::Point pt2 = points[points.size()-1];

                    int xMin = min(pt1.getX(), pt2.getX());// - layerMinWidth[decoder.getLayer()]/2;
                    int xMax = max(pt1.getX(), pt2.getX());// + layerMinWidth[decoder.getLayer()]/2;
                    int yMin = min(pt1.getY(), pt2.getY());// - layerMinWidth[decoder.getLayer()]/2;
                    int yMax = max(pt1.getY(), pt2.getY());// + layerMinWidth[decoder.getLayer()]/2;
                    int dist = (xMax-xMin) + (yMax-yMin);
                    
                    if(dist != 0) {
                        wl_ += dist;
                        bgSeg wireSeg( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                        rwireRtree.insert( make_pair( wireSeg, net ) );
                    }
                }
                opcode = decoder.next();
            }
        }
        
        for(dbSWire* swire : net->getSWires()) {
            
            for(dbSBox* sbox : swire->getWires()) {
                int xMin = sbox->xMin();
                int yMin = sbox->yMin();
                int xMax = sbox->xMax();
                int yMax = sbox->yMax();


                bgBox wireBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                swireRtree.insert( make_pair( wireBox, net ) );
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
        double &clockPeriod,
        unordered_map<dbInst*, double> &absoluteSlack, 
        unordered_map<dbInst*, bool> &isCrit) {

    sta_->updateTiming(false);
    sta::dbNetwork* network = sta_->getDbNetwork();
    sta::Graph* graph = sta_->ensureGraph();

    sta::Vertex* wstVertex;
    sta::Slack wstSlack, totNegSlack;
    sta_->worstSlack(sta::MinMax::max(), wstSlack, wstVertex);
    totNegSlack = sta_->totalNegativeSlack(sta::MinMax::max());
    cout << "Worst Slack    : " << wstSlack << endl;
    cout << "Total NegSlack : " << totNegSlack << endl;

    double wstSlk = (double)sta::delayAsFloat(wstSlack);


    string clockName = "clk";
    const sta::Clock* staClock = sta_->findClock(clockName.c_str());
    if(staClock == NULL) {
        cout << "Cannot find clock (" << clockName << ")" << endl;
        clockPeriod = 1.0;
    } else {
        cout << staClock->name() << " -period " << staClock->period() << endl;
        clockPeriod = staClock->period();
    }

    dbBlock* block = db_->getChip()->getBlock();
    // Get slack info.
    for(dbInst* inst : block->getInsts()) {
        absoluteSlack[inst] = 0.0;
        for(dbITerm* iterm : inst->getITerms()) {
            if(iterm->getIoType() == dbSigType::POWER || iterm->getIoType() == dbSigType::GROUND)
                continue;
            uint32_t vertexId = iterm->staVertexId();
            sta::Vertex* vertex = graph->vertex(vertexId);
            if(vertex == NULL)
                continue;
            sta::Slack staSlk = sta_->vertexSlack(vertex, sta::MinMax::max());
            double slack = (double)sta::delayAsFloat(staSlk);
            absoluteSlack[inst] = min(absoluteSlack[inst], slack);
        }

        isCrit[inst] = (absoluteSlack[inst] < 0.9* wstSlk)? true : false;
        //
        if(isCrit[inst]) {
            cout << inst->getName() << " is in a critical path (" << absoluteSlack[inst] << ")"  << endl;
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

            trackSupply+= supply;
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

    // Init Rtrees
    SegRtree<dbNet*> rWireRtree;
    BoxRtree<dbNet*> sWireRtree;
    BoxRtree<dbInst*> instRtree;
    BoxRtree<Gcell*> gcellRtree;
    SegRtree<RSMT*> rsmtRtree;
    initWireRtree(db_, rWireRtree, sWireRtree);
    initInstRtree(db_, instRtree);
    initGcellRtree((Grid*)grid_, gcellRtree);
    initRsmtRtree((Grid*)grid_, rsmtRtree);
    // 
    //
    unordered_map<dbInst*, bool> isCrit;
    unordered_map<dbInst*, double> absoluteSlack;
    double clockPeriod;
    // Get timing info.
    getTimingInfo(db_, sta_, clockPeriod, absoluteSlack, isCrit);


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
        gcell->updateTimingInfo(absoluteSlack);
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
    unordered_map<dbInst*, int> instAccPoints;
    unordered_map<dbInst*, int> instBlkPoints;
    unordered_map<dbInst*, int> instBndPoints;
    unordered_map<dbITerm*, int> termAccPoints;
    // Get left/right white space of inst
    unordered_map<dbInst*, double> whiteSpaceL;
    unordered_map<dbInst*, double> whiteSpaceR;
    unordered_map<dbInst*, double> whiteSpaceT;
    unordered_map<dbInst*, double> whiteSpaceD;

    // Get overlaps by PDN stripes
    unordered_map<dbInst*, double> sWireOverlap;

    unordered_map<dbInst*, double> bboxSize;
    unordered_map<dbInst*, int> numCutEdges;
    unordered_map<dbInst*, int> cellType;


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
        cellType[inst] = typeEncoder[tarMaster];
        
        // Calculate # of points
        instAccPoints[inst] = 0;
        instBndPoints[inst] = 0;
        instBlkPoints[inst] = 0;

        Rect instBBox;
        dbBox* instBox = inst->getBBox();
        instBox->getBox(instBBox);
        dbTechLayer* techM1Layer = getDb()->getTech()->findRoutingLayer(1);
        vector<int>::iterator xMinIter = lower_bound(xGrid[techM1Layer].begin(), xGrid[techM1Layer].end(), instBBox.xMin());
        vector<int>::iterator xMaxIter = upper_bound(xGrid[techM1Layer].begin(), xGrid[techM1Layer].end(), instBBox.xMax());
        vector<int>::iterator yMinIter = lower_bound(yGrid[techM1Layer].begin(), yGrid[techM1Layer].end(), instBBox.yMin());
        vector<int>::iterator yMaxIter = upper_bound(yGrid[techM1Layer].begin(), yGrid[techM1Layer].end(), instBBox.yMax());
        //cout << inst->getName() << " (" << instBBox.xMin() << " " << instBBox.yMin() << ") (" << instBBox.xMax() << " " << instBBox.yMax() << ")" << endl;
       
        //cout << xMaxIter - xMinIter << " " << yMaxIter - yMinIter << endl;
        instBndPoints[inst] = (yMaxIter-yMinIter) * (xMaxIter-xMinIter);
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
                    vector<int>::iterator xMinIter = lower_bound(xGrid[techLayer].begin(), xGrid[techLayer].end(), pinBBox.xMin());
                    vector<int>::iterator xMaxIter = upper_bound(xGrid[techLayer].begin(), xGrid[techLayer].end(), pinBBox.xMax());
                    vector<int>::iterator yMinIter = lower_bound(yGrid[techLayer].begin(), yGrid[techLayer].end(), pinBBox.yMin());
                    vector<int>::iterator yMaxIter = upper_bound(yGrid[techLayer].begin(), yGrid[techLayer].end(), pinBBox.yMax());
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
                instBlkPoints[inst] += numPoints;
            } else {
                
                if(iterm->getSigType() == dbSigType::SIGNAL) {
                    instAccPoints[inst] += numPoints;
                } else {
                    instBlkPoints[inst] += numPoints;
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

        
        bboxSize[inst] = 1.0/(gcellWidth)*(xMax-xMin)*(yMax-yMin);

        
        bgBox tarBox(bgPoint(instBBox.xMin(), instBBox.yMin()), bgPoint(instBBox.xMax(), instBBox.yMax()));
        vector<pair<bgBox, dbInst*>> queryResults;
        // Calculate white space (horizontal)
        xMin = instBBox.xMin() - gcellWidth;
        xMax = instBBox.xMax() + gcellWidth;
        yMin = instBBox.yMin()+1;
        yMax = instBBox.yMax()-1;

        bgBox horQueryBox(bgPoint(xMin, yMin), bgPoint(xMax,yMax));
        instRtree.query(bgi::intersects(horQueryBox), back_inserter(queryResults));
        whiteSpaceL[inst] = 1.0;
        whiteSpaceR[inst] = 1.0;
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
                whiteSpaceL[inst] = min(whiteSpaceL[inst], dist);
            } else if (adjCenX > instBBox.xMax()) {
                // rightside
                whiteSpaceR[inst] = min(whiteSpaceR[inst], dist);
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
        // Calculate white space (horizontal)
        xMin = instBBox.xMin()+1;
        xMax = instBBox.xMax()-1;
        yMin = instBBox.yMin()-gcellHeight;
        yMax = instBBox.yMax()+gcellHeight;

        bgBox verQueryBox(bgPoint(xMin, yMin), bgPoint(xMax,yMax));
        instRtree.query(bgi::intersects(verQueryBox), back_inserter(queryResults));
        whiteSpaceT[inst] = 1.0;
        whiteSpaceD[inst] = 1.0;
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
                whiteSpaceD[inst] = min(whiteSpaceD[inst], dist);
            } else if (adjCenY > instBBox.yMax()) {
                // topside
                whiteSpaceT[inst] = min(whiteSpaceT[inst], dist);
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

        sWireOverlap[inst] = 0.0;
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
            
            sWireOverlap[inst] += overlapRatio;
        }
    }

    
    cout << "Start to extract clip graphs" << endl; 
    for(Gcell* gcell : grid->getGcells()) {


        set<dbInst*> instSet = gcell->getInstSet();

        unordered_map<dbInst*, double> relPosX;
        unordered_map<dbInst*, double> relPosY;
        unordered_map<dbInst*, int> cutEdges;        

        Rect boundBox = gcell->getBBox();

        for(dbInst* tarInst : instSet) {
            dbBox* instBox = tarInst->getBBox();
            double cenX = 1.0*(instBox->xMin() + instBox->xMax()) /2;
            double cenY = 1.0*(instBox->yMin() + instBox->yMax()) /2;

            cenX = (cenX - boundBox.xMin())/boundBox.dx();
            cenY = (cenY - boundBox.yMin())/boundBox.dy();
            relPosX[tarInst] = cenX;
            relPosY[tarInst] = cenY;

            cutEdges[tarInst] = 0;
            for(dbITerm* tarITerm : tarInst->getITerms()) {
                dbNet* tarNet = tarITerm->getNet();
                if(tarNet != NULL) {
                    // Calculate cutEdges
                    if(tarITerm->getIoType() == dbIoType::OUTPUT) {
                        for(dbITerm* sinkITerm : tarNet->getITerms()) {
                            dbInst* sinkInst = sinkITerm->getInst();
                            if(instSet.find(sinkInst) == instSet.end()) {
                                cutEdges[tarInst]++;
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
                                cutEdges[tarInst]++;    
                            }
                        }
                    }
                }
            }
        }
        Graph* instGraph = new Graph;
        instGraph->setDb(db_);
        instGraph->setSta(sta_);
        instGraph->setGraphModel(graphModel_);
        instGraph->setEdgeWeightModel(edgeWeightModel_);
        instGraph->init(instSet);

        // for timing
        instGraph->setIsCrit(isCrit);
        instGraph->setSlack(clockPeriod, absoluteSlack);
        // for pin accessibility
        instGraph->setNumPoints(instAccPoints,
                                instBlkPoints,
                                instBndPoints);
        instGraph->setWhiteSpace(whiteSpaceL,
                                 whiteSpaceR,
                                 whiteSpaceT,
                                 whiteSpaceD);
        instGraph->setSWireOverlap(sWireOverlap);
        instGraph->setCutEdges(cutEdges);
        instGraph->setBBoxSize(bboxSize);
        instGraph->setCellType(cellType);
        instGraph->setIsCrit(isCrit); 
        instGraph->setRelPos(relPosX, relPosY);
        gcell->setGraph(instGraph);

   
        //instGraph->print();

        //for(dbInst* tarInst: instSet) {
        //    // For debug
        //    cout << tarInst->getName() << endl;
        //    cout << "   - relpos (x,y) : " << relPosX[tarInst] << " " << relPosY[tarInst] << endl;
        //    cout << "   - # acc points : " << instAccPoints[tarInst] << endl;
        //    cout << "   - # blk points : " << instBlkPoints[tarInst] << endl;
        //    cout << "   - # bnd points : " << instBndPoints[tarInst] << endl;
        //    cout << "   - # cut edges  : " << cutEdges[tarInst] << endl;
        //    cout << "   - wspace (L/R) : " << whiteSpaceL[tarInst] << " " << whiteSpaceR[tarInst] << endl;
        //    cout << "   - wspace (T/D) : " << whiteSpaceT[tarInst] << " " << whiteSpaceD[tarInst] << endl;
        //    cout << "   - encoded type : " << cellType[tarInst] << endl;
        //    cout << "   - size of bbox : " << bboxSize[tarInst] << endl;
        //    cout << "   - is critical  : " << isCrit[tarInst] << endl;
        //    cout << "   - wire overlap : " << sWireOverlap[tarInst] << endl;
        //    cout << endl; 
        //}


    
    }
    cout << "Done!" << endl;
}








}
