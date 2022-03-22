#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include "opendb/geom.h"
#include "sta/Sta.hh"
#include "sta/Network.hh"
#include "sta/DelayFloat.hh"
#include "sta/Graph.hh"
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


using namespace odb;
using namespace std;
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::make_pair;
using std::pair;
using std::set;



namespace ClipGraphExtract {
ClipGraphExtractor::ClipGraphExtractor() : db_(nullptr), sta_(nullptr),
    numRows_(5), maxRouteLayer_(7),
    grid_(nullptr),
    graphModel_(Star), edgeWeightModel_(A), fileName_("") {};

ClipGraphExtractor::~ClipGraphExtractor() {
  clear(); 
}

void
ClipGraphExtractor::clear() {
  db_ = nullptr; 
  sta_ = nullptr;
  numRows_ =5;          // default: x5
  maxRouteLayer_ = 7;  // default: Metal7
  graphModel_ = Star;
  edgeWeightModel_ = A;
  fileName_ = "";
}



void ClipGraphExtractor::initRtree1() {


    wireRtree_ = (void*) (new SegRtree<dbNet*>);
    instRtree_ = (void*) (new BoxRtree<dbInst*>);
    SegRtree<dbNet*> *wireRtree = (SegRtree<dbNet*>*) wireRtree_;
    BoxRtree<dbInst*>* instRtree = (BoxRtree<dbInst*>*) instRtree_;
    
    // make wireRtree
    int x, y, ext;
    dbBlock* block = db_->getChip()->getBlock();
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
                        wireRtree->insert( make_pair( wireSeg, net ) );
                    }
                }
                opcode = decoder.next();
            }
        }
    }
    // make instRtree
    for( dbInst* inst : block->getInsts() ) {
        dbBox* bBox = inst->getBBox();
        bgBox b (bgPoint(bBox->xMin(), bBox->yMin()),
                    bgPoint(bBox->xMax(), bBox->yMax()));
        instRtree->insert( make_pair(b, inst) );
    }


}

void ClipGraphExtractor::initRtree2() {

    Grid* grid = (Grid*) grid_;
    rsmtRtree_ = (void*) (new SegRtree<RSMT*>);
    gcellRtree_ = (void*) (new BoxRtree<Gcell*>);
    BoxRtree<Gcell*>* gcellRtree = (BoxRtree<Gcell*>*) gcellRtree_;
    SegRtree<RSMT*>* rsmtRtree = (SegRtree<RSMT*>*) rsmtRtree_;

 
    // make gcellRtree
    for( Gcell* gcell : grid->getGcells() ) {
        bgBox gcellBox = gcell->getQueryBox();
        gcellRtree->insert( make_pair(gcellBox, gcell) );
    }

  
    // make rsmtRtree
    for( dbNet* net : db_->getChip()->getBlock()->getNets()) {
        RSMT* rsmt = grid->getRSMT(net);
        vector<Rect> segments = rsmt->getSegments();
        // insert segments into rtree
        for(Rect& seg : segments) {
            // update (1) #cut-nets (2) wire utilization
            bgSeg bgseg(bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()));
            rsmtRtree->insert( make_pair( bgseg, rsmt ) );
        }
    }


}

void
ClipGraphExtractor::init() {
    sta_->updateTiming(false);
    Flute::readLUT();


    initGrid();
    cout << "InitGrid is done" << endl;
   

    initGraph();
    cout << "InitGraph is done" << endl;

}

void ClipGraphExtractor::initGrid() {
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


    cout << "TrackSupply    : " << trackSupply << endl;
    cout << "WireCapacity   : " << wireCapacity << endl;

    // Get core area
    odb::Rect blockArea;
    block->getBBox()->getBox(blockArea);

    // Initialize Gcell Grid
    grid->setDb(getDb());
    grid->setBoundary(blockArea);
    grid->setGcellWidth(gcellWidth);
    grid->setGcellHeight(gcellHeight);
    grid->setWireCapacity(wireCapacity);
    grid->setTrackSupply(trackSupply);
    grid->setNumLayers(numLayer);
    grid->setWireMinWidth(minWidth);
    grid->init();
    //  

    cout << "Grid initialization finished" << endl;


    initRtree1();
    // has to finish grid->init() before calling initRtree()
    initRtree2();

    SegRtree<dbNet*> *wireRtree = (SegRtree<dbNet*>*) wireRtree_;
    BoxRtree<Gcell*>* gcellRtree = (BoxRtree<Gcell*>*) gcellRtree_;
    BoxRtree<dbInst*>* instRtree = (BoxRtree<dbInst*>*) instRtree_;
    SegRtree<RSMT*>* rsmtRtree = (SegRtree<RSMT*>*) rsmtRtree_;
    
    //
    for( dbNet* net : block->getNets() ) {
        if(net->isSpecial()) {
        }
        // search overlapping gcells
        RSMT* rsmt = grid->getRSMT(net);
        rsmt->searchOverlaps(gcellRtree);
    }

    cout << "RSMT construction finished" << endl;

    // add Instance in gcell 
    for( Gcell* gcell : grid->getGcells() ) {
        gcell->extractPlaceFeature(instRtree);
        gcell->extractPlaceFeature(rsmtRtree);
		gcell->extractRouteFeature(wireRtree);
    }
    
    
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
}


void ClipGraphExtractor::initGraph() {

    sta_->updateTiming(false);
    sta::dbNetwork* network = sta_->getDbNetwork();
    sta::Graph* graph = sta_->ensureGraph();
    Grid* grid = (Grid*)grid_;

    dbBlock* block = getDb()->getChip()->getBlock();

    sta::Vertex* wstVertex;
    sta::Slack wstSlack, totNegSlack;
    sta_->worstSlack(sta::MinMax::max(), wstSlack, wstVertex);
    totNegSlack = sta_->totalNegativeSlack(sta::MinMax::max());

    cout << "Worst Slack    : " << wstSlack << endl;
    cout << "Total NegSlack : " << totNegSlack << endl;



    // Get slack info.
    unordered_map<dbInst*, double> minSlack;

    for(dbInst* inst : block->getInsts()) {
    
        minSlack[inst] = DBL_MAX;

        for(dbITerm* iterm : inst->getITerms()) {

            if(iterm->getIoType() == dbSigType::POWER || iterm->getIoType() == dbSigType::GROUND)
                continue;

            uint32_t vertexId = iterm->staVertexId();
            sta::Vertex* vertex = graph->vertex(vertexId);
            if(vertex == NULL)
                continue;

            sta::Slack staSlk = sta_->vertexSlack(vertex, sta::MinMax::max());
            //sta::Slack minSlack = sta_->vertexSlack(vertex, sta::MinMax::min());
        
            double slack = (double)sta::delayAsFloat(staSlk);
            minSlack[inst] = min(minSlack[inst], slack);
        }

        //cout << inst->getName() << "'s min slack is " << minSlack[inst] << endl;

    }


    // Rtree for M1-M2 via access point
    dbSet<dbTechLayer> techLayers = getDb()->getTech()->getLayers();
    unordered_map<dbTechLayer*, vector<int>> xGrid;
    unordered_map<dbTechLayer*, vector<int>> yGrid;
    dbTech* tech = getDb()->getTech();

    cout << "Init x-y Grid" << endl;
    for(dbTechLayer* techLayer : tech->getLayers()) {
        if(techLayer->getType() == dbTechLayerType::ROUTING) {
            dbTrackGrid* trackGrid = block->findTrackGrid(techLayer);
            trackGrid->getGridX(xGrid[techLayer]);
            trackGrid->getGridY(yGrid[techLayer]);
        }
    }
    cout << "Init trackgrid is done" << endl;
    /*
    dbTechLayer* techLayer1 = getDb()->getTech()->findRoutingLayer(1);
    dbTechLayer* techLayer2 = getDb()->getTech()->findRoutingLayer(2);
    if(techLayer1->getDirection() == techLayer2->getDirection()) {
        cout << "??" << endl;
        exit(0);
    }
    if(techLayer1->getDirection() == dbTechLayerDir::HORIZONTAL &&
            techLayer2->getDirection() == TechLayerDir::VERTICAL) {
        techLayer1->getTrackGrid()->getGridY(yGrid);
        techLayer2->getTrackGrid()->getGridX(xGrid);

    } else if(techLayer2->getDirection() == dbTechLayerDir::HORIZONTAL &&
            techLayer1->getDirection() == TechLayerDir::VERTICAL) {
        techLayer2->getTrackGrid()->getGridY(yGrid);
        techLayer1->getTrackGrid()->getGridX(xGrid);
    } else {
        cout << "exception case!" << endl;
        exit(0);
    }
    */
    
    // Get # of access points info.
    unordered_map<dbInst*, int> instAccPoints;
    unordered_map<dbInst*, int> instBlkPoints;
    unordered_map<dbInst*, int> instBndPoints;

    unordered_map<dbITerm*, int> termAccPoints;


    for(dbInst* inst : block->getInsts()) {
        instAccPoints[inst] = 0;

        Rect instBBox;
        dbBox* instBox = inst->getBBox();
        instBox->getBox(instBBox);


        dbTechLayer* techM1Layer = getDb()->getTech()->findRoutingLayer(1);

        vector<int>::iterator xMinIter = lower_bound(xGrid[techM1Layer].begin(), xGrid[techM1Layer].end(), instBBox.xMin());
        vector<int>::iterator xMaxIter = upper_bound(xGrid[techM1Layer].begin(), xGrid[techM1Layer].end(), instBBox.xMax());
        vector<int>::iterator yMinIter = lower_bound(yGrid[techM1Layer].begin(), yGrid[techM1Layer].end(), instBBox.yMin());
        vector<int>::iterator yMaxIter = upper_bound(yGrid[techM1Layer].begin(), yGrid[techM1Layer].end(), instBBox.yMin());
        instBndPoints[inst] = (yMaxIter-yMinIter) * (xMaxIter-xMinIter);


        int xOrig, yOrig;
        inst->getOrigin(xOrig, yOrig);
        dbTransform transform;
        transform.setOrient(inst->getOrient());
        transform.setOffset( Point(xOrig, yOrig) );

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
                    
                    int xMin = pinBBox.xMin();
                    int yMin = pinBBox.yMin();
                    int xMax = pinBBox.xMax();
                    int yMax = pinBBox.yMax();


                    vector<int>::iterator xMinIter = lower_bound(xGrid[techLayer].begin(), xGrid[techLayer].end(), xMin);
                    vector<int>::iterator xMaxIter = upper_bound(xGrid[techLayer].begin(), xGrid[techLayer].end(), xMax);
                    vector<int>::iterator yMinIter = lower_bound(yGrid[techLayer].begin(), yGrid[techLayer].end(), yMin);
                    vector<int>::iterator yMaxIter = upper_bound(yGrid[techLayer].begin(), yGrid[techLayer].end(), yMax);
                    
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



            dbNet* net = iterm->getNet();
            if(net == NULL) {
                instBlkPoints[inst] += numPoints;
            } else {
                if(net->getSigType() == dbSigType::POWER) {
                    instBlkPoints[inst] += numPoints; break;
                } else if(net->getSigType() == dbSigType::GROUND) {
                    instBlkPoints[inst] += numPoints; break;
                } else {
                    instAccPoints[inst] += numPoints; break;
                }
            }
        }
    }


    // Get left/right white space of inst
    unordered_map<dbInst*, int> whiteSpaceL;
    unordered_map<dbInst*, int> whiteSpaceR;



    for(Gcell* gcell : grid->getGcells()) {

        set<dbInst*> instSet = gcell->getInstSet();
        Graph* instGraph = new Graph;
        instGraph->setDb(db_);
        instGraph->setSta(sta_);
        instGraph->setGraphModel(graphModel_);
        instGraph->setEdgeWeightModel(edgeWeightModel_);
        instGraph->init(instSet);

        // for timing
        instGraph->setMinSlack(minSlack);
        // for pin accessibility
        instGraph->setNumAccPoints(instAccPoints);
        instGraph->setNumBndPoints(instBndPoints);
        instGraph->setNumBlkPoints(instBlkPoints);
        instGraph->setWhiteSpaceL(whiteSpaceL);
        instGraph->setWhiteSpaceR(whiteSpaceR);
        gcell->setGraph(instGraph);
    }
    cout << "Done!" << endl;
}



void
ClipGraphExtractor::extract(int lx, int ly, int ux, int uy) {

}


void
ClipGraphExtractor::setDb(odb::dbDatabase* db) {
  db_ = db;
}

void
ClipGraphExtractor::setSta(sta::dbSta* sta) {
  sta_ = sta;
}


void ClipGraphExtractor::setGcellSize(int numRows) {
    numRows_ = numRows;
}

void ClipGraphExtractor::setMaxRouteLayer(int maxRouteLayer) {
    maxRouteLayer_ = maxRouteLayer;
}



void
ClipGraphExtractor::setGraphModel(const char* graphModel) {
  if( strcmp(graphModel, "star") == 0 ) {
    graphModel_ = Star; 
  }
  else if( strcmp(graphModel, "clique") == 0 ) {
    graphModel_ = Clique;
  }
}

void
ClipGraphExtractor::setSaveFileName(const char* fileName) {
  fileName_ = fileName;
}

void
ClipGraphExtractor::setSaveFilePrefix(const char* prefix) {
    prefix_ = prefix;
}


void
ClipGraphExtractor::setEdgeWeightModel( const char* edgeWeightModel ) {
  if( strcmp(edgeWeightModel, "a") == 0 ) {
    edgeWeightModel_ = A;
  }
  else if( strcmp(edgeWeightModel, "b") == 0 ) {
    edgeWeightModel_ = B;
  }
  else if( strcmp(edgeWeightModel, "c") == 0 ) {
    edgeWeightModel_ = C;
  }
  else if( strcmp(edgeWeightModel, "d") == 0 ) {
    edgeWeightModel_ = D;
  }
  else if( strcmp(edgeWeightModel, "e") == 0 ) {
    edgeWeightModel_ = E;
  }
  else {
    cout << "ERROR: edgeWeight is wrong: " << edgeWeightModel << endl;
    exit(1);
  }
}

}


  //inst_RTree* inst_rTree = (inst_RTree*) inst_rTree_;
    //// Print worst slack
    //unordered_map<dbInst*, float> slack;


        //for(dbInst* inst : block->getInsts()) {

    //    for(dbITerm* iterm : inst->getITerms()) {

    //        if(iterm->getIoType() == dbSigType::POWER ||
    //                iterm->getIoType() == dbSigType::GROUND)
    //            continue;

    //        uint32_t vertexId = iterm->staVertexId();
    //        sta::Vertex* vertex = graph->vertex(vertexId);
    //        sta::Slack maxSlack = sta_->vertexSlack(vertex, sta::MinMax::max());
    //        sta::Slack minSlack = sta_->vertexSlack(vertex, sta::MinMax::min());

    //        if(inst->getName() == "f_permutation__out_reg_902_") {
    //            cout << iterm->getMTerm()->getName() << " " <<  maxSlack << " " << minSlack << endl;
    //        }
    //    }
    //}

    //cout << "Done" << endl;


    //for(dbNet* net : block->getNets()) {
    //    
    //    for(dbITerm* iterm : net->getITerms()) {
    //        uint32_t vertexId = iterm->staVertexId();
    //        //cout << iterm->getName() << " " << vertexId << endl;
    //    }
    //}
