#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include "opendb/geom.h"
#include "sta/Sta.hh"
#include "sta/Network.hh"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "instGraph.h"
#include "grid.h"

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



void
ClipGraphExtractor::init() {
    sta_->updateTiming(false);

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
    // init rtree
    BoxRtree<dbInst*> instRtree;
    BoxRtree<Gcell*> gcellRtree;
    SegRtree<dbNet*> drRtree;
    SegRtree<RSMT*> rsmtRtree;
    
    // make gcellRtree
    for( Gcell* gcell : grid->getGcells() ) {
        bgBox gcellBox = gcell->getQueryBox();
        gcellRtree.insert( make_pair(gcellBox, gcell) );
    }

    // make instRtree
    for( dbInst* inst : block->getInsts() ) {
        dbBox* bBox = inst->getBBox();
        bgBox b (bgPoint(bBox->xMin(), bBox->yMin()),
                    bgPoint(bBox->xMax(), bBox->yMax()));
        instRtree.insert( make_pair(b, inst) );
    }

    // init FLUTE
    Flute::readLUT();
    // wireRtree (eGR result)
    

    for( dbNet* net : block->getNets() ) {
		
		// Segmentation fault occurs because Net is not entered here.
        if(net->isSpecial()) {
//          cout << net->getName() << " is SpecialNet" << endl;
//			continue;
        }


       
        RSMT* myRSMT = grid->createRSMT(net);
        vector<pair<bgBox, Gcell*>> queryResults;
        vector<odb::Rect> segments = myRSMT->getSegments();
        // insert segments into rtree
        for(odb::Rect& seg : segments) {
            // update (1) #cut-nets (2) wire utilization
            bgSeg bgseg(bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()));
            rsmtRtree.insert( make_pair( bgseg, myRSMT ) );
        }


        // make wireRtree
        dbWire* wire = net->getWire();
       
        if( wire && wire->length() ) {

            //cout << net->getName() << " has wire" << endl;
            int wl = wire->getLength();
            int wl_ = 0;
            int x,y, ext;
            
            dbWireDecoder decoder;
            decoder.begin(wire);
            
            vector<odb::Point> points;

            dbWireDecoder::OpCode opcode = decoder.peek();

            while(opcode != dbWireDecoder::END_DECODE) {
                bool hasPath=false;

                switch(opcode) {
                    case dbWireDecoder::PATH: {
                        //cout << "PATH" << endl; 
                        points.clear(); break;
                    }
                    case dbWireDecoder::JUNCTION: {
                        //cout << "JUNcTION" << endl; 
                        break;
                    }
                    case dbWireDecoder::SHORT: {
                        //cout << "SHORT" << endl;  points.clear();            
                        break;
                    }
                    case dbWireDecoder::TECH_VIA: {
                        //cout << "TECH_VIA" << endl;
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        hasPath=true;
                        break;
                    }
                    case dbWireDecoder::VIA: 
                        //cout << "VIA" << endl; 
                        break;
                    case dbWireDecoder::VWIRE:
                        //cout << "VWIRE" << endl; 
                        break;
                    case dbWireDecoder::POINT_EXT:
                        decoder.getPoint(x,y,ext);
                        //cout << "POINT_EXT " << ext << endl; 
                        points.push_back(odb::Point(x,y));
                        hasPath=true;
                        break;

                    case dbWireDecoder::BTERM:
                        //cout <<"BTERM" << endl; 
                        break;
                    case dbWireDecoder::ITERM:
                        //cout <<"ITERM" << endl; 
                        break;
                    case dbWireDecoder::RULE:
                        //cout << "RULE" << endl; 
                        break;
                    case dbWireDecoder::POINT: {
                        ////cout << "POINT" << endl;
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        hasPath=true;

                        break;
                    }
                    default:
                        //cout << "DEFAULT" << endl; 
                        break;
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
                        //if (xMin != xMax && yMin != yMax)
                        //    cout <<"(" <<  xMin << " " << yMin << ") (" << xMax << " " << yMax << ")**" << endl;
                        //else
                        //    cout <<"(" <<  xMin << " " << yMin << ") (" << xMax << " " << yMax << ") " << dist << endl;

                        bgSeg wireSeg( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                        drRtree.insert( make_pair( wireSeg, net ) );
                    }
                }

                opcode = decoder.next();
            }


            //if(wl != wl_) {
            //    cout << "different wirelength..." << endl;
            //    cout << wl << " " << wl_ << endl;
            //    exit(0);
            //}


        }
        // search overlapping gcells
        myRSMT->searchOverlaps(gcellRtree);
    }

    cout << "dr rtree size : " << drRtree.size() << endl;

    cout << "RSMT construction finished" << endl;

    // add Instance in gcell 
    for( Gcell* gcell : grid->getGcells() ) {
        
        gcell->extractFeaturePL(instRtree);
        gcell->extractFeatureRSMT(rsmtRtree);
		gcell->extractFeatureDR(drRtree);

        //if(gcell->getNumMarkers() > 0)
            // for debug
            //gcell->print();
    }
    ////// FOR DEBUG
    cout << "[REP] Max RUDY : " << grid->getMaxRUDY() << endl;

    
    double maxWireDenDR = 0;
    double maxWireDenPL = 0;
    double avgRUDY = 0;
    for(Gcell* gcell : grid->getGcells()) { 
       avgRUDY += gcell->getRUDY();
        maxWireDenPL = max(maxWireDenPL, gcell->getWireDensity(ModelType::PL));
        maxWireDenDR = max(maxWireDenDR, gcell->getWireDensity(ModelType::DR));
    }
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

    for(Gcell* gcell : grid->getGcells()) {

        set<dbInst*> instSet = gcell->getInstSet();
        Graph* instGraph = new Graph;
        instGraph->setDb(db_);
        instGraph->setSta(sta_);
        instGraph->setGraphModel(graphModel_);
        instGraph->setEdgeWeightModel(edgeWeightModel_);
        instGraph->init(instSet);
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
    //sta::Vertex* wstVertex;
    //sta::Slack wstSlack;
    //unordered_map<dbInst*, float> slack;


    //sta_->worstSlack(sta::MinMax::max(), wstSlack, wstVertex);
    //cout << "Worst slack : " << wstSlack << endl;

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
