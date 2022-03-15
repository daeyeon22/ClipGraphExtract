#include "opendb/geom.h"
#include "grid.h"
#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>


using namespace std;
using namespace odb;
using namespace feature_extractor;
using namespace ClipGraphExtract;


void
ClipGraphExtractor::saveMapImages( const char* imgDir ) {
    Grid* grid = (Grid*)(grid_);
    grid->saveMapImages(string(imgDir));
}



void ClipGraphExtractor::initGcellGrid(int numRows, int maxLayer) {
    grid_ = (void*) new Grid();
    Grid* grid = (Grid*) grid_;
    dbBlock* block = getDb()->getChip()->getBlock();
    
    // Calculate grid size = rowHeight * numRows
    dbSite* site = block->getRows().begin()->getSite(); 

    int gcellWidth = site->getHeight() * numRows;
    int gcellHeight = site->getHeight() * numRows;

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

            if(numLayer == maxLayer)
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


namespace feature_extractor {

void Grid::init() {
    assert(gcellWidth_ == 0);
    assert(gcellHeight_ == 0);
    assert(bbox_.dx() * bbox_.dy() == 0);


    cout << "Initialize gcell grid (" << gcellWidth_ << " " << gcellHeight_ << ")" << endl;


    numCols_ = bbox_.dx() / gcellWidth_;
    numRows_ = bbox_.dy() / gcellHeight_;

    int x1, y1, x2, y2;
    for(x1=0; x1 < bbox_.xMax()-gcellWidth_; x1+= gcellWidth_) {
        for(y1=0; y1 < bbox_.yMax()-gcellHeight_; y1+=gcellHeight_) {
            x2 = min(bbox_.xMax(), x1 + gcellWidth_);
            y2 = min(bbox_.yMax(), y1 + gcellHeight_);
            Gcell* gcell = createGcell(x1,y1,x2,y2);
            gcell->setTrackSupply(trackSupply_);
            gcell->setWireCapacity(wireCapacity_);
            gcell->setNumLayers(numLayers_);
        }
    }
}


double Grid::getMaxRUDY() {
    double maxRUDY = 0;
    for(Gcell* gcell : gcells_) 
        maxRUDY = max(maxRUDY, gcell->getRUDY());

    return maxRUDY;
}



vector<Gcell*> Grid::getGcells() {
    return gcells_;
}

Gcell* Grid::createGcell(int x1, int y1, int x2, int y2) {
    Gcell* gcell = new Gcell();
    gcell->setBoundary(Rect(x1,y1, x2,y2));
    gcells_.push_back(gcell);
    return gcell;
}


RSMT* Grid::getRSMT(odb::dbNet* net) {
    return net2rsmt_[net];
}


RSMT* Grid::createRSMT(odb::dbNet* net) {
    RSMT* myRSMT = new RSMT(net);
    dbSet<dbITerm> iterms = net->getITerms();

    //
    net2rsmt_[net] = myRSMT;


    // add terminals
    int x,y;
    for(dbITerm* iterm : net->getITerms()) {
        iterm->getAvgXY(&x, &y);
        myRSMT->addTerminal(x,y);
    }
    for(dbBTerm* bterm : net->getBTerms()) {
        //cout << bterm->getName() << endl;
        if(bterm->getFirstPinLocation(x,y)) {
            myRSMT->addTerminal(x,y);
        }
    }

    //cout << "#Terminals : " << myRSMT->getNumTerminals() << endl;

    // create RSMT
    myRSMT->createTree();
    myRSMT->setWireWidth(minWidth_);


    // DEBUG
    //double w_den = myRSMT->getWireUniformDensity();
    //cout << net->getName() << endl;
    //cout << "   - wire length (RSMT) : " << myRSMT->getWireLengthRSMT() << endl;
    //cout << "   - wire area (RSMT)   : " << myRSMT->getWireLengthRSMT() * minWidth_ << endl;
    //cout << "   - bbox area          : " << myRSMT->getBBox().area() << endl;
    //cout << "   - wire uniform den   : " << w_den << endl;

    //if(myRSMT->getBBox().area() == 0){
    //    cout << "BBox is 0" << endl;
    //}
    assert(w_den <0 || w_den > 1);
    rsmts_.push_back(myRSMT);

    return myRSMT;
}


void Grid::setWireMinWidth(int width) {
    minWidth_ = width;
}


void Grid::setDb(dbDatabase* db) {
    db_ = db;
}
void Grid::setBoundary(odb::Rect rect) {
    bbox_ = rect;
}

void Grid::setGcellWidth(int width) {
    gcellWidth_ = width;
}

void Grid::setGcellHeight(int height) {
    gcellHeight_ = height;
}

void Grid::setWireCapacity(int wcap) {
    wireCapacity_ = wcap;
}

void Grid::setTrackSupply(int tsup) {
    trackSupply_ = tsup;
}

void Grid::setNumLayers(int nlyr) {
    numLayers_ = nlyr;
}

Rect Grid::getBoundary() {
    return bbox_;
}


};


    /*
        
        bgBox queryBox = gcell->getQueryBox();
        vector< pair<bgBox, dbInst*> > queryResults;
        instRtree.query(bgi::intersects(queryBox), std::back_inserter(queryResults));
        for(auto& val : queryResults) {
            dbInst* inst = val.sceond;
            gcell->addInst(inst);
        }
    }






            queryResults.clear();
            bgSeg querySeg(bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()));
            gcellRtree->query(bgi::intersects(querySeg), std::back_inserter(queryResults));

            for(pair<bgBox,Gcell*> val : queryResults) {
                Gcell* gcell = val.second;
                Rect rect = gcell->getBBox();
                gcell->updateResourceModelRSMT(seg);
            }
        // search intersecting gcells 
        bgBox qeuryBox = myRSMT->getQueryBBox();
        queryResults.clear();
        gcellRtree->query(bgi::intersects(queryBox), std::back_inserter(queryResults));
        
        // update parital RUDY
        for(pair<bgBox,Gcell*> val : queryResults) {
            Gcell* gcell = val.second;

            odb::Rect gcellBBox = gcell->getRect();
            odb::Rect netBBox = myRSMT->getBBox();
            
            
            uint64 gcellArea = gcellBBox.area();
            uint64 intersectArea = gcellBBox.intersect(netBBox);

            double dn = myRSMT->getWireUniformDensity();
            double R = 1.0 * intersectArea / gcellArea;
            double partial_RUDY = dn*R;

            gcell->addRUDY( parital_RUDY );
        }
    */


