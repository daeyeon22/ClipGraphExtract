#include "opendb/geom.h"
#include "grid.h"
#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include <iostream>
#include <string>
#include <vector>
#include <cassert>


using namespace std;
using namespace odb;
using namespace feature_extractor;
using namespace ClipGraphExtract;

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
    for(dbTechLayer* layer : techLayers) {
        if(layer->getType() == dbTechLayerType::ROUTING) {
            numLayer++;
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
    // Get core area
    odb::Rect blockArea;
    block->getBBox()->getBox(blockArea);

    // Initialize Gcell Grid
    grid->setBoundary(blockArea);
    grid->setGcellWidth(gcellWidth);
    grid->setGcellHeight(gcellHeight);
    grid->setWireCapacity(wireCapacity);
    grid->setTrackSupply(trackSupply);
    grid->init();
    //  

    // init rtree
    BoxRtree<dbInst*> instRtree;
    BoxRtree<Gcell*> gcellRtree;
    SegRtree<dbNet*> egrRtree;
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
    for( dbNet* net : block->getNets()) {
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
            dbWireDecoder decoder;
            decoder.begin(wire);
            
            vector<odb::Point> points;

            while(decoder.peek() != dbWireDecoder::END_DECODE) {
                uint layerNum;
                uint layerMinWidth;
                uint layerMinSpacing;

                switch(decoder.peek()) {
                    case dbWireDecoder::SHORT: {
                        decoder.next();
                        points.clear();            
                        break;
                    }
                    case dbWireDecoder::TECH_VIA: {
                        decoder.next();                   
                        break;
                    }
                    case dbWireDecoder::POINT: {
                        decoder.next();                   
                        int x,y;
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        if(points.size() >1) {
                            odb::Point pt1 = points[points.size()-2];
                            odb::Point pt2 = points[points.size()-1];
                            
                            int xMin = min(pt1.getX(), pt2.getX());// - layerMinWidth[decoder.getLayer()]/2;
                            int xMax = max(pt1.getX(), pt2.getX());// + layerMinWidth[decoder.getLayer()]/2;
                            int yMin = min(pt1.getY(), pt2.getY());// - layerMinWidth[decoder.getLayer()]/2;
                            int yMax = max(pt1.getY(), pt2.getY());// + layerMinWidth[decoder.getLayer()]/2;

                            bgSeg wireSeg( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                            egrRtree.insert( make_pair( wireSeg, net ) );
                        }
                        break;
                    }
                    default:
                        decoder.next();
                        break;
                }

            }
        }
        // search overlapping gcells
        myRSMT->searchOverlaps(gcellRtree);
    }
    // add Instance in gcell 
    for( Gcell* gcell : grid->getGcells() ) {
        //
        gcell->extractFeatureEGR(egrRtree);
        gcell->extractFeaturePL(instRtree);
        gcell->extractFeatureRSMT(rsmtRtree);
    }
}


namespace feature_extractor {

void Grid::init() {
    assert(gcellWidth_ == 0);
    assert(gcellHeight_ == 0);
    assert(bbox_.dx() * bbox_.dy() == 0);

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
        }
    }
}

vector<Gcell*> Grid::getGcells() {
    return gcells_;
}

Gcell* Grid::createGcell(int x1, int y1, int x2, int y2) {
    Gcell* gcell = new Gcell;
    gcell->setBoundary(Rect(x1,y1, x2,y2));
    gcells_.push_back(gcell);
    return gcell;
}


RSMT* Grid::createRSMT(odb::dbNet* net) {
    RSMT* myRSMT = new RSMT(net);
    dbSet<dbITerm> iterms = net->getITerms();

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

    // create RSMT
    myRSMT->createTree();
    rsmts_.push_back(myRSMT);
    return myRSMT;
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


