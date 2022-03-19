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
using namespace ClipGraphExtract;
using namespace ClipGraphExtract;


void
ClipGraphExtractor::saveGridImages( const char* imgDir ) {
    Grid* grid = (Grid*)(grid_);
    grid->saveGridImages(string(imgDir));
}


namespace ClipGraphExtract {

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


