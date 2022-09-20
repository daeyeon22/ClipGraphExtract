#include "grid.h"
#include <set>
#include <vector>
#include <assert.h>
#include <iostream>

namespace ClipGraphExtract {

using namespace std;
using namespace odb;

RSMT::RSMT() : net_(nullptr), width_(0), wl_(0), hpwl_(0) {}

RSMT::RSMT(dbNet* net) :
    net_(net), width_(0), wl_(0), hpwl_(0) {}

Rect RSMT::getBBox() {
    return bbox_;
}

bgBox RSMT::getQueryBox() {
    return bgBox( bgPoint( bbox_.xMin(), bbox_.yMin()), bgPoint( bbox_.xMax(), bbox_.yMax() ) );
}

void RSMT::addTerminal(int x, int y) {
    terminals_.push_back(odb::Point(x,y));
}

bool RSMT::isGlobalNet() {
    return !isLocalNet();
}

bool RSMT::isLocalNet() {
	return bboxOverlaps_.size() < 2 ? true : false;
    //return true;
}

void RSMT::setWireWidth(int width) {
    width_ = width;
}

int RSMT::getNumTerms() {
    return (int)terminals_.size();
}


vector<Rect> RSMT::getSegments() {

    
    // NEED TO DEBUG
    //
    vector<Rect> segments;
    
    if(terminals_.size() > 1) {
        int n, x1, y1, x2, y2;
        Flute::Tree *tree = &rsmt_;
        for(int i=0; i < 2*tree->deg -2; i++) {
            n = tree->branch[i].n;
            x1 = tree->branch[i].x;
            y1 = tree->branch[i].y;
            x2 = tree->branch[n].x;
            y2 = tree->branch[n].y;

            // if barnch has a L shape, decomposes into 2 semgnets
            if(x1 != x2 && y1 != y2) {
                Rect seg1, seg2;
                if(rand() % 2 == 0) {
                    seg1 = Rect( x1, y1, x1, y2 );
                    seg2 = Rect( x1, y2, x2, y2 );

                } else {
                    seg1 = Rect( x1, y1, x2, y1 );
                    seg2 = Rect( x2, y1, x2, y2 );
                }
                segments.push_back(seg1);
                segments.push_back(seg2);
            } else {
                segments.push_back( Rect( x1, y1, x2, y2 ) );
            }
        }
    }

    return segments;
}


void RSMT::createTree() {

    int xMin = INT_MAX;
    int yMin = INT_MAX;
    int xMax = INT_MIN;
    int yMax = INT_MIN;
    
    int deg = terminals_.size();
    int xs[deg] = {0};
    int ys[deg] = {0};

    for(int i=0; i < deg; i++) {
        xs[i] = terminals_[i].getX();
        ys[i] = terminals_[i].getY();

        // get BBox
        xMin = min(xMin, xs[i]);
        yMin = min(yMin, ys[i]);
        xMax = max(xMax, xs[i]);
        yMax = max(yMax, ys[i]);

    }

    if(deg > 1) {
        // RSMT
        rsmt_ = Flute::flute(deg, xs, ys, FLUTE_ACCURACY);
        wl_ = Flute::wirelength(rsmt_); 
    }
    // BBOX
    bbox_ = Rect(xMin, yMin, xMax, yMax);
}


int 
RSMT::getWireLengthRSMT() {
    return wl_;
    //if(terminals_.size() > 1) 
    //    return Flute::wirelength(rsmt_);   
    //else
    //    return 0;
}


int
RSMT::getWireLengthHPWL() {
    return bbox_.dx() + bbox_.dy(); 
}

double
RSMT::getWireUniformUtil() {
    uint totalArea = bbox_.dx() * bbox_.dy();
    uint wireArea = width_ * getWireLengthRSMT();
    double u_den = (totalArea == 0)? 1.0 : min(1.0, 1.0*wireArea/totalArea);
    //double u_den = 1.0 * wireArea / totalArea;
    assert(u_den < 0);
    return u_den;
}


void
RSMT::searchOverlaps(BoxRtree<Gcell*> *tree) {
    
    //
    vector<pair<bgBox, Gcell*>> queryResults;
    set<Gcell*> overlaps;

    // RSMT
    for(Rect& seg : getSegments()) {
        queryResults.clear();
        bgSeg bg_seg ( bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()) );
        
        tree->query(bgi::intersects(bg_seg), back_inserter(queryResults));
        for(auto& val : queryResults) {
            Gcell* gcell = val.second;
            overlaps.insert(gcell);
        }
    }

    rsmtOverlaps_ = vector<Gcell*>(overlaps.size());
    copy(overlaps.begin(), overlaps.end(), rsmtOverlaps_.begin());

    // BBox
    overlaps.clear();
    queryResults.clear();
    tree->query(bgi::intersects(getQueryBox()), back_inserter(queryResults));
	
    for(auto& val : queryResults) {
        Gcell* gcell = val.second;
        overlaps.insert(gcell);
    }

	bboxOverlaps_ = vector<Gcell*>(overlaps.size());
	copy(overlaps.begin(), overlaps.end(), bboxOverlaps_.begin());
}







};

