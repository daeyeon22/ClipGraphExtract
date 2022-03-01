#include "grid.h"
#include <set>
#include <vector>

namespace feature_extractor {

using namespace std;
using namespace odb;

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
}

void RSMT::setMinWidth(int width) {
    minWidth_ = width;
}

vector<Rect> RSMT::getSegments() {
    vector<Rect> segments;
    int n, x1, y1, x2, y2;
    Flute::Tree *tree = &rsmt_;
    for(int i=0; i < tree->deg -2; i++) {
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

    // RSMT
    rsmt_ = Flute::flute(deg, xs, ys, FLUTE_ACCURACY);
    // BBOX
    bbox_ = Rect(xMin, yMin, xMax, yMax);
}


int 
RSMT::getWireLengthRSMT() {
    return Flute::wirelength(rsmt_);   
}


int
RSMT::getWireLengthHPWL() {
    return bbox_.dx() + bbox_.dy(); 
}

double
RSMT::getWireUniformDensity() {
    int totalArea = bbox_.dx() * bbox_.dy();
    int wireArea = minWidth_ * getWireLengthRSMT();
    return 1.0* wireArea / totalArea;
}


void
RSMT::searchOverlaps(BoxRtree<Gcell*> &tree) {
    
    //
    vector<pair<bgBox, Gcell*>> queryResults;
    set<Gcell*> overlaps;

    // RSMT
    for(Rect& seg : getSegments()) {
        queryResults.clear();
        bgSeg bg_seg ( bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()) );
        
        tree.query(bgi::intersects(bg_seg), back_inserter(queryResults));
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
    tree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));

    for(auto& val : queryResults) {
        Gcell* gcell = val.second;
        overlaps.insert(gcell);
    }


    bboxOverlaps_ = vector<Gcell*>(overlaps.size());
    copy(overlaps.begin(), overlaps.end(), bboxOverlaps_.begin());
}







};

