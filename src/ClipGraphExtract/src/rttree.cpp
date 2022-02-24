#include "grid.h"





odb::Rect RtTree::getBBox() {
    return bbox_;
}


bgBox RtTree::getQueryBox() {
    return bgBox( bgPoint( bbox_.xMin(), bbox_.yMin()), bgPoint( bbox_.xMax(), bbox_.yMax() ) );
}

void RtTree::addTerminal(int x, int y) {
    terminals_.push_back(odb::Point(x,y));
}

bool RtTree::isGlobalNet() {
    return !isLocalNet();
}


bool RtTree::isLocalNet() {
    if(rsmtOverlaps_.size() > 1)
        return false;
    else
        return true;
}


vector<odb::Rect> RtTree::decomposeRSMT() {
    vector<odb::Rect> segments;

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
            if(rand() % 2 == 0) {
                seg1 = odb::Rect( x1, y1, x1, y2 );
                seg2 = odb::Rect( x1, y2, x2, y2 );

            } else {
                seg1 = odb::Rect( x1, y1, x2, y1 );
                seg2 = odb::Rect( x2, y1, x2, y2 );
            }
            segments.push_back(seg1);
            segments.push_back(seg2);
        } else {
            segments.push_back( odb::Rect( x1, y1, x2, y2 ) );
        }
    }

    return segments;
}


void RtTree::createRSMT() {

    int xMin = INT_MAX;
    int yMin = INT_MAX;
    int xMax = INT_MIN;
    int yMax = INT_MIN;
    
    int deg = terminals_.size();
    int xs[n] = {0};
    int ys[n] = {0};

    for(int i=0; i < n; i++) {
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
    bbox_ = odb::Rect(xMin, yMin, xMax, yMax);
}


int 
RtTree::getWireLengthRSMT() {
    return Flute::wirelength(rsmt_);   
}


int
RtTree::getWireLengthHPWL() {
    return bbox_.dx() + bbox_.dy(); 
}



