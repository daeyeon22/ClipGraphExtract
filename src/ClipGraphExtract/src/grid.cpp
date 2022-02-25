#include "opendb/geom.h"






typedef bg::model::point<int, 2, bg::cs::cartesian> bgPoint;
typedef bg::model::box<bgPoint> bgBox;
typedef bg::model::segment<bgPoint> bgSeg;

void ClipGraphExtractor::initGrid(int numRows, int maxLayer) {

    grid = new Grid();
    dbBlock* block = getDb()->getChip()->getBlock();
    
    // Calculate grid size = rowHeight * numRows
    dbSite* site = block->getRows().begin()->getSite(); 

    int gcellWidth = site->getHeight() * numRows();
    int gcellHeight = site->getHeight() * numRows();

    // get routing supply / capacity for each gcell
    dbSet<dbTechLayer> techLayers = getDb()->getTech()-getLayers();
    int numLayer=0;
    int numSupply=0;
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
            
            numSupply+= supply;
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
    grid->initGrid();
    //  
   
    template<typename A> 
    using bgi::rtree<pair<bgBox, A>, bgi::quadratic<6>> BoxRtree;
    template<typename A> 
    using bgi::rtree<pair<bgSeg, A>, bgi::quadratic<6>> SegRtree;

    BoxRtree<dbInst*> instRtree;
    BoxRtree<Gcell*> gcellRtree;
    SegRtree<dbNet*> egrRtree;
    SegRtree<RSMT*> rsmtRtree;

    //bgi::rtree<pair<bgBox, dbInst*>, bgi::quadratic<6>> instRtree;
    //bgi::rtree<pair<bgBox, Gcell*>, bgi::quadratic<6>> gcellRtree;
    //bgi::rtree<pair<bgSeg, dbNet*>, bgi::quadratic<6>> egrWireRtree;
    //bgi::rtree<pair<bgSeg, RSMT*>, bgi::quadratic<6>> rtSegRtree;
    //bgi::rtree<pair<bgPoint, dbInst*>, bgi::quadratic<6>> viaRtree;
    
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
        dbSet<dbITerm> iterms = net->getITerms();
        
        // add terminals
        int x,y;
        RSMT* myRSMT = grid->createRSMT(net);
        for(dbITerm* iterm : net->getITerms()) {
            iterm->getAvgXY(&x, &y);
            myRSMT->addTerminal(x,y);
        }
        for(dbBTerm* bterm : net->getBTerms()) {
            bterm->getAvgXY(&x, &y);
            myRSMT->addTerminal(x,y);
        }

        // create RSMT
        myRSMT->createRSMT();
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
            dcoder.begin(wire);
            
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
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        if(point.size() >1) {
                            odb::Point pt1 = points[points.size()-2];
                            odb::Point pt2 = points[points.size()-1];
                            
                            int xMin = min(pt1.getX(), pt2.getX()) - layerMinWidth[decoder.getLayer()]/2;
                            int xMax = max(pt1.getX(), pt2.getX()) + layerMinWidth[decoder.getLayer()]/2;
                            int yMin = min(pt1.getY(), pt2.getY()) - layerMinWidth[decoder.getLayer()]/2;
                            int yMax = max(pt1.getY(), pt2.getY()) + layerMinWidth[decoder.getLayer()]/2;

                            bgSeg wireSeg( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                            egrWireRtree->insert( make_pair( wireSeg, net ) );
                        }
                        
                        break;
                    }

                    default:
                        decoder.next();
                        break;

                }

            }
        }
    }




    // add Instance in gcell 
    for( Gcell* gcell : grid->getGcells() ) {
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




}


