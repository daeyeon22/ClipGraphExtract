#include "opendb/geom.h"







void
ClipGraphExtractor::initGrid(int numRows, int maxRtLayer) {

    grid = new Grid();
    dbBlock* block = getDb()->getChip()->getBlock();
    
    // Calculate grid size = rowHeight * numRows
    dbSite* site = block->getRows().begin()->getSite(); 

    int gcellWidth = site->getHeight() * numRows();
    int gcellHeight = site->getHeight() * numRows();
    
    odb::Rect blockArea;
    block->getBBox()->getBox(blockArea);

    grid->setBoundary(blockArea);
    grid->setGcellWidth(gcellWidth);
    grid->setGcellHeight(gcellHeight);
    grid->initGrid();
    
    //  
    


    // get layer pitches
    unordered_map<dbTechLayer*, uint> layerMinWidth;
    unordered_map<dbTechLayer*, uint> layerMinSpacing;
    dbSet<dbTechLayer> techLayers = getDb()->getTech()-getLayers();
    for(dbTechLayer* layer : techLayers) {
        layerMinWidth[layer] = layer->getWidth();
        layerMinSpacing[layer] = layer->getSpacing();
        //if(layer->getType() == dbTechLayerType::ROUTING) {
        //}
    }

    typedef bg::model::point<int, 2, bg::cs::cartesian> bgPoint;
    typedef bg::model::box<bgPoint> bgBox;
    typedef bg::model::segment<bgPoint> bgSeg;


    bgi::rtree<pair<bgBox, dbInst*>, bgi::quadratic<6>> instRtree;
    bgi::rtree<pair<bgSeg, dbNet*>, bgi::quadratic<6>> egrWireRtree;
    bgi::rtree<pair<bgBox, Gcell*>, bgi::quadratic<6>> gcellRtree;
    bgi::rtree<pair<bgSeg, RtTree*>, bgi::quadratic<6>> rtSegRtree;

    bgi::rtree<pair<bgPoint, dbInst*>, bgi::quadratic<6>> viaRtree;

    
    // make gcellRtree
    for( Gcell* gcell : grid->getGcells() ) {
        odb::Rect rect;
        gcell->getBox(rect);

        bgBox gcellBox(bgPoint(rect.xMin(), rect.yMin()), bgPoint(rect.xMax(), rect.yMax()));
        gcellRtree.insert( make_pair(gcellBox, gcell) );
    }


    // make instRtree
    for( dbInst* inst : block->getInsts() ) {
        dbBox* bBox = inst->getBBox();
        bgBox b (bgPoint(bBox->xMin(), bBox->yMin()),
                    bgPoint(bBox->xMax(), bBox->yMax()));

        instRtree.insert( make_pair(b, inst) );
    };


    // init FLUTE
    Flute::readLUT();

    // wireRtree (eGR result)
    for( dbNet* net : block->getNets()) {
        dbSet<dbITerm> iterms = net->getITerms();
        
        // add terminals
        int x,y;
        RtTree* myRtTree = grid->createRoutingTree(net);
        for(dbITerm* iterm : net->getITerms()) {
            iterm->getAvgXY(&x, &y);
            myRtTree->addTerminal(x,y);
        }
        for(dbBTerm* bterm : net->getBTerms()) {
            bterm->getAvgXY(&x, &y);
            myRtTree->addTerminal(x,y);
        }

        // create RSMT
        myRtTree->createRSMT();


        vector<pair<bgBox, Gcell*>> queryResults;
        vector<odb::Rect> segments = myRtTree->decomposeRSMT();

        for(odb::Rect& seg : segments) {

            // update #cut-nets
            queryResults.clear();
            bgSeg querySeg(bgPoint(seg.xMin(), seg.yMin()), bgPoint(seg.xMax(), seg.yMax()));
            gcellRtree->query(bgi::intersects(querySeg), std::back_inserter(queryResults));

            for(pair<bgBox,Gcell*> val : queryResults) {
                Gcell* gcell = val.second;
                Rect rect = gcell->getBBox();

                int xMin = gcell->getBBox()->xMin();
                int xMax = gcell->getBBox()->xMax();
                int yMin = gcell->getBBox()->yMin();
                int yMax = gcell->getBBox()->yMax();

                if(seg.xMin() == seg.xMax() && seg.yMin() != seg.yMax()) {
                    // vertical
                    if(seg.yMin() < yMin) {

                    } else {

                    }
                } else if (seg.xMin() != seg.xMax() && seg.yMin() == seg.yMax()) {
                    // horizontal
                } else {    
                    // point
                }




                if( seg.xMin() <= xMin ) {

                } else if (seg.xMin() > xMin && seg.xMin() < xMax) {

                } else {

                }







            }



        }

        



        // search intersecting gcells 
        bgBox qeuryBox = myRtTree->getQueryBBox();


        gcellRtree->query(bgi::intersects(queryBox), std::back_inserter(queryResults));
        
        // update parital RUDY
        for(pair<bgBox,Gcell*> val : queryResults) {
            Gcell* gcell = val.second;

            odb::Rect gcellBBox = gcell->getRect();
            odb::Rect netBBox = myRtTree->getBBox();
            
            
            uint64 gcellArea = gcellBBox.area();
            uint64 intersectArea = gcellBBox.intersect(netBBox);

            double dn = myRtTree->getWireUniformDensity();
            double R = 1.0 * intersectArea / gcellArea;
            double partial_RUDY = dn*R;

            gcell->addRUDY( parital_RUDY );
        }

        // update cut-nets (RSMT)
        







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

                            bgBox wireBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                            wireRtree->insert( make_pair( wireBox, net ) );
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


    //
    for( Gcell* gcell : grid->getGcells() ) {
        odb::Rect rect;
        gcell->getBox(rect);

        bgBox queryBox( bgPoint(rect.xMin(), rect.yMin()), bgPoint(rect.xMax(), rect.yMax()) );
        vector< pair<bgBox, dbInst*> > foundInsts;









    }





}


