#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"

#include "sta/Graph.hh"
#include "sta/Sta.hh"
#include "sta/Network.hh"
#include "sta/Liberty.hh"
#include "sta/Sdc.hh"
#include "sta/PortDirection.hh"
#include "sta/Corner.hh"
#include "sta/PathExpanded.hh"
#include "sta/PathEnd.hh"
#include "sta/PathRef.hh"
#include "sta/Search.hh"
#include "sta/Bfs.hh"

#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"



#include "CImg.h"
#include "instGraph.h"
#include "binGraph.h"
#include <iostream>
#include <string>
#include <sstream>
#include <regex>
#include <fstream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bgi::rtree< std::pair<box, inst_value>, 
                    bgi::quadratic<6> > inst_RTree;

typedef bgi::rtree< std::pair<box, wire_value>, 
                    bgi::quadratic<6> > wire_RTree;

typedef bgi::rtree< std::pair<box, via_value>, 
                    bgi::quadratic<6> > via_RTree;

typedef bgi::rtree< std::pair<box, pin_value>, 
                    bgi::quadratic<6> > pin_RTree;

typedef bgi::rtree< std::pair<box, drc_value>,
					bgi::quadratic<6> > drc_RTree;

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
    inst_rTree_(nullptr), wire_rTree_(nullptr), via_rTree_(nullptr), pin_rTree_(nullptr),
	drc_rTree_(nullptr), graphModel_(Star), edgeWeightModel_(A), fileName_("") {};

ClipGraphExtractor::~ClipGraphExtractor() {
  clear(); 
}

void
ClipGraphExtractor::clear() {
  db_ = nullptr; 
  sta_ = nullptr;
  if( inst_rTree_ ) {
    delete (inst_RTree*) inst_rTree_;
  }

  if( wire_rTree_ ) {
    delete (wire_RTree*) wire_rTree_;
  }

  if( via_rTree_ ) {
    delete (via_RTree*) via_rTree_;
  }

  if( pin_rTree_ ) {
    delete (pin_RTree*) pin_rTree_;
  }

  if( drc_rTree_ ) {
    delete (pin_RTree*) drc_rTree_;
  }

  if( binGraph_ ) {
    delete (bingraph::Graph*) binGraph_;
  }

  binGraph_ = nullptr;
  rTree_ = nullptr;
  inst_rTree_ = nullptr;
  wire_rTree_ = nullptr;
  via_rTree_ = nullptr;
  pin_rTree_ = nullptr;
  drc_rTree_ = nullptr;
  graphModel_ = Star;
  edgeWeightModel_ = A;
  fileName_ = "";
}

void 
ClipGraphExtractor::init() {  
    using namespace odb;
    inst_rTree_ = (void*) (new inst_RTree);
    wire_rTree_ = (void*) (new wire_RTree);
    via_rTree_ = (void*) (new via_RTree);
    pin_rTree_ = (void*) (new pin_RTree);
    drc_rTree_ = (void*) (new drc_RTree);
    
    binGraph_ = (void*) (new bingraph::Graph); 
    
	inst_RTree* inst_rTree = (inst_RTree*) inst_rTree_;
    wire_RTree* wire_rTree = (wire_RTree*) wire_rTree_;
    via_RTree* via_rTree = (via_RTree*) via_rTree_;
    pin_RTree* pin_rTree = (pin_RTree*) pin_rTree_;
    drc_RTree* drc_rTree = (drc_RTree*) drc_rTree_;

    unordered_map<unsigned int, int> layerPitches;
    dbSet<dbTechLayer> layers = db_->getTech()->getLayers();
    for(dbTechLayer* layer : layers){
        if(layer->getName()[0] == 'M' || layer->getName()[0] == 'm'){
            unsigned int layerNum = (unsigned int) layer->getName()[1] - '0';
            unsigned int layerWidth = layer->getWidth();
            unsigned int layerSpacing = layer->getSpacing();

            layerPitches[layerNum] = layerWidth + layerSpacing;
        }
    }

    // DB Query to fill in inst_RTree
	dbBlock* block = db_->getChip()->getBlock();
    for( dbInst* inst : block->getInsts() ) {
        dbBox* bBox = inst->getBBox();
        box b (point(bBox->xMin(), bBox->yMin()), 
            point(bBox->xMax(), bBox->yMax()));
        inst_value temp{inst, b};
        inst_rTree->insert( make_pair(b, temp) );
    }
    
    // DB Query to fill in pin_RTree
    dbSet<dbITerm> iterms = block->getITerms();
    for(dbITerm* iterm : iterms){
        int x, y;
        pin_value p;

        iterm->getAvgXY(&x, &y);
        p.name_ = iterm->getMTerm()->getMaster()->getName()+"/"+iterm->getMTerm()->getName();
        p.setBox(x, y);
        pin_rTree->insert( make_pair(p.box_, p) );
    }

    // DB Query to fill in wire_RTree and via_RTree
    dbSet<dbNet> nets = block->getNets();
    for ( dbNet* net : nets ){ // net
        dbWire* wire = net->getWire();

        if ( wire && wire->length() ){ // wire
            dbWireDecoder decoder;
            decoder.begin(wire);
            
            vector< pair<int, int> > point;
            wire_value w;
            w.pitches_ = layerPitches;
            via_value v;
            while(decoder.peek() != dbWireDecoder::END_DECODE){
                unsigned int  layerNum;
                unsigned int layerWidth;
                unsigned int layerSpacing;

                switch(decoder.peek()){
                    case dbWireDecoder::PATH:
                    case dbWireDecoder::JUNCTION:
                    case dbWireDecoder::SHORT:{
                        decoder.next();
                        dbTechLayer* layer = decoder.getLayer();
                        layerNum = (unsigned int) layer->getName()[1] - '0';
                        layerWidth = layer->getWidth();
                        layerSpacing = layer->getSpacing();

                        point.clear();

                        w.layerNum_ = layerNum;
                        w.net_ = net;
                        w.width_ = layerWidth;
                        w.spacing_ = layerSpacing;

                        break;
                                                 }
                    case dbWireDecoder::VIA:
                    case dbWireDecoder::TECH_VIA:{
                        decoder.next();
                        dbTechLayer* layer = decoder.getLayer();
                        layerNum = (unsigned int) layer->getName()[1] - '0';
                        layerWidth = layer->getWidth();
                        
                        pair<int, int> viaPoint = point.back();
                        
                        v.bottomLayer_ = layerNum;
                        v.net_ = net;
                        v.width_ = layerWidth;
                        v.setBox(viaPoint.first, viaPoint.second);
                        via_rTree->insert( make_pair(v.box_, v) );
                        break;
                                                 }
                    case dbWireDecoder::POINT:{
                        int x, y;
                        decoder.next();
                        decoder.getPoint(x, y);
                        point.push_back( make_pair(x, y) );
                        if(point.size() == 2){
                            int fx = point[0].first;
                            int fy = point[0].second;
                            int tx = point[1].first;
                            int ty = point[1].second;
                            
                            w.setBox(fx, fy, tx, ty);
							if(tx-fx == 0) 
								w.type_ = 'V';
							else if(ty-fy == 0) 
								w.type_ = 'H';
							wire_rTree->insert( make_pair(w.box_, w) );
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
}

void
ClipGraphExtractor::extract(int lx, int ly, int ux, int uy) {
  inst_RTree* inst_rTree = (inst_RTree*) inst_rTree_;
  sta_->updateTiming(false);
  sta::dbNetwork* network = sta_->getDbNetwork();
  sta::Graph* graph = sta_->ensureGraph();
  
  box queryBox( point(lx, ly), point(ux, uy) );

  vector< pair<box, inst_value> > foundInsts; 
  inst_rTree->query(bgi::intersects(queryBox), 
      std::back_inserter(foundInsts));

  cout << "NumFoundInsts: " << foundInsts.size() << endl;

  set<odb::dbInst*> instSet;
  for(pair<box, inst_value>& val : foundInsts) {
    odb::dbInst* inst = val.second.inst_;
    instSet.insert( inst ); 
  }
  
  Graph instGraph;
  instGraph.setDb(db_);
  instGraph.init(instSet, graphModel_, edgeWeightModel_);
  instGraph.saveFile(fileName_);
  cout << "Done!" << endl;
}

// added by dykim
void

ClipGraphExtractor::extractBinGraph(int numRows, int maxLayer) {

    cout << "Extract bin graph" << endl;
    inst_RTree* inst_rTree = (inst_RTree*) inst_rTree_;
    wire_RTree* wire_rTree = (wire_RTree*) wire_rTree_;
    via_RTree* via_rTree = (via_RTree*) via_rTree_;
    pin_RTree* pin_rTree = (pin_RTree*) pin_rTree_;
    
    bingraph::Graph* binGraph = (bingraph::Graph*) binGraph_;
    dbTech* tech = db_->getTech();
    dbChip* chip = db_->getChip();
    dbBlock* block = chip->getBlock();

    int xMin = block->getBBox()->xMin();
    int xMax = block->getBBox()->xMax();
    int yMin = block->getBBox()->yMin();
    int yMax = block->getBBox()->yMax();


    dbSite* site = block->getRows().begin()->getSite();
    uint blockWidth = block->getBBox()->getDX();
    uint blockHeight = block->getBBox()->getDY();
    uint siteHeight = site->getHeight();
    uint siteWidth = site->getWidth();


    int unitSize = numRows * siteHeight;

    cout << xMin << " " << yMin << " " << xMax << " " << yMax << endl;
    cout << "unit size : " << unitSize << endl;

    
    int numVertices = std::ceil( 1.0* blockWidth / unitSize ) * std::ceil( 1.0*blockHeight / unitSize );
   
    int blockNum = 0;

    for(int lx=xMin; lx <= xMax-unitSize; lx+=unitSize) {
        for(int ly=yMin; ly <= yMax-unitSize; ly+=unitSize) {
            
            blockNum++;

            int ux = lx + unitSize;
            int uy = ly + unitSize;
            //cout << "Block Num: " << blockNum << "\t";
            //cout << lx/2000.0 << " " << ly/2000.0 << " " 
            //     << ux/2000.0 << " " << uy/2000.0 << endl;

            box queryBox( point(lx, ly), point(ux, uy) );
            vector< pair<box, inst_value> > foundInsts;
            vector< pair<box, wire_value> > foundWires;
            vector< pair<box, via_value> > foundVias;
            vector< pair<box, pin_value> > foundPins;

            inst_rTree->query(bgi::intersects(queryBox), 
                    std::back_inserter(foundInsts));
            
            wire_rTree->query(bgi::intersects(queryBox), 
                    std::back_inserter(foundWires));
            
            via_rTree->query(bgi::intersects(queryBox), 
                    std::back_inserter(foundVias));

            pin_rTree->query(bgi::intersects(queryBox), 
                    std::back_inserter(foundPins));

            vector<odb::dbInst*> insts;
            vector<wire_value> wireValues;
            vector<via_value> viaValues;
            vector<pin_value> pinValues;

            for(pair<box, inst_value>& val : foundInsts) {
                insts.push_back(val.second.inst_);
            }
            
            for(pair<box, wire_value>& val : foundWires) {
                wireValues.push_back(val.second);
                //cout << "wire: " << val.second.layerNum_ << "\t" 
                //                 << val.second.lx_/2000.0 << " " 
                //                 << val.second.ly_/2000.0 << " " 
                //                 << val.second.ux_/2000.0 << " " 
                //                 << val.second.uy_/2000.0 << endl;

            }
            
            for(pair<box, via_value>& val : foundVias) {
                viaValues.push_back(val.second);
                //cout << "via: " << val.second.lx_/2000.0 << " " 
                //                << val.second.ly_/2000.0 << " " 
                //                << val.second.ux_/2000.0 << " " 
                //                << val.second.uy_/2000.0 << endl;

            }

            for(pair<box, pin_value>& val : foundPins){
                pinValues.push_back(val.second);
                //cout << "pin: " << val.second.name_ << " " 
                //                << val.second.x_/2000.0 << " " 
                //                << val.second.y_/2000.0 << endl;
            
            }

            binGraph->addVertex(lx, ly, ux, uy, maxLayer, insts, 
                    wireValues, viaValues, pinValues);        
        }
    }

    //binGraph->initEdges();
    //binGraph->saveFile(fileName_.c_str());
    // temporal block
    cout << "Done" << endl;
}

vector<string> splitAsTokens(string str, string delim){
    vector<string> _tokens;
    size_t start, end=0;
    while(end < str.size()){
        start = end;
        while(start < str.size() && (delim.find(str[start]) != string::npos)){
            start++;
        }
        end = start;
        while(end < str.size() && (delim.find(str[end]) == string::npos)){
            end++;
        }
        if(end-start != 0){
            _tokens.push_back(string(str, start, end-start));
        }
    }
    return _tokens;
}

void
ClipGraphExtractor::labelingBinGraph(const char* invRoutingReport) {
    drc_RTree* drc_rTree = (drc_RTree*) drc_rTree_;

    ifstream inFile(invRoutingReport);
 
    const std::regex r1("Bounds[ \t\r\n\v\f]:[ \t\r\n\v\f]");
    string line;

	dbBlock* block = db_->getChip()->getBlock();
    int dbUnitMicron = block->getDbUnitsPerMicron();

    while(getline(inFile, line)) {
        std::smatch match;
        if(regex_search(line, match, r1)) {
            string tail = match.suffix();
            
            string delim = " (),";
            vector<string> tokens = splitAsTokens(tail, delim);
            if(tokens.size() != 4) {
            
            }

            int lx = dbUnitMicron * atof(tokens[0].c_str());
            int ly = dbUnitMicron * atof(tokens[1].c_str());
            int ux = dbUnitMicron * atof(tokens[2].c_str());
            int uy = dbUnitMicron * atof(tokens[3].c_str());
    
            //cout << "(" << lx << " " << ly << ") (" << ux << " " << uy << ")" << endl;
            box b (point(lx, ly), point(ux, uy));
			drc_value d;
		
			d.setBox((int)lx, (int)ly, (int)ux, (int)uy);
            drc_rTree->insert( make_pair(b, d) );
        }
    }

    bingraph::Graph* binGraph = (bingraph::Graph*)binGraph_;
    for(bingraph::Vertex* vert : binGraph->getVertices()) {
        int lx = vert->getLx();
        int ux = vert->getUx();
        int ly = vert->getLy();
        int uy = vert->getUy();

        box queryBox( point(lx, ly), point(ux, uy) );
        vector< pair<box, drc_value> > foundDrcs;
        
		drc_rTree->query(bgi::intersects(queryBox), 
				std::back_inserter(foundDrcs));

		for(pair<box, drc_value>& val : foundDrcs) {
			vert->addDrcValue(val.second);
		}

        int label = foundDrcs.size() > 0 ? 1 : 0;
        vert->setLabel(label);
    }
}


void
ClipGraphExtractor::saveBinGraph() {
    bingraph::Graph* binGraph = (bingraph::Graph*) binGraph_;   
    binGraph->saveFile(prefix_.c_str());
}

void
ClipGraphExtractor::showCongestionMap() {
    bingraph::Graph* binGraph = (bingraph::Graph*) binGraph_;   
    binGraph->showCongestion();

    


}



void
ClipGraphExtractor::setDb(odb::dbDatabase* db) {
  db_ = db;
}

void
ClipGraphExtractor::setSta(sta::dbSta* sta) {
  sta_ = sta;
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
