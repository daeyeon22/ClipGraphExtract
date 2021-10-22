#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"

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

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "instGraph.h"
#include "binGraph.h"
#include <iostream>
#include <string>
#include <sstream>
#include <regex>
#include <fstream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// for easier coding with boost
typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;

// save point and dbInst* pointer 
typedef std::pair<box, odb::dbInst*> value;
typedef bgi::rtree< value, bgi::quadratic<6> > RTree;

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

ClipGraphExtractor::ClipGraphExtractor() : 
  db_(nullptr), sta_(nullptr), rTree_(nullptr),
  graphModel_(Star), edgeWeightModel_(A), 
  fileName_("") {};

ClipGraphExtractor::~ClipGraphExtractor() {
  clear(); 
}

void
ClipGraphExtractor::clear() {
  db_ = nullptr; 
  sta_ = nullptr;
  if( rTree_ ) {
    delete (RTree*) rTree_;
  }

  if( binGraph_ ) {
    delete (bingraph::Graph*) binGraph_;
  }

  binGraph_ = nullptr;
  rTree_ = nullptr;
  graphModel_ = Star;
  edgeWeightModel_ = A;
  fileName_ = "";
}

void 
ClipGraphExtractor::init() {  
    using namespace odb;
    rTree_ = (void*) (new RTree);
    binGraph_ = (void*) (new bingraph::Graph); 
    RTree* rTree = (RTree*)rTree_;

  // DB Query to fill in RTree
	dbBlock* block = db_->getChip()->getBlock();
  for( dbInst* inst : block->getInsts() ) {
    dbBox* bBox = inst->getBBox();
    box b (point(bBox->xMin(), bBox->yMin()), 
        point(bBox->xMax(), bBox->yMax()));
    rTree->insert( make_pair(b, inst) );
  }

}

void
ClipGraphExtractor::extract(int lx, int ly, int ux, int uy) {
  RTree* rTree = (RTree*)rTree_;
  sta_->updateTiming(false);
  sta::dbNetwork* network = sta_->getDbNetwork();
  sta::Graph* graph = sta_->ensureGraph();
  
  box queryBox( point(lx, ly), point(ux, uy) );

  vector<value> foundInsts; 
  rTree->query(bgi::intersects(queryBox), 
      std::back_inserter(foundInsts));

  cout << "NumFoundInsts: " << foundInsts.size() << endl;

  set<odb::dbInst*> instSet;
  for(value& val : foundInsts) {
    odb::dbInst* inst = val.second;
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
ClipGraphExtractor::extractBinGraph(int numRows) {

    cout << "Extract bin graph" << endl;
    //bingraph::Graph binGraph;
    RTree* rTree = (RTree*)rTree_;
    bingraph::Graph* binGraph = (bingraph::Graph*) binGraph_;
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


    
    for(int lx=xMin; lx <= xMax-unitSize; lx+=unitSize) {
        for(int ly=yMin; ly <= yMax-unitSize; ly+=unitSize) {

            int ux = lx + unitSize;
            int uy = ly + unitSize;

            //cout << "(" << lx << " " << ly << ")" << endl;

            box queryBox( point(lx, ly), point(ux, uy) );
            vector<value> foundInsts;

            rTree->query(bgi::intersects(queryBox), std::back_inserter(foundInsts));

            vector<odb::dbInst*> insts;

            for(value& val : foundInsts) {
                insts.push_back(val.second);
            }

            binGraph->addVertex(lx, ly, ux, uy, insts);        
        }
    }

    binGraph->initEdges();
    //binGraph->saveFile(fileName_.c_str());

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


    typedef std::pair<box, int> marker;
    bgi::rtree<marker, bgi::quadratic<6>> drcRtree;
    //RTree* drcMarkerRtree = new RTree;

    ifstream inFile(invRoutingReport);
 
    const std::regex r1("Bounds[ \t\r\n\v\f]:[ \t\r\n\v\f]");
    string line;

	dbBlock* block = db_->getChip()->getBlock();
    int dbUnitMicron = block->getDbUnitsPerMicron();


    while(getline(inFile, line)) {
        std::smatch match;
        if(regex_search(line, match, r1)) {
            //cout << line << endl;
            string tail = match.suffix();
            //cout << tail << endl;
            
            string delim = " (),";
            vector<string> tokens = splitAsTokens(tail, delim);
            //for(auto& token : tokens)
            //    cout << token << endl;
            if(tokens.size() != 4) {
            
            }

            int lx = dbUnitMicron * atof(tokens[0].c_str());
            int ly = dbUnitMicron * atof(tokens[1].c_str());
            int ux = dbUnitMicron * atof(tokens[2].c_str());
            int uy = dbUnitMicron * atof(tokens[3].c_str());
    
            cout << "(" << lx << " " << ly << ") (" << ux << " " << uy << ")" << endl;
            box b (point(lx, ly), point(ux, uy));
            drcRtree.insert( make_pair(b, 1) );
        }
    }

    bingraph::Graph* binGraph = (bingraph::Graph*)binGraph_;
    for(bingraph::Vertex* vert : binGraph->getVertices()) {
        int lx = vert->getLx();
        int ux = vert->getUx();
        int ly = vert->getLy();
        int uy = vert->getUy();

        vector<marker> foundMarkers;
        box queryBox( point(lx, ly), point(ux, uy) );
        drcRtree.query(bgi::intersects(queryBox), std::back_inserter(foundMarkers));


        int label = foundMarkers.size() > 0 ? 1 : 0;
        vert->setLabel(label);
    }
}


void
ClipGraphExtractor::saveBinGraph() {
    bingraph::Graph* binGraph = (bingraph::Graph*) binGraph_;   
    binGraph->saveFile(prefix_.c_str());
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
