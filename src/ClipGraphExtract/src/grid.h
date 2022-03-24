#ifndef __GCELLGRID__
#define __GCELLGRID__
#include <vector>
#include <string>
#include <boost/geometry.hpp>
#include <unordered_map>
#include "flute.h"
#include "opendb/db.h"
#include "opendb/geom.h"
#include "instGraph.h"


// Typedef for boost geometry
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<int, 2, bg::cs::cartesian> bgPoint;
typedef bg::model::box<bgPoint> bgBox;
typedef bg::model::segment<bgPoint> bgSeg;

template<typename A>
using BoxRtree = bgi::rtree<std::pair<bgBox, A>, bgi::quadratic<6>>;
template<typename A>
using SegRtree = bgi::rtree<std::pair<bgSeg, A>, bgi::quadratic<6>>;

namespace ClipGraphExtract {

class Marker;
class Gcell;
class RSMT;


enum Orient {
    LEFT=0,
    RIGHT=1,
    TOP=2,
    BOTTOM=3
};

enum ModelType {
    TREE, ROUTE, DR,EGR,PL
};

class Resource {
  private:
    int trackSupply_[4] = {0};
    int trackDemand_[4] = {0};
    int wireCapacity_;
    int wireLength_;

  public:
    
    Resource() : wireCapacity_(0), wireLength_(0) {}

    double getChanUtil(Orient type) { return 1.0 * trackDemand_[type] / trackSupply_[type]; }
    double getWireUtil() { return 1.0 * wireLength_ / wireCapacity_; }
    int getTrackSupply(Orient type) { return trackSupply_[type]; }
    int getTrackDemand(Orient type) { return trackDemand_[type]; }
    int getWireCapacity() { return wireCapacity_; }
    int getWireLength() { return wireLength_; }
    void incrTrackDemand(Orient type) { trackDemand_[type]++; }
    void addTrackDemand(Orient type, int dem) { trackDemand_[type] += dem; }
    void setTrackSupply(int tSup) {  for(int i=0; i < 4; i++) trackSupply_[i] = tSup; }
    void setWireCapacity(int wCap) {  wireCapacity_ = wCap;  }
    void setWireLength(int wl) { wireLength_ = wl; }
    void addWireLength(int wl) { wireLength_ += wl; }
};


class Gcell {
  private:   

    int col_, row_;

    // Gcell features
    odb::Rect bbox_;

    //Resource rmEGR_; // using GR results
    Resource rmDR_; // using DR results
    Resource rmPL_; // using PLACE results

    // using placement, RSMT results
    int numInsts_;
    int numTerms_;
    //int numLocalNets_;
    //int numGlobalNets_;

    int totalCellArea_;
    int totalPinArea_;
    int numLayers_;

    double cellUtil_;
    double pinUtil_;
    double RUDY_;
    double lNetRUDY_;
    double gNetRUDY_;
    double sNetRUDY_;



    Graph* graph_;

    std::vector<odb::dbInst*> insts_;
    std::vector<Marker*> markers_;
    std::vector<RSMT*> rsmts_;


  public:
    Gcell(int col, int row);
    
    bgBox getQueryBox();
    odb::Rect getBBox(){ return bbox_; }
    //void updateResourceRSMT(odb::Rect seg);
    //void updateResourceGR(odb::Rect seg);
    //void addInst(odb::dbInst* inst);
    void createTree();
    std::vector<odb::Rect> getSegments(); // available after createTree()
    std::vector<odb::dbInst*> getInsts();
    std::set<odb::dbInst*> getInstSet();

    //
    void extractPlaceFeature(BoxRtree<odb::dbInst*> *rtree);
    void extractPlaceFeature(SegRtree<RSMT*> *rtree);
    void extractRouteFeature(SegRtree<odb::dbNet*> *rtree);

    // 
    void annotateLabel(BoxRtree<Marker*> &rtree);

    //void extractFeature((void*)rtree, ModelType type);

    // helper
    int getCol();
    int getRow();
    int getNumInsts();
    int getNumTerms();
    int getNumNets();
    int getNumGNets();
    int getNumLNets();
    uint getArea();
    uint getCellArea();
    uint getPinArea();
    int getNumMarkers();

    void getNumMarkers(int &lnet, int &gnet, int &inst);

    double getAvgTerms();

    double getRUDY();
    double getLNetRUDY();
    double getSNetRUDY();
    double getGNetRUDY();
    double getPinUtil();
    double getCellUtil();
    double getLNetUtil(ModelType type);
    double getGNetUtil(ModelType type);
    double getWireUtil(ModelType type);
    double getChanUtil(ModelType type);
    double getChanUtilV(ModelType type);
    double getChanUtilH(ModelType type);
    double getChanUtil(Orient orient, ModelType = ModelType::TREE);
    double getBufferUtil();
    double getTNS();
    double getClkRatio();


    int getTrackDemand(Orient orient, ModelType type = ModelType::TREE);
    int getTrackSupply(Orient orient, ModelType type = ModelType::TREE);
    int getWireCapacity(ModelType type = ModelType::TREE);
    void setTrackSupply(int tSup);
    void setWireCapacity(int wCap);
    void setNumLayers(int nLyr);
    void setBoundary(odb::Rect rect);
    void print();

    // initGraph() in ClipGraphExtractor
    void setGraph(Graph* graph);
    void saveGraph(std::string dirPath, std::string fileName);
};

class Marker {
  public:
    
    enum Tag { BoC, PoC, RWoN, NONE };  
    enum Category { L2L, L2I, L2G, G2I, G2G, I2I, ERR, SELF };

    Marker();
    void print();
    void setType(std::string type);
    void setRule(std::string rule);
    void setBoundary(odb::Rect rect);
    void setFromTag(Tag tag);
    void setToTag(Tag tag);
    void setFromNet(RSMT* rsmt);
    void setToNet(RSMT* rsmt);
    void setFromInst(odb::dbInst* inst);
    void setToInst(odb::dbInst* inst);


    bool isFromNet();
    bool isToNet();
    bool isFromInst();
    bool isToInst();


    RSMT* getFromNet();
    RSMT* getToNet();
    odb::dbInst* getToInst();

    bgBox getQueryBox();
    odb::Rect getBBox();
    odb::Point getCentor();

    Category getCategory();
    Tag getFromTag();
    Tag getToTag();

    std::string getType();
    std::string getRule();


  private:
    std::string type_;
    std::string rule_;
    Tag fromTag_;
    Tag toTag_;
    
    RSMT* fromNet_;
    RSMT* toNet_;
    odb::dbInst* fromInst_;
    odb::dbInst* toInst_;


    odb::Rect bbox_;

};


class RSMT {
  private:
    odb::dbNet* net_;
    odb::Rect bbox_;
    std::vector<odb::Point> terminals_;
    Flute::Tree rsmt_;
    int width_;


    std::vector<Gcell*> rsmtOverlaps_;
    std::vector<Gcell*> bboxOverlaps_;
    std::vector<Marker*> markers_;
   
    
  public:
    RSMT();
    RSMT(odb::dbNet* net);
    odb::dbNet* getNet() { return net_; }
    std::vector<odb::Rect> getSegments();

    bgBox getQueryBox();
    odb::Rect getBBox();
    void addTerminal(int x, int y);
    void searchOverlaps(BoxRtree<Gcell*> *tree);
    void setWireWidth(int width);
    

    bool isLocalNet();
    bool isGlobalNet();
    bool hasDRV();
    void createTree();
    int getNumTerms();
    int getWireLengthRSMT();
    int getWireLengthHPWL();
    double getWireUniformUtil();

};

class Grid {

    odb::Rect bbox_;
    int numCols_, numRows_;
    int gcellWidth_, gcellHeight_;
    int wireCapacity_;
    int trackSupply_;
    int numLayers_;
    int minWidth_;


    odb::dbDatabase* db_;
    std::vector<Gcell*> gcells_;
    std::vector<RSMT*> rsmts_;
    // After read drc.rpt
    std::vector<Marker*> markers_;
   
    std::unordered_map<odb::dbNet*, RSMT*> net2rsmt_;

  public:
    odb::dbDatabase* getDb() { return db_; }
    odb::Rect getBoundary();
    std::vector<Gcell*> getGcells();

    Gcell* createGcell(int col, int row, int width, int height); //int x1, int y1, int x2, int y2);
    RSMT* createRSMT(odb::dbNet* net);
    Marker* createMarker(int x1, int y1, int x2, int y2);
    
    
    RSMT* getRSMT(odb::dbNet* net);

    int getNumRows() { return numCols_; }
    int getNumCols() { return numRows_; }
    void init();
    void setWireMinWidth(int width);
    void setDb(odb::dbDatabase* db);
    void setBoundary(odb::Rect rect);
    void setGcellWidth(int width);
    void setGcellHeight(int height);
    void setWireCapacity(int wCap);
    void setTrackSupply(int tSup);
    void setNumLayers(int nLyr);
    void saveGridImages(std::string dirPath);

    
    


    
    //
    double getMaxRUDY();
    double getMaxCellUtil();

    void reportDRC();

};


};

#endif
