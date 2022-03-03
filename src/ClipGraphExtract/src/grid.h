#ifndef __GCELLGRID__
#define __GCELLGRID__
#include <vector>
#include <string>
#include <boost/geometry.hpp>
#include "flute.h"
#include "opendb/db.h"
#include "opendb/geom.h"

//#include "bgTypedef.h"
//#include "opendb/geom.h"

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

//namespace Flute {
//    class Tree;
//};

//namespace odb {
//    class dbDatabase;
//    class dbNet;
//    class dbInst;
//    class Rect;
//    class Point;
//};

namespace feature_extractor {

class DrcMarker;
class Gcell;
class RSMT;


enum Orient {
    LEFT=0,
    RIGHT=1,
    TOP=2,
    BOTTOM=3
};

enum ModelType {
    DR,EGR,PL
};

class Resource {
  private:
    uint trackSupply_[4] = {0};
    uint trackDemand_[4] = {0};
    uint wireCapacity_;
    uint wireLength_;

  public:
    
    Resource() : wireCapacity_(0), wireLength_(0) {}

    double getChannelDensity(Orient type) { return 1.0 * trackDemand_[type] / trackSupply_[type]; }
    double getWireDensity() { return 1.0 * wireLength_ / wireCapacity_; }
    uint getTrackSupply(Orient type) { return trackSupply_[type]; }
    uint getTrackDemand(Orient type) { return trackDemand_[type]; }
    uint getWireCapacity() { return wireCapacity_; }
    uint getWireLength() { return wireLength_; }
    void incrTrackDemand(Orient type) { trackDemand_[type]++; }
    void addTrackDemand(Orient type, uint dem) { trackDemand_[type] += dem; }
    void setTrackSupply(uint tSup) {  for(int i=0; i < 4; i++) trackSupply_[i] = tSup; }
    void setWireCapacity(uint wCap) {  wireCapacity_ = wCap;  }
    void setWireLength(uint wl) { wireLength_ = wl; }
    void addWireLength(uint wl) { wireLength_ += wl; }
};


class Gcell {
  private:   
    // Gcell features
    odb::Rect bbox_;

    Resource rmEGR_; // using GR results
    Resource rmDR_; // using DR results
    Resource rmPL_; // using PLACE results

    // using placement, RSMT results
    uint numInstances_;
    uint numTerminals_;
    uint numLocalNets_;
    uint numGlobalNets_;

    uint totalCellArea_;
    uint totalPinArea_;

    double cellDensity_;
    double pinDensity_;
    double RUDY_;

    std::vector<odb::dbInst*> insts_;


  public:
    Gcell();
    
    bgBox getQueryBox();
    odb::Rect getBBox(){ return bbox_; }
    //void updateResourceRSMT(odb::Rect seg);
    //void updateResourceGR(odb::Rect seg);
    //void addInst(odb::dbInst* inst);
    void createTree();
    std::vector<odb::Rect> getSegments(); // available after createTree()

    //
    void extractFeaturePL(BoxRtree<odb::dbInst*> &rtree);
    void extractFeatureEGR(SegRtree<odb::dbNet*> &rtree);
    void extractFeatureRSMT(SegRtree<RSMT*> &rtree);

    //void extractFeature((void*)rtree, ModelType type);

    // helper
    uint getNumInstances();
    uint getNumTerminals();
    uint getNumGlobalNets();
    uint getNumLocalNets();
    uint getArea();
    uint getCellArea();
    uint getPinArea();
    double getRUDY();
    double getPinDensity();
    double getCellDensity();
    double getWireDensity(ModelType type);
    double getChannelDensity(ModelType type);
    double getChannelDensity(Orient orient, ModelType = ModelType::PL);
    uint getTrackDemand(Orient orient, ModelType type = ModelType::PL);
    uint getTrackSupply(Orient orient, ModelType type = ModelType::PL);
    uint getWireCapacity(ModelType type = ModelType::PL);
    void setTrackSupply(int tSup);
    void setWireCapacity(int wCap);
    void setBoundary(odb::Rect rect);
    void print();

};


enum ObjectType {
    NET, INSTANCE, BLOCKAGE, SNET
};

enum DrvTag {
    LOCAL2LOCAL,
    LOCAL2GLOBAL,
    GLOBAL2GLOBAL,
    PIN2LOCAL,
    PIN2GLOBAL
};

class DrcMarker {
  private:
    std::string type_;
    std::string rule_;
    ObjectType objType1_, objType2_;

    void* obj1_;
    void* obj2_;
    odb::Rect bbox_;

  public:
    void setType(std::string type);
    void setDesignRule(std::string rule);
    void setBoundary(odb::Rect rect);
    void setObject1(odb::dbNet* net);
    void setObject2(odb::dbNet* net);
    void setObject1(odb::dbInst* inst);
    void setObject2(odb::dbInst* inst);

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
    std::vector<DrcMarker*> markers_;
   
    
  public:
    RSMT();
    RSMT(odb::dbNet* net);
    odb::dbNet* getNet() { return net_; }
    std::vector<odb::Rect> getSegments();

    bgBox getQueryBox();
    odb::Rect getBBox();
    void addTerminal(int x, int y);
    void searchOverlaps(BoxRtree<Gcell*> &tree);
    void setWireWidth(int width);
    

    bool isLocalNet();
    bool isGlobalNet();
    bool hasDRV();
    void createTree();
    int getNetDegree(); 
    uint getWireLengthRSMT();
    uint getWireLengthHPWL();
    double getWireUniformDensity();

};

class Grid {

    odb::Rect bbox_;
    int numCols_, numRows_;
    int gcellWidth_, gcellHeight_;
    int wireCapacity_;
    int trackSupply_;
    int minWidth_;


    odb::dbDatabase* db_;
    std::vector<Gcell*> gcells_;
    std::vector<RSMT*> rsmts_;
    // After read drc.rpt
    std::vector<DrcMarker*> markers_;
    

  public:
    odb::dbDatabase* getDb() { return db_; }
    odb::Rect getBoundary();
    std::vector<Gcell*> getGcells();

    Gcell* createGcell(int x1, int y1, int x2, int y2);
    RSMT* createRSMT(odb::dbNet* net);
    void init();
    void setWireMinWidth(int width);
    void setDb(odb::dbDatabase* db);
    void setBoundary(odb::Rect rect);
    void setGcellWidth(int width);
    void setGcellHeight(int height);
    void setWireCapacity(int wCap);
    void setTrackSupply(int tSup);


    void saveMapImages(std::string dirPath);
};



};

#endif
