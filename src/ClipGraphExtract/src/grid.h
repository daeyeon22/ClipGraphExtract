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

struct ResourceModel {
    int trackSupply[4] = {0};
    int trackDemand[4] = {0};
    int wireCapacity;
    int wireLength;
    
    double getWireDensity();
    double getTrackSupply(Orient type);
    double getTrackDemand(Orient type);
    double getWireCapaicty();
    double getWireLength();

    void setTrackSupply(int tSup) {
        for(int i=0; i < 4; i++) trackSupply[i] = tSup;
    }

    void setWireCapacity(int wCap) { 
        wireCapacity = wCap;
    }
};


class Gcell {
  private:   
    // Gcell features
    odb::Rect bbox_;

    ResourceModel rmEGR; // using GR results
    ResourceModel rmDR; // using DR results
    ResourceModel rmRSMT; // using PLACE results

    // using placement, RSMT results
    int numInstances_;
    int numTerminals_;
    int numLocalNets_;
    int numGlobalNets_;

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
    //void updateResourceModelRSMT(odb::Rect seg);
    //void updateResourceModelGR(odb::Rect seg);
    //void addInst(odb::dbInst* inst);
    void createTree();
    std::vector<odb::Rect> getSegments(); // available after createTree()

    //
    void extractFeaturePL(BoxRtree<odb::dbInst*> &rtree);
    void extractFeatureEGR(SegRtree<odb::dbNet*> &rtree);
    void extractFeatureRSMT(SegRtree<RSMT*> &rtree);

    // helper
    uint getArea();
    uint getCellArea();
    uint getPinArea();
    double getRUDY();
    double getPinDensity();
    double getCellDensity();


    void setTrackSupply(int tSup);
    void setWireCapacity(int wCap);
    void setBoundary(odb::Rect rect);

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
