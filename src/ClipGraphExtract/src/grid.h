#include "opendb/geom.h"



namespace Flute {
    class Tree;
};

namespace odb {
    class dbDatabase;
};



struct MyPoint {
    int x, y;
};

struct MyBox {
    MyPoint p1, p2;
};

struct MyLine {
    MyPoint p1, p2;
};


enum ObjectType {
    NET, INSTANCE, BLOCKAGE, SNET
};

class DrcMarker {
    std::string type_;
    std::string rule_;
    ObjectType objType1_, objType2_;

    void* obj1_;
    void* obj2_;
    odb::Rect bbox_;
};


class RtTree {
  private:
    odb::dbNet* net_;
    odb::Rect bbox_;
    std::vector<odb::Point> terminals_;
    Flute::Tree rsmt_;

    std::vector<Gcell*> rsmtOverlaps_;
    std::vector<Gcell*> bboxOverlaps_;
    std::vector<DrcMarker*> markers_;
    
  public:
    RtTree(odb::dbNet* net) : net_(net) {}

    std::vector<odb::Rect> decomposeRSMT();

    bgBox getQueryBox();
    odb::Rect getBBox();
    void addTerminal(int x, int y);

    bool isLocalNet();
    bool isGlobalNet();
    bool hasDRV();
    void createRSMT();
    int getNetDegree() { return terminals_.size(); }
    int getWireLengthRSMT();
    int getWireLengthHPWL();
    double getWireUniformDensity();


    
};

enum Orient {
    LEFT=0,
    RIGHT=1,
    TOP=2,
    BOTTOM=3
};

struct ResourceModel {

    int numSupply[4] = {0};
    int numDemand[4] = {0};
    int wireCapacity;
    int wireLength;
    
    double getWireUtilization();
    double getNumSupply(Orient type);
    double getNumDemand(Orient type);
    double getWireCapaicty();
    double getWireLength();
};


class Gcell {
  private:   
    // Gcell features
    MyBox bbox_;
    odb::Rect rect_;

    ResourceModel rmGR; // using GR results
    ResourceModel rmDR; // using DR results
    ResourceModel rmRSMT; // using PLACE results

    // using placement, RSMT results
    int numInstances_;
    int numTerminals_;
    int numLocalNets_;
    int numGlobalNets_;

    double utilization_;
    double RUDY_;

    std::vector<dbInst*> insts_;

  public:
    odb::Rect getBBox(){ return rect_; }
    void updateResourceModelRSMT(odb::Rect seg);
    void updateResourceModelGR(odb::Rect seg);




};

class Grid {

    odb::Rect rect_;
    int numCols_, numRows_;
    int gridWidth_, gridHeight_;

    odb::dbDatabase* db_;
    std::vector<Gcell> gcells_;
    std::vector<RtTree> rtTrees_;
    // After read drc.rpt
    std::vector<DrcMarker> markers_;
    

};




