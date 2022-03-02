#include "grid.h"
#include "CImg.h"
#include <iostream>
#include <tinycolormap.hpp>


namespace feature_extractor {
using namespace std;
using namespace cimg_library;
using namespace odb;
using namespace tinycolormap;
typedef cimg_library::CImg< unsigned char > CImgObj;
static const unsigned char yellow[] = {255, 255, 0}, white[] = {255, 255, 255},
                           green[] = {0, 255, 0}, blue[] = {120, 200, 255},
                           darkblue[] = {69, 66, 244},
                           purple[] = {255, 100, 255}, black[] = {0, 0, 0},
                           red[] = {255, 0, 0};


void Grid::saveMapImages(string dirPath) {

    int dbUnitMicron = getDb()->getChip()->getBlock()->getDbUnitsPerMicron();

    // img scaling factor
    float sf = 0.01;
    int imgWidth = (int) ( sf * getBoundary().dx() );
    int imgHeight = (int) ( sf * getBoundary().dy() );
    int imgDepth = 1;   
    int imgChannel = 3; // RGB
    int imgSpectrum = 255; // RGB
    int xOrigin = getBoundary().xMin();
    int yOrigin = getBoundary().yMin();
    float opacity = 1.0;
    

    cout << "imageSize (" << imgWidth << " " << imgHeight << ")" << endl;
    // RUDY Map
    //
    CImgObj img1(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum);
    for(Gcell* gcell : gcells_) {

        gcell->getBBox().print();
        cout << "   - cell density  : " << gcell->getCellDensity() << endl;
        cout << "   - pin density   : " << gcell->getPinDensity() << endl;
        cout << "   - RUDY          : " << gcell->getRUDY() << endl;


        double value = gcell->getRUDY(); // + 0.5;
        value = min(1.0, value);
        
        int lev = min((int)(100 * gcell->getRUDY()) % 10, 8);
        //const unsigned char *color = RLEV[lev];
        int x1 = gcell->getBBox().xMin();
        int y1 = gcell->getBBox().yMin();
        int x2 = gcell->getBBox().xMax();
        int y2 = gcell->getBBox().yMax();
        x1 = (int) ( sf * (x1 - xOrigin) );
        x2 = (int) ( sf * (x2 - xOrigin) );
        y1 = (int) ( sf * (y1 - yOrigin) );
        y2 = (int) ( sf * (y2 - yOrigin) );
        
        const Color tinyColor = GetColor(value, ColormapType::Heat);
        unsigned char color[] = {tinyColor.ri(), tinyColor.gi(), tinyColor.bi()};
        img1.draw_rectangle(x1,y1,x2,y2,color,opacity);
    }
    string imgName = "RUDY";
    string imgPath = dirPath + "/RUDY.jpg";
    img1.draw_text(50,50,imgName.c_str(), black, NULL, 1, 30);
    img1.save_jpeg(imgPath.c_str(), 200);
}




};
//static const unsigned char  RLEV[9][3] = {{255, 245, 235}, {254, 230, 206},{253, 208, 162},
//                                        {253, 174, 107}, {253, 141, 60}, {241, 105, 19},
//                                        {217, 72, 1},  {166, 54, 3},{127, 39, 4}},
//                            BLEV[9][3] = {{252, 251, 253}, {239, 237, 245}, {218, 218, 235},
//                                        {188, 189, 220}, {158, 154, 200}, {128, 125, 186},
//                                        {106, 81, 163}, {84, 39, 143}, {63, 0, 125}};
//

