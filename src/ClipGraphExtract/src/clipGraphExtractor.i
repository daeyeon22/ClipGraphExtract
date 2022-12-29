%module ClipGraphExtractor

%{
#include "openroad/OpenRoad.hh"
#include "clip_graph_ext/clipGraphExtractor.h"

namespace ord {
ClipGraphExtract::ClipGraphExtractor*
getClipGraphExtractor(); 
odb::dbDatabase*
getDb();
sta::dbSta*
getSta();
}



using ord::getClipGraphExtractor;
using ord::getDb;
using ord::getSta;
using ClipGraphExtract::ClipGraphExtractor;

%}

%inline %{

void
parse_drc_report_cmd(const char* file_name) {

    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->parseDrcReport(file_name);

}

void
set_gcell_size_cmd(int numRows) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setGcellSize(numRows);
}

void
set_max_route_layer_cmd(int maxRouteLayer) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setMaxRouteLayer(maxRouteLayer);
}


void
set_graph_extract_save_file_name_cmd(const char* file)
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->setSaveFileName(file);
}

void
set_graph_extract_save_file_prefix_cmd(const char* prefix)
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setSaveFilePrefix(prefix);
}

void
set_drc_report_cmd(const char* fileName) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setDrcReport(fileName);
}



void 
graph_extract_init_cmd()
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->setDb(getDb());
  graphExt->setSta(getSta());
  graphExt->init();
}


void
graph_extract_cmd() 
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->extract();

}

void
save_features_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveFeatures(dirPath);
}


/*
void 
save_inst_features_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveInstFeatures(dirPath);
}

void
save_inst_labels_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveInstLabels(dirPath);
}


void
save_graphs_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveGraphs(dirPath);
}
*/

void
save_labels_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveLabels(dirPath);
}

void
save_grid_images_cmd(const char* imgDir, const char* prefix) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveGridImages(imgDir, prefix);
}


/*
void
graph_extract_cmd(int lx, int ly, int ux, int uy) 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->extract(lx, ly, ux, uy);
}
*/
void
graph_extract_clear_cmd() 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->clear();
}

%}
