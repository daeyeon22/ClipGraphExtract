#include <tcl.h>
#include "sta/StaMain.hh"
#include "openroad/OpenRoad.hh"
#include "clip_graph_ext/MakeClipGraphExtractor.h"
#include "clip_graph_ext/clipGraphExtractor.h"


namespace sta {
extern const char *graph_extractor_tcl_inits[];
}

extern "C" {
extern int Clipgraphextractor_Init(Tcl_Interp* interp);
}

namespace ord {

ClipGraphExtract::ClipGraphExtractor * 
makeClipGraphExtractor() 
{
  return new ClipGraphExtract::ClipGraphExtractor; 
}

void 
initClipGraphExtractor(OpenRoad *openroad) 
{
  Tcl_Interp* tcl_interp = openroad->tclInterp();
  Clipgraphextractor_Init(tcl_interp);
  sta::evalTclInit(tcl_interp, sta::graph_extractor_tcl_inits);
  openroad->getClipGraphExtractor()->setDb(openroad->getDb());
  openroad->getClipGraphExtractor()->setSta(openroad->getSta());
}

void
deleteClipGraphExtractor(ClipGraphExtract::ClipGraphExtractor *graphext)
{
  delete graphext;
}


}
