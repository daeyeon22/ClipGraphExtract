#ifndef MAKE_GRAPH_EXTRACTOR
#define MAKE_GRAPH_EXTRACTOR

namespace ClipGraphExtract {
class ClipGraphExtractor;
}

namespace ord {

class OpenRoad;

ClipGraphExtract::ClipGraphExtractor*
makeClipGraphExtractor();

void
initClipGraphExtractor(OpenRoad *openroad);

void
deleteClipGraphExtractor(ClipGraphExtract::ClipGraphExtractor *graphext);

}

#endif
