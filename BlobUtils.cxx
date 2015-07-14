#include "BlobUtils.h"
#include "TCutG.h"
#include <iostream>

// Returns 1 if the seg2 overlaps in TanB, with a higher mass
// Returns -1 if the seg2 overlaps in TanB, with a lower mass
int Segment::overlaps( const Segment &seg2) const {
  if ( (seg2.tb_end() < tb_start()) || (seg2.tb_start()>tb_end())) return 0;
  if (seg2.mass()>mass()) return 1;
  else return -1;
}


void Segment::Print() const {
  std::cout << "Segment Mass=" << mass() << " TanB start="<< tb_start() << " end=" << tb_end() << std::endl;
}

// blob with mass deltaM less than the segment mass?
bool Blob::overlaps( const Segment &seg ) const {
  Mass lastmass = seg.mass() - m_deltaM;
  SegMap::const_iterator segments = m_segments.find(lastmass);
  if (segments == m_segments.end()) return 0;
  std::cout << " lastmass: " << lastmass << std::endl;
  seg.Print();
  Print();
  for( std::vector<const Segment*>::const_iterator it= (*segments).second.begin(); it<(*segments).second.end(); it++ ) {
    if ((*it)->overlaps(seg)) return 1;
  }
  return 0;
}

void Blob::add( const Segment * seg) {
  if (m_segments.count(seg->mass())==0) m_segments[seg->mass()] = std::vector<const Segment*>();
  m_segments[seg->mass()].push_back(seg);
}

Blob& Blob::operator+=(const Blob& rhs) {
  if (rhs.m_deltaM != m_deltaM)  std::cout << "WARNING - trying to add two blobs with different delaM" << std::endl;
  // For each mass, merge the vectors of segments
  for (SegMap::iterator kv = m_segments.begin(); kv != m_segments.end(); kv++) {
    SegMap::const_iterator rhs_segs = rhs.m_segments.find((*kv).first);
    if (rhs_segs == rhs.m_segments.end()) continue;
    (*kv).second.insert( (*kv).second.end(),  (*rhs_segs).second.begin(), (*rhs_segs).second.end() );
  }
  return *this;
}

void Blob::Print() const {
  std::cout << "=== Blob ==="  << std::endl;
  for (SegMap::const_iterator kv = m_segments.begin(); kv != m_segments.end(); kv++) {
    //std::cout << "Mass " << (*kv).first << " " << (*kv).second.size() << " segments:" << std::endl; 
   for( std::vector<const Segment*>::const_iterator it= (*kv).second.begin(); it<(*kv).second.end(); it++ ) (*it)->Print();
  }
}



std::vector<TCutG*> Blob::TCutGs(TString name, float min_tb) const {
  // First find outline. 
  // Find top left point
  BlobStep start(1, (*m_segments.begin()).second.back());
  BlobStep curr=start;
  BlobStep last=start;
  int iP=0;
  TCutG* cutg = new TCutG(name+"_excl", 0);
  do {
    float last_tb = (last.pos().second>min_tb) ? last.pos().second : min_tb;
    float curr_tb = (curr.pos().second>min_tb) ? curr.pos().second : min_tb;
    cutg->SetPoint(iP++, last.pos().first, last_tb);
    if (curr.dir()!=last.dir()) {
      // Change of direction - add extra point
      cutg->SetPoint(iP++, last.pos().first+last.dir()*m_deltaM*0.25, (last_tb+curr_tb)/2);
    }
    last=curr;
    // Find the next point
    curr = this->jump(last);
    //float last_tb = (last.pos().second>min_tb) ? last.pos().second : min_tb;
    //float curr_tb = (curr.pos().second>min_tb) ? curr.pos().second : min_tb;
    // Add to cutg
    /*
    if (last.pos().second>=min_tb) {
      cutg->SetPoint(iP++, last.pos().first, last.pos().second);
      if (curr.dir()!=last.dir()) {
        // Change of direction - add extra point
        cutg->SetPoint(iP++, last.pos().first+last.dir()*m_deltaM*0.5, (last.pos().second+curr.pos().second)/2);
      }
    }
    */
  } while (curr!=start);
  // stop if we've got back to the start
  // Add first point again, and a helper point
  if (last.pos().second<=min_tb) 
  {
    cutg->SetPoint(iP++, curr.pos().first-2*m_deltaM, min_tb);
    cutg->SetPoint(iP++, curr.pos().first-m_deltaM, min_tb);
    cutg->SetPoint(iP++, curr.pos().first, min_tb);
    cutg->SetPoint(iP++, curr.pos().first, min_tb*1.001);
  }
  cutg->SetPoint(iP++, curr.pos().first, curr.pos().second);
  for (int i=0; i<cutg->GetN(); i++) {
    double x,y;
    cutg->GetPoint(i, x, y);
    //std::cout  << i << " : " << x << " : " << y << std::endl;
  }
  std::vector<TCutG*> v;
  v.push_back(cutg);
  return v;
}

Blob::BlobStep Blob::jump(const BlobStep& curr) const {
  int dir = curr.dir();
  float mass = curr.pos().first;
  float next_mass = mass + dir*m_deltaM;
  const Segment* seg = curr.seg(); 
  // See if there are any segments at next mass
  auto next_segs = m_segments.find(next_mass);
  const Segment* next_seg=0;
  if (next_segs!=m_segments.end()) {
    // There are segments - now find the last/first one which overlaps 
    // Here assuming segments are ordered.
    // this should always be the case
    for( auto tseg : (*next_segs).second ) {
      if (!tseg->overlaps(*seg)) continue;
      next_seg = tseg;
      // if going backwards, stop on first seg (will be lowest in tanb)
      if (dir<0) break;
      //otherwsie continue on till last (will be highest in tanb)
    }
  }
  if (next_seg) {
    // Check there isn't another segment at this mass that
    // overlaps with the next_seg
    auto segs = (*m_segments.find(mass)).second;
    auto it = std::find(segs.begin(), segs.end(), seg);
    if (dir>0) it++;
    else it--;
    if (it!=segs.end() and it>=segs.begin() and  (*it)->overlaps(*next_seg)) {
      // found another segment that overlaps next_seg
      // switch direction
      return BlobStep(dir*-1, *it);
    } else {
      return BlobStep(dir, next_seg);
    }
  } else {
    return BlobStep(dir*-1, seg);
  }
}

Blob::BlobStep::BlobStep( int dir, const Segment* seg ) {
  m_dir = dir;
  float tb = (dir>0) ? seg->tb_end() : seg->tb_start();
  m_pos = Point(seg->mass(), tb);
  m_seg = seg;
};

// Update the blobs with the limit info in std::vector<TLimit_Info>
// the deltaM should give the difference in mass between this layer and the last
// minTb is the dummy tanB value to use for start segments that are out of range
void update_blobs( std::vector<TLimit_Info> v_info, std::vector<Blob> &v_blobs, float deltaM, float minTb ) {
  std::vector<const Segment*> segments;
  if (v_info.size()==0) return;
  // If start on a downcrossing, add a pseudo-upcorssing at start
  if (v_info[0].upcross==0) {
    TLimit_Info info;
    info.mA = v_info[0].mA;
    info.tanb = minTb;
    info.upcross = 1;
    info.extrap = 1;
    v_info.insert( v_info.begin(), info);
  }
  // Sanity checks
  if ( v_info.size()%2 !=0) {
    std::cout << "WARNING - odd number of TLimit_Info - skipping" << std::endl;
    for (auto info : v_info) std::cout << " mA " << info.mA << " tanb " << info.tanb 
      << " upcross " << info.upcross << " extrap " << info.extrap << std::endl;
    return;
  }
  // Loop over pairs of TLimit_Info to form Segments
  int iSeg=0;
  for (int i=0; (i+1)<v_info.size(); i+=2) {
    // Skip Segment if end is below minimum tan(beta)
    if (v_info[i+1].tanb < minTb) continue;
    Segment *s = new Segment(v_info[i].mA, v_info[i].tanb, v_info[i+1].tanb);
    iSeg++;
    // See if we can odd blob to an existing segment
    bool found_blob=false;
    for (std::vector<Blob>::iterator blob_it = v_blobs.begin(); blob_it != v_blobs.end(); blob_it++) {
      if ((*blob_it).overlaps(*s)) {
        (*blob_it).add(s);
        found_blob=true;
        break;
      }
    }
    if (!found_blob) {
      // Otherwise make a new blob
      v_blobs.push_back( Blob(deltaM) );
      std::cout << s << std::endl;
      v_blobs.back().add(s);
    }
  }
  std::cout << "found " << iSeg << " segments. There are  " << v_blobs.size() << " blobs " << std::endl;
}
