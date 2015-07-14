#ifndef BlobUtils_h
#define BlobUtils_h

#include <utility>
#include <vector>
#include <map>
#include "LimitUtilities.h"

typedef std::pair<float, float> Point;
typedef float Mass;
typedef float TanB;

class TCutG;

class Segment {

  private:
    Mass m_mass;
    TanB m_tb_start;
    TanB m_tb_end;

  public:
    Segment( Mass mass, TanB tb_start, TanB tb_end ) {
      m_mass = mass; m_tb_start=tb_start; m_tb_end=tb_end;
    };

    // Returns 1 if the seg2 overlaps in TanB, with a higher mass
    // Returns -1 if the seg2 overlaps in TanB, with a lower mass
    int overlaps( const Segment &seg2) const;

    Mass mass() const { return m_mass; };
    TanB tb_start() const { return m_tb_start; };
    TanB tb_end() const { return m_tb_end; };

    void Print() const;
};

class Blob {

  typedef std::map<Mass, std::vector<const Segment*> > SegMap;

  private:
    SegMap m_segments;
    Mass m_deltaM;

    class BlobStep {
      private:
        int m_dir;
        Point m_pos;
        const Segment* m_seg;
      public:
        BlobStep( int dir, const Segment* seg );
        Point pos() const { return m_pos; };
        int dir() const { return m_dir; };
        const Segment* seg() const { return m_seg; };
        bool operator==( const BlobStep& other ) const {
          return (m_pos == other.m_pos);
        };
        bool operator!=( const BlobStep& other ) const {
          return (!(*this == other));
        };
    };
    Blob::BlobStep jump(const BlobStep& curr) const;

  public:
  Blob(Mass deltaM) { m_deltaM = deltaM; };

  // Does this segment overla one of the segments in the
  // blob with mass deltaM less than the segment mass?
  bool overlaps( const Segment &seg ) const;
  void add( const Segment * seg);
  Blob& operator+=(const Blob& rhs);
  friend Blob operator+(Blob lhs, const Blob& rhs) { return lhs += rhs; };
  void Print() const;
  std::vector<TCutG*> TCutGs(TString name, float min_tb) const;
};



// Update the blobs with the limit info in std::vector<TLimit_Info>
// the deltaM should give the difference in mass between this layer and the last
// minTb is the dummy tanB value to use for start segments that are out of range
void update_blobs( std::vector<TLimit_Info> v_info, std::vector<Blob> &v_blobs, float deltaM, float minTb );

#endif
