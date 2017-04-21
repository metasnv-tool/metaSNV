//==========================================================================
//============================== CODONS ====================================
const std::pair<std::string, char> codons[] = {
  std::make_pair("TAA",'X'),std::make_pair("TGA",'X'),std::make_pair("TAG",'X'),//STOP
  std::make_pair("GCT",'A'),std::make_pair("GCC",'A'),std::make_pair("GCA",'A'),std::make_pair("GCG",'A'),//Ala
  std::make_pair("CGT",'R'),std::make_pair("CGC",'R'),std::make_pair("CGA",'R'),std::make_pair("CGG",'R'),std::make_pair("AGA",'R'),std::make_pair("AGG",'R'),//Arg
  std::make_pair("AAT",'N'),std::make_pair("AAC",'N'),//Asn
  std::make_pair("GAT",'D'),std::make_pair("GAC",'D'),//Asp
  std::make_pair("TGT",'C'),std::make_pair("TGC",'C'),//Cys
  std::make_pair("CAA",'Q'),std::make_pair("CAG",'Q'),//Gln
  std::make_pair("GAA",'E'),std::make_pair("GAG",'E'),//Glu
  std::make_pair("GGT",'G'),std::make_pair("GGC",'G'),std::make_pair("GGA",'G'),std::make_pair("GGG",'G'),//Gly
  std::make_pair("CAT",'H'),std::make_pair("CAC",'H'),//His
  std::make_pair("ATT",'I'),std::make_pair("ATC",'I'),std::make_pair("ATA",'I'),//Ile
  std::make_pair("TTA",'L'),std::make_pair("TTG",'L'),std::make_pair("CTT",'L'),std::make_pair("CTC",'L'),std::make_pair("CTA",'L'),std::make_pair("CTG",'L'),//Leu
  std::make_pair("AAA",'K'),std::make_pair("AAG",'K'),//Lys 
  std::make_pair("ATG",'M'),//Met and START     
  std::make_pair("TTT",'F'),std::make_pair("TTC",'F'),//Phe
  std::make_pair("CCT",'P'),std::make_pair("CCC",'P'),std::make_pair("CCA",'P'),std::make_pair("CCG",'P'),//Pro
  std::make_pair("TCT",'S'),std::make_pair("TCC",'S'),std::make_pair("TCA",'S'),std::make_pair("TCG",'S'),std::make_pair("AGT",'S'),std::make_pair("AGC",'S'),//Ser
  std::make_pair("ACT",'T'),std::make_pair("ACC",'T'),std::make_pair("ACA",'T'),std::make_pair("ACG",'T'),//Thr
  std::make_pair("TGG",'W'),//Trp
  std::make_pair("TAT",'Y'),std::make_pair("TAC",'Y'),//Tyr
  std::make_pair("GTA",'V'),std::make_pair("GTG",'V'),std::make_pair("GTT",'V'),std::make_pair("GTC",'V'),//Val
};

//==========================================================================
//Define the 3 bit to base mapping and reverse
char intToBase[5] = {'A','T','C','G','N'};

const std::pair<char,short> baseToShort[] = {
  std::make_pair('A',0),
  std::make_pair('T',1),
  std::make_pair('C',2),
  std::make_pair('G',3),
  std::make_pair('N',4)
};
//==========================================================================
/**
   @brief Optimal genome representation
*/
class Genome{

public:
  /**
     @brief Standard string constructor
   */
  Genome(std::string seq) {//Get this giant string into a manageble form!
    //Initialize maps
    std::map<char,short> map_baseToShort(baseToShort, baseToShort + sizeof baseToShort / sizeof baseToShort[0]);

    length = seq.length();
    //Create internal sequence structure. Each int holds 10 letters.
    sequence.resize((length/10) + 1);
    unsigned int repr = 0;
    int pos = 0;
    long insertPosition = 0;
    for (int i=0; i<length; ++i,++pos) {
        if (pos == 10) {
            sequence[insertPosition] = repr;
            ++insertPosition;
            repr = 0;
            pos = 0;
        }
        //Note that this is going to represent the furthest away base on the most significant bits.
        //I can only hope that makes sense. Should be kept coherent throughout the class.
        repr |= (map_baseToShort[seq[i]] << 3*pos);
    }
    if (pos!=0) {//Add the end too
      sequence[insertPosition] = repr;
    }
    //And done!
  }

  /**
     @brief Return sequence between coordinates, including end
   */

  std::string getSequence(long start, long end) {
    //Start with some sanity checks
    if (end < start) {
      return "";
    } else if (end > length) {
      return "";
    }
    std::string res = "";
    //Now, actually get what the user wants
    for (long i=start; i<=end; ++i) {
      res += intToBase[(sequence[i/10] >> 3*(i%10)) & 7];
    }
    return res;
  };

  long getLength() {
    return length;
  };

private:

  long length;
  std::vector<unsigned int> sequence;
};

//==========================================================================
/**
   @brief Gene definition structure
*/
struct Gene{

  std::string name;
  long start,end;
  char strand;

  Gene(long s, long e, std::string n, char st)
    :name(n)
    ,start(s)
    ,end(e)
    ,strand(st) { }

  ~Gene() { }

};

class GeneDef {
 public:
  GeneDef() {
    m_Genes.clear();
  };

  GeneDef(Gene& g) {
    m_Genes.clear();
    m_Genes.push_back(g);
  };

  bool hasGene() {
    return m_Genes.size() != 0;
  };

  const Gene& getGene() const {
    return m_Genes.front();
  };

  GeneDef& operator+=(const GeneDef &g) {
    m_Genes.push_back(g.getGene());
    return *this;
  };

  bool operator==(const GeneDef &g) const {
    return false;
  };

private:
  std::list<Gene> m_Genes;
};
