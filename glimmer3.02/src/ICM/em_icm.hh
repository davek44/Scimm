//  A. L. Delcher
//
//  File:  icm.hh
//
//  Last Modified:  3 May 2004
//
//  Declarations for the Interpolated Context Model (ICM) used by Glimmer



#ifndef _EM_ICM_H_INCLUDED
#define _EM_ICM_H_INCLUDED


#include  "delcher.hh"
#include  "gene.hh"
#include "icm.hh"

using namespace std;


// Temporary so it will compile

struct  em_ICM_Training_Node_t
  {
   short int  mut_info_seq;
     // The base that is to be restricted in mut_info_pos
   double  (* count) [ALPHA_SQUARED];
     // count [x] [y] (where y is a pair of letters (e.g., aa, ac, ag, .. tt))
     // is the number of times the first base of y occurs at position x
     // and the last base of y occurs at position (model_len - 1)
  };


struct ICM_Training_Data_t
{
  char * seq;
  double p;
};

class  em_ICM_t
  {
  protected:
   bool  empty;
   int  model_len;
     // Number of bases in window
   int  model_depth;
     // Most levels in model tree, i.e., most number of positions
     // to be dependent on
   int  periodicity;
     // Number of different models, alternating cyclicly
   int  num_nodes;
     // Number of nodes in tree
   ICM_Score_Node_t  * * score;

  public:
   em_ICM_t
       (int m = DEFAULT_MODEL_LEN,
        int d = DEFAULT_MODEL_DEPTH,
        int p = DEFAULT_PERIODICITY);
   ~ em_ICM_t
       ();

   int  Get_Model_Depth
       (void)
     { return  model_depth; }
   int  Get_Model_Len
       (void)
     { return  model_len; }
   int  Get_Periodicity
       (void)
     { return  periodicity; }

   void  Build_Indep_WO_Stops
       (double gc_frac, const vector <const char *> & stop_codon);
   void  Build_Reverse_Codon_WO_Stops
    (double codon_prob [64], const vector <const char *> & stop_codon);
   void  Cumulative_Score
       (const string & s, vector <double> & score, int frame)  const;
   void  Cumulative_Score_String
       (char * string, int len, int frame, double * score);
   void  Display
       (FILE * fp);
   void  Full_Window_Distrib
       (char * string, int frame, float * dist);
   double  Full_Window_Prob
       (const char * string, int frame)  const;
   void  Input
       (FILE * fp);
   void  Output
       (FILE * fp, bool binary_form);
   void  Output_Node
       (FILE * fp, ICM_Score_Node_t * node, int id, int frame, bool binary_form);
   double  Partial_Window_Prob
       (int predict_pos, const char * string, int frame)  const;
   void  Read
       (char * path);
   double  Score_String
       (const char * string, int len, int frame)  const;
   void  Set_Label_String
       (char * label, int id, int frame);
   void  Write_Header
       (FILE * fp, bool binary_form);
  };


class  em_ICM_Training_t  :  public ICM_t
  {
  private:
   em_ICM_Training_Node_t  * * train;
   double  (* count_memory) [ALPHA_SQUARED];
     // Holds dynamically allocated block for all counts to avoid
     // loss at tail from separate callocs.  Saved here to allow
     // freeing it in the destructor.

   void  Complete_Tree
       (const vector <ICM_Training_Data_t> & data);
   void  Count_Char_Pairs_Restricted
       (const ICM_Training_Data_t string, int level);
   em_ICM_Training_Node_t *  Get_Training_Node
       (const char * w, int frame, int level);
   void  Interpolate_Probs
       (int frame, int sub, double ct []);
   void  Take_Logs
       (void);

  public:
   em_ICM_Training_t
       (int m = DEFAULT_MODEL_LEN,
        int d = DEFAULT_MODEL_DEPTH,
        int p = DEFAULT_PERIODICITY);
   ~ em_ICM_Training_t
       ();

   void  Train_Model
       (const vector <ICM_Training_Data_t> & data);
  };


void  em_Count_Char_Pairs
    (double ct [] [ALPHA_SQUARED], char * string, double p, int w, int period);
void  em_Count_Single_Chars
    (double ct [ALPHABET_SIZE], char * string, double p, int w, int period);
double  em_Get_Mutual_Info
    (double ct [], int n, double sum);

#endif

