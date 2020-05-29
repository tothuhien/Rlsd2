#include "stdarg.h"
#include "utils.h"
#include "date.h"
#include "lsd.h"

using namespace std;

int tree2data(istream& tree,Pr* pr,int & s,Node** &nodes);
int readInputDate(InputOutputStream* io, Pr* pr,Node** &nodes,bool& constraintConsistent);
int readPartitionFile(istream &partFile, Pr* pr);
int tree2dataS(FILE *,Pr*,Node**);
int extrait_outgroup(InputOutputStream *io, Pr* pr);
int getBranchOut(Pr* pr,Node** nodes,list<string> &outgroups,bool &keepBelow,int &r);